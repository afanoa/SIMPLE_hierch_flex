!@descr: supported physical-score tail candidates and conditional endpoint diagnostics
module simple_flex_supported_extrema
use, intrinsic :: iso_fortran_env, only: real64
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
private
public :: flex_supported_tail_labels, flex_fission_channels, flex_extrema_cross_statistic
public :: EXTREMA_OK, EXTREMA_INVALID, EXTREMA_POPULATION_LIMITED, EXTREMA_HALF_LIMITED
public :: EXTREMA_CALIBRATION_PENDING, EXTREMA_HEURISTIC, EXTREMA_INDEPENDENT_NOISE, EXTREMA_REFITTED_NULL

integer, parameter :: EXTREMA_OK = 0, EXTREMA_INVALID = 1, EXTREMA_POPULATION_LIMITED = 2
integer, parameter :: EXTREMA_HALF_LIMITED = 3, EXTREMA_CALIBRATION_PENDING = 4
integer, parameter :: EXTREMA_HEURISTIC = 0, EXTREMA_INDEPENDENT_NOISE = 1, EXTREMA_REFITTED_NULL = 2

contains

    !> Candidate only: score must have one physical meaning across halves; counts/ESS do not validate it.
    !! Global sorting breaks exact ties by original row ID, never by EO.
    !! All rows survive as A/middle/B (1/2/3). Native physical acceptance remains external.
    subroutine flex_supported_tail_labels( score, row_ids, eo, ntail, minpop, minhalf, labels, counts, status, message )
        use simple_srch_sort_loc, only: hpsort
        real(real64), intent(in) :: score(:)
        integer, intent(in) :: row_ids(:), eo(:), ntail, minpop, minhalf
        integer, allocatable, intent(out) :: labels(:)
        integer, intent(out) :: counts(3,2), status
        character(len=*), intent(out) :: message
        integer, allocatable :: ord(:), ids(:), candidate(:)
        integer :: n, i, s, h
        counts = 0
        status = EXTREMA_INVALID
        message = 'invalid score or original row/half identity'
        n = size(score)
        if( n == 0 .or. size(row_ids) /= n .or. size(eo) /= n ) return
        if( minpop < 1 .or. minhalf < 1 .or. ntail < 1 ) return
        if( any(.not. ieee_is_finite(score)) ) return
        if( any(row_ids <= 0) .or. any(eo < 0) .or. any(eo > 1) ) return
        allocate(ids(n), source=row_ids)
        call hpsort(ids)
        if( any(ids(2:n) == ids(1:n-1)) ) return
        status = EXTREMA_POPULATION_LIMITED
        message = 'population_limited: every candidate state must meet minpop'
        if( ntail < minpop .or. ntail > n/2 ) return
        if( n - 2*ntail < minpop ) return
        status = EXTREMA_INVALID
        message = 'physical score has no range'
        if( maxval(score) == minval(score) ) return
        allocate(ord(n), candidate(n))
        ord = [(i,i=1,n)]
        call hpsort(ord, score_less)
        candidate = 2
        candidate(ord(1:ntail)) = 1
        candidate(ord(n-ntail+1:n)) = 3
        do s = 1, 3
            do h = 0, 1
                counts(s,h+1) = count(candidate == s .and. eo == h)
            end do
        end do
        status = EXTREMA_HALF_LIMITED
        message = 'half_support_limited: candidate retained as parent'
        if( any(counts < minhalf) ) return
        call move_alloc(candidate, labels)
        status = EXTREMA_OK
        message = 'supported initializer only; physical validation is separate'

    contains

        logical function score_less( i1, i2 )
            integer, intent(in) :: i1, i2
            if( score(i1) == score(i2) )then
                score_less = row_ids(i1) < row_ids(i2)
            else
                score_less = score(i1) < score(i2)
            endif
        end function score_less

    end subroutine flex_supported_tail_labels

    !> Arithmetic only: independent eta with the exact original Gaussian covariance gives independent S/M noise.
    !! Both covariances double. The caller must establish that noise/operator provenance.
    !! Fission of full-y fitted scores removes no basis/pose/mask leakage and creates no particles.
    subroutine flex_fission_channels( observed, independent_noise, selection, measurement, status )
        real(real64), intent(in) :: observed(:,:), independent_noise(:,:)
        real(real64), allocatable, intent(out) :: selection(:,:), measurement(:,:)
        integer, intent(out) :: status
        status = EXTREMA_INVALID
        if( any(shape(observed) /= shape(independent_noise)) .or. size(observed) == 0 ) return
        if( any(.not. ieee_is_finite(observed)) .or. any(.not. ieee_is_finite(independent_noise)) ) return
        allocate(selection, source=observed + independent_noise)
        allocate(measurement, source=observed - independent_noise)
        if( any(.not. ieee_is_finite(selection)) .or. any(.not. ieee_is_finite(measurement)) )then
            deallocate(selection, measurement)
            return
        endif
        status = EXTREMA_OK
    end subroutine flex_fission_channels

    !> Columns are frozen candidates, rows common physical coefficients; null maxima include the same search.
    !! INDEPENDENT_NOISE requires S-only candidates/conditional null; REFITTED_NULL a complete selection/refit null.
    !! Caller provenance is required; map agreement implies no independence or automatic state acceptance.
    subroutine flex_extrema_cross_statistic( delta_even, delta_odd, null_maxima, validation_mode, &
        &jscore, selected, p_search, status )
        real(real64), intent(in) :: delta_even(:,:), delta_odd(:,:), null_maxima(:)
        integer, intent(in) :: validation_mode
        real(real64), allocatable, intent(out) :: jscore(:)
        integer, intent(out) :: selected, status
        real(real64), intent(out) :: p_search
        integer :: c
        selected = 0
        p_search = -1._real64
        status = EXTREMA_INVALID
        if( any(shape(delta_even) /= shape(delta_odd)) .or. size(delta_even) == 0 ) return
        if( any(.not. ieee_is_finite(delta_even)) .or. any(.not. ieee_is_finite(delta_odd)) ) return
        if( validation_mode < EXTREMA_HEURISTIC .or. validation_mode > EXTREMA_REFITTED_NULL ) return
        allocate(jscore(size(delta_even,2)))
        do c = 1, size(jscore)
            jscore(c) = dot_product(delta_even(:,c), delta_odd(:,c))
        end do
        if( any(.not. ieee_is_finite(jscore)) )then
            deallocate(jscore)
            return
        endif
        selected = maxloc(jscore, dim=1)
        status = EXTREMA_CALIBRATION_PENDING
        if( validation_mode == EXTREMA_HEURISTIC .or. size(null_maxima) == 0 ) return
        status = EXTREMA_INVALID
        if( any(.not. ieee_is_finite(null_maxima)) ) return
        p_search = real(1 + count(null_maxima >= jscore(selected)),real64) / real(size(null_maxima)+1,real64)
        status = EXTREMA_OK
    end subroutine flex_extrema_cross_statistic

end module simple_flex_supported_extrema
