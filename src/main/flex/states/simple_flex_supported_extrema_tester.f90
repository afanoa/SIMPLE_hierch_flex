!@descr: support, identity, fission arithmetic and conditional endpoint statistic tests
module simple_flex_supported_extrema_tester
use, intrinsic :: iso_fortran_env, only: real64
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
use simple_flex_supported_extrema, only: flex_supported_tail_labels, flex_fission_channels, &
    &flex_extrema_cross_statistic, EXTREMA_OK, EXTREMA_INVALID, EXTREMA_POPULATION_LIMITED, &
    &EXTREMA_HALF_LIMITED, EXTREMA_CALIBRATION_PENDING, EXTREMA_HEURISTIC, &
    &EXTREMA_INDEPENDENT_NOISE, EXTREMA_REFITTED_NULL
use simple_test_utils, only: assert_true, assert_int
implicit none
private
public :: test_flex_supported_extrema

contains

    subroutine test_flex_supported_extrema()
        real(real64) :: score(18), reversed(18), y(2,3), eta(2,3), de(2,3), do_(2,3), nulls(4), p
        real(real64), allocatable :: sel(:,:), meas(:,:), js(:)
        integer :: ids(18), eo(18), order(18), counts(3,2), status, i, selected
        integer, allocatable :: labels(:), permuted(:)
        character(len=160) :: message
        ids = [(10*i,i=1,18)]
        eo = [(mod(i,2),i=1,18)]
        score = [(real(i,real64),i=1,18)]
        call flex_supported_tail_labels(score, ids, eo, 6, 6, 3, labels, counts, status, message)
        call assert_int(status, EXTREMA_OK, 'supported extrema accepts exact boundary populations')
        if( status /= EXTREMA_OK ) return
        call assert_true(all(labels(1:6) == 1), 'lower supported tail')
        call assert_true(all(labels(7:12) == 2), 'all middle rows retained')
        call assert_true(all(labels(13:18) == 3), 'upper supported tail')
        call assert_true(all(counts == 3), 'original half populations are retained')
        call assert_int(sum(counts), 18, 'no particle deletion or duplicated support')
        call flex_supported_tail_labels(score, ids, 1-eo, 6, 6, 3, permuted, counts, status, message)
        call assert_int(status, EXTREMA_OK, 'EO code swap remains valid')
        if( status == EXTREMA_OK ) call assert_true(all(permuted == labels), 'EO code swap does not change labels')
        score = [(real((i-1)/2,real64),i=1,18)]
        call flex_supported_tail_labels(score, ids, eo, 5, 5, 2, labels, counts, status, message)
        call assert_int(status, EXTREMA_OK, 'tied boundary fixture remains supported')
        if( status /= EXTREMA_OK ) return
        order = [(19-i,i=1,18)]
        reversed = score(order)
        call flex_supported_tail_labels(reversed, ids(order), eo(order), 5, 5, 2, permuted, counts, status, message)
        call assert_int(status, EXTREMA_OK, 'permuted input remains supported')
        if( status == EXTREMA_OK ) call assert_true(all(permuted == labels(order)), 'ties use original IDs, not row order')
        call flex_supported_tail_labels(score, ids, eo, 5, 6, 2, labels, counts, status, message)
        call assert_int(status, EXTREMA_POPULATION_LIMITED, 'undersized endpoint is rejected')
        call assert_true(.not. allocated(labels), 'failed population gate publishes no labels')
        call flex_supported_tail_labels(score, ids, eo, 7, 5, 2, labels, counts, status, message)
        call assert_int(status, EXTREMA_POPULATION_LIMITED, 'undersized middle cannot be published as a state')
        call flex_supported_tail_labels(score, ids, 0*eo, 6, 6, 3, labels, counts, status, message)
        call assert_int(status, EXTREMA_HALF_LIMITED, 'half loss is rejected without rebalancing')
        call assert_true(.not. allocated(labels), 'failed half gate publishes no labels')
        ids(2) = ids(1)
        call flex_supported_tail_labels(score, ids, eo, 6, 6, 3, labels, counts, status, message)
        call assert_int(status, EXTREMA_INVALID, 'duplicate physical row IDs are rejected')
        ids(2) = 20
        score = 0._real64
        call flex_supported_tail_labels(score, ids, eo, 6, 6, 3, labels, counts, status, message)
        call assert_int(status, EXTREMA_INVALID, 'constant score cannot manufacture extrema by ID')
        score(1) = ieee_value(0._real64, ieee_quiet_nan)
        call flex_supported_tail_labels(score, ids, eo, 6, 6, 3, labels, counts, status, message)
        call assert_int(status, EXTREMA_INVALID, 'nonfinite physical score is rejected')

        y = reshape([1._real64,2._real64,3._real64,4._real64,5._real64,6._real64],[2,3])
        eta = reshape([.5_real64,-.5_real64,1._real64,-1._real64,2._real64,-2._real64],[2,3])
        call flex_fission_channels(y, eta, sel, meas, status)
        call assert_int(status, EXTREMA_OK, 'fission arithmetic with supplied independent draw')
        if( status == EXTREMA_OK )then
            call assert_true(all(sel == y+eta), 'selection channel uses the positive draw')
            call assert_true(all(meas == y-eta), 'measurement channel uses the same negative draw')
            call assert_true(all(.5_real64*(sel+meas) == y), 'fission preserves original observations')
        endif
        call flex_fission_channels(y, eta(:,1:2), sel, meas, status)
        call assert_int(status, EXTREMA_INVALID, 'fission rejects mismatched arrays')
        call assert_true(.not. allocated(sel) .and. .not. allocated(meas), 'invalid fission leaves no channels')

        de(:,1) = [2._real64,0._real64]
        de(:,2) = [0._real64,1._real64]
        de(:,3) = [1._real64,1._real64]
        do_(:,1) = [-1._real64,0._real64]
        do_(:,2) = [0._real64,2._real64]
        do_(:,3) = [1._real64,1._real64]
        nulls = [-1._real64,1._real64,2._real64,3._real64]
        call flex_extrema_cross_statistic(de, do_, nulls, EXTREMA_HEURISTIC, js, selected, p, status)
        call assert_int(status, EXTREMA_CALIBRATION_PENDING, 'same-feature/unknown provenance never returns a p value')
        call assert_true(p < 0._real64, 'unvalidated contrast is not an acceptance probability')
        if( allocated(js) ) call assert_true(all(js == [-2._real64,2._real64,2._real64]), 'negative J is not clipped')
        call assert_int(selected, 2, 'candidate ties retain frozen input order')
        call flex_extrema_cross_statistic(de, do_, nulls, EXTREMA_INDEPENDENT_NOISE, js, selected, p, status)
        call assert_int(status, EXTREMA_OK, 'declared independent-noise null statistic')
        call assert_true(abs(p-.6_real64) <= 1.e-15_real64, 'null equality included and plus-one correction retained')
        call flex_extrema_cross_statistic(de, do_, nulls(:0), EXTREMA_REFITTED_NULL, js, selected, p, status)
        call assert_int(status, EXTREMA_CALIBRATION_PENDING, 'empty refitted null remains uncalibrated')
        nulls(2) = ieee_value(0._real64, ieee_quiet_nan)
        call flex_extrema_cross_statistic(de, do_, nulls, EXTREMA_REFITTED_NULL, js, selected, p, status)
        call assert_int(status, EXTREMA_INVALID, 'nonfinite null values are rejected')
        call assert_true(p < 0._real64, 'invalid null cannot publish a p value')
        call test_state_set_consumer
    end subroutine test_flex_supported_extrema

    subroutine test_state_set_consumer()
        use simple_flex_pca_records, only: flex_latent, flex_state_set
        use simple_flex_pca_weights, only: build_supported_extrema_state_weights
        class(flex_latent), allocatable :: latent
        class(flex_state_set), allocatable :: states
        integer :: ids(18), eo(18), counts(3,2), status, calibration, i
        real(real64) :: score(18)
        character(len=160) :: message
        allocate(latent, states)
        allocate(latent%z(18,2))
        ids = [(10*i,i=1,18)]
        eo = [(mod(i,2),i=1,18)]
        score = [(real(i,real64),i=1,18)]
        latent%z(:,1) = score
        latent%z(:,2) = 2._real64 * score
        call build_supported_extrema_state_weights(latent, score, ids, eo, 6, 6, 3, .true., &
            &states, counts, status, calibration, message)
        call assert_int(status, EXTREMA_OK, 'canonical state-set consumer accepts supported tails')
        call assert_int(calibration, EXTREMA_CALIBRATION_PENDING, 'candidate support never certifies physical calibration')
        if( status == EXTREMA_OK )then
            call assert_int(states%nstates, 3, 'canonical state set retains A, middle and B')
            call assert_true(all(sum(states%weights,dim=2) == 1.), 'canonical unit weights preserve every row')
            call assert_true(all(states%neff == 6.), 'unit-weight ESS equals distinct supported counts')
            call assert_true(all(states%targets(1,:) == [3.5,9.5,15.5]), 'targets describe delivered group means')
            call assert_true(all(states%bandwidths == 0.), 'hard partition has no smoothing bandwidth')
            call assert_true(.not. allocated(states%kdist), 'kernel CV distances are not fabricated')
        endif
        call build_supported_extrema_state_weights(latent, score, ids, eo, 6, 6, 3, .false., &
            &states, counts, status, calibration, message)
        call assert_int(status, EXTREMA_INVALID, 'unverified physical score frame blocks consumer')
        call assert_int(states%nstates, 3, 'failed proposal leaves preceding state set intact')
        call build_supported_extrema_state_weights(latent, score, ids, eo, 7, 6, 3, .true., &
            &states, counts, status, calibration, message)
        call assert_int(status, EXTREMA_POPULATION_LIMITED, 'consumer refuses unsupported middle')
        call states%kill
        call latent%kill
    end subroutine test_state_set_consumer

end module simple_flex_supported_extrema_tester
