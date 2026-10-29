!@descr: closed-form and identity checks for same-factor ancestral initialization
module simple_flex_ancestral_init_tester
use, intrinsic :: iso_fortran_env, only: real64
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
use simple_flex_ancestral_init, only: flex_ancestral_mix
use simple_test_utils,          only: assert_true, assert_int
implicit none
private
public :: test_flex_ancestral_init

contains

    subroutine test_flex_ancestral_init()
        real(real64) :: leaf(6,2), ancestor(6,2), u(6), expected(6,2), exact_w(6,2), expected_t(2,2)
        real(real64), allocatable :: p(:,:), w(:,:), tr(:,:), prior(:)
        integer, allocatable :: labels(:)
        integer :: rows(6), parent_rows(6), eo(6), parent_eo(6), targets(2), parent_targets(2), status
        logical :: training(6), verified
        character(len=160) :: message
        rows = [11,17,23,29,31,47]
        parent_rows = rows
        eo = [0,1,0,1,0,1]
        parent_eo = eo
        targets = [1,2]
        parent_targets = [1,2]
        training = [.true.,.true.,.true.,.true.,.false.,.false.]
        verified = .true.
        leaf(:,1) = [1._real64,1._real64,0._real64,0._real64,1._real64,0._real64]
        leaf(:,2) = 1._real64 - leaf(:,1)
        ancestor(:,1) = [1._real64,0._real64,0._real64,0._real64,1._real64,0._real64]
        ancestor(:,2) = 1._real64 - ancestor(:,1)
        u = [.2_real64,.8_real64,.1_real64,.9_real64,.3_real64,.99_real64]
        expected(:,1) = [1._real64,26._real64/41,6._real64/41,6._real64/41,1._real64,6._real64/41]
        expected(:,2) = 1._real64 - expected(:,1)
        exact_w(:,1) = [19._real64/60,13._real64/60,1._real64/20,1._real64/20,19._real64/60,1._real64/20]
        exact_w(:,2) = [0._real64,1._real64/8,7._real64/24,7._real64/24,0._real64,7._real64/24]
        expected_t(1,:) = [1._real64,0._real64]
        expected_t(2,:) = [1._real64/3,2._real64/3]
        call mix(.5_real64)
        call assert_int(status, 0, 'ancestral blend accepts a verified partition pair')
        if( status /= 0 ) return
        call assert_true(maxval(abs(p - expected)) <= 2.e-14_real64, 'ancestral membership closed form')
        call assert_true(maxval(abs(w - exact_w)) <= 2.e-14_real64, 'ancestral column normalization')
        call assert_true(maxval(abs(tr - expected_t)) <= 2.e-14_real64, 'training-only target correspondence')
        call assert_true(all(labels == [1,2,1,2,1,2]), 'paired categorical labels use supplied draws')
        call mix(0._real64)
        call assert_int(status, 0, 'leaf-only control is valid')
        if( status == 0 ) call assert_true(maxval(abs(p - leaf)) <= 2.e-14_real64, 'leaf control is unchanged')
        leaf(5:6,:) = leaf(6:5:-1,:)
        call mix(.5_real64)
        call assert_int(status, 0, 'changed holdout rows are valid')
        if( status == 0 ) call assert_true(maxval(abs(tr - expected_t)) <= 2.e-14_real64, 'holdout cannot train mapping')
        leaf(5:6,:) = leaf(6:5:-1,:)
        parent_rows(2) = 18
        call mix(.5_real64)
        call assert_true(status /= 0 .and. .not. allocated(p), 'mismatched original rows return no initializer')
        parent_rows = rows
        parent_eo(1) = 1
        call mix(.5_real64)
        call assert_true(status /= 0, 'mismatched halves are rejected')
        parent_eo = eo
        verified = .false.
        call mix(.5_real64)
        call assert_true(status /= 0, 'unverified factor is rejected')
        verified = .true.
        targets = [1,1]
        call mix(.5_real64)
        call assert_true(status /= 0, 'duplicate target identities are rejected')
        targets = [1,2]
        leaf(1,:) = 0._real64
        call mix(.5_real64)
        call assert_true(status /= 0, 'zero-mass row is rejected without imputation')
        leaf(1,:) = [-1._real64,2._real64]
        call mix(.5_real64)
        call assert_true(status /= 0, 'signed weights are rejected')
        leaf(1,:) = [ieee_value(0._real64,ieee_quiet_nan),0._real64]
        call mix(.5_real64)
        call assert_true(status /= 0, 'non-finite weights are rejected')
        leaf = .5_real64
        ancestor = .5_real64
        u(1) = .5_real64
        u(2) = nearest(.5_real64, -1._real64)
        call mix(.5_real64)
        call assert_int(status, 0, 'exact categorical boundary fixture is valid')
        if( status == 0 )then
            call assert_int(labels(1), 2, 'equality at the CDF boundary selects the next target')
            call assert_int(labels(2), 1, 'draw below the CDF boundary selects the first target')
        endif

    contains

        subroutine mix( fraction )
            real(real64), intent(in) :: fraction
            call flex_ancestral_mix(leaf, ancestor, rows, parent_rows, eo, parent_eo, targets, parent_targets, &
                &training, u, fraction, 'fixture_site', 'fixture_site', verified, p, w, tr, prior, labels, status, message)
        end subroutine mix

    end subroutine test_flex_ancestral_init

end module simple_flex_ancestral_init_tester
