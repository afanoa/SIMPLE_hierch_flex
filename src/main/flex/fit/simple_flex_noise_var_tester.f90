!@descr: flex noise-variance kernels under a forced team against one thread: the residual moments and the band whitening
module simple_flex_noise_var_tester
use simple_core_module_api, only: dp, fplane_type, logfhandle
use simple_flex_noise_var,  only: accumulate_moments
use simple_flex_pca_planes, only: flex_plane_store
use simple_test_utils,      only: assert_true
!$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
implicit none
private
public :: run_all_flex_noise_var_tests

integer, parameter :: BOX = 16, NH = BOX / 2 + 1, KB = 6, NP = 37, NTEAM = 4

contains

    subroutine run_all_flex_noise_var_tests()
        write(logfhandle,'(A)') '-- test_noise_var_moments_team'
        call test_noise_var_moments_team()
        write(logfhandle,'(A)') '-- test_noise_var_whiten_team'
        call test_noise_var_whiten_team()
    end subroutine run_all_flex_noise_var_tests

    subroutine test_noise_var_moments_team()
        type(fplane_type), allocatable :: fpls(:)
        complex(dp),       allocatable :: mu(:)
        real(dp),          allocatable :: f2(:,:,:), aa(:,:,:), bb(:,:,:)
        real    :: rmats(3,3,NP)
        integer :: hs((2*KB+1)**2), ks((2*KB+1)**2), rs((2*KB+1)**2), ns, h, k, i, run, nthr_prev
        call make_planes(fpls)
        allocate(mu(NH*BOX*BOX))
        do i = 1, size(mu)
            mu(i) = cmplx(sin(0.37d0 * i), cos(0.91d0 * i), kind=dp)
        end do
        do i = 1, NP
            rmats(:,:,i) = rotation(0.7 * i, 1.3 * i, 2.1 * i)
        end do
        ns = 0
        do k = -KB, 0
            do h = -KB, KB
                if( k == 0 .and. h <= 0 ) cycle
                if( h * h + k * k < 4 .or. h * h + k * k > KB * (KB + 1) ) cycle
                ns = ns + 1
                hs(ns) = h; ks(ns) = k; rs(ns) = nint(sqrt(real(h * h + k * k)))
            end do
        end do
        allocate(f2(0:KB,NP,2), aa(0:KB,NP,2), bb(0:KB,NP,2), source=0.d0)
        nthr_prev = 1
        !$ nthr_prev = omp_get_max_threads()
        do run = 1, 2
            !$ call omp_set_num_threads(merge(1, NTEAM, run == 1))
            call accumulate_moments(fpls, rmats, hs(:ns), ks(:ns), rs(:ns), mu, NH, BOX, f2(:,:,run), aa(:,:,run), bb(:,:,run))
        end do
        !$ call omp_set_num_threads(nthr_prev)
        call assert_true(all(f2(:,:,1) == f2(:,:,2)) .and. all(aa(:,:,1) == aa(:,:,2)) .and. all(bb(:,:,1) == bb(:,:,2)), &
            &'residual moments: a team of 4 gives the one-thread sums bit for bit')
        call assert_true(all(f2(2:,:,1) > 0.d0) .and. all(bb(2:,:,1) > 0.d0), 'residual moments: every ring holds data and model')
    end subroutine test_noise_var_moments_team

    subroutine test_noise_var_whiten_team()
        class(flex_plane_store), allocatable :: store
        type(fplane_type),       allocatable :: raw(:), one(:), team(:)
        real    :: band_var(6,NP+3), f
        integer :: band_of(0:2*BOX), rows(NP), i, b, h, nthr_prev
        logical :: l_same
        allocate(store)
        do i = 1, NP + 3
            do b = 1, 6
                band_var(b,i) = 0.5 + abs(sin(0.61 * real(b + 6 * i)))
            end do
        end do
        do h = 0, 2 * BOX
            band_of(h) = max(1, min(6, h / 2))
        end do
        rows = [(NP + 4 - i, i = 1, NP)]
        call store%set_noise_bands(band_var, band_of)
        call make_planes(raw)
        call make_planes(one)
        call make_planes(team)
        nthr_prev = 1
        !$ nthr_prev = omp_get_max_threads()
        !$ call omp_set_num_threads(1)
        call store%whiten_bands(rows, one)
        !$ call omp_set_num_threads(NTEAM)
        call store%whiten_bands(rows, team)
        !$ call omp_set_num_threads(nthr_prev)
        l_same = .true.
        do i = 1, NP
            l_same = l_same .and. all(team(i)%cmplx_plane == one(i)%cmplx_plane) .and. &
                &all(team(i)%transfer_plane == one(i)%transfer_plane) .and. all(team(i)%ctfsq_plane == one(i)%ctfsq_plane)
        end do
        call assert_true(l_same, 'band whitening: a team of 4 gives the one-thread planes bit for bit')
        f = 1. / sqrt(band_var(band_of(4),rows(1)))
        call assert_true(abs(one(1)%cmplx_plane(3,-2) - f * raw(1)%cmplx_plane(3,-2)) < 1.e-6 .and. &
            &abs(one(1)%ctfsq_plane(3,-2) - f * f * raw(1)%ctfsq_plane(3,-2)) < 1.e-6, &
            &'band whitening: a sample is scaled by 1 / sqrt(v) of its row and band')
        call store%kill
    end subroutine test_noise_var_whiten_team

    subroutine make_planes( fpls )
        type(fplane_type), allocatable, intent(out) :: fpls(:)
        integer :: i, h, k, n
        allocate(fpls(NP))
        do i = 1, NP
            allocate(fpls(i)%cmplx_plane(-KB:KB,-KB:0), fpls(i)%transfer_plane(-KB:KB,-KB:0), fpls(i)%ctfsq_plane(-KB:KB,-KB:0))
            do k = -KB, 0
                do h = -KB, KB
                    n = i * 1000 + (k + KB) * 100 + h + KB
                    fpls(i)%cmplx_plane(h,k)    = cmplx(sin(0.13 * n), cos(0.29 * n))
                    fpls(i)%transfer_plane(h,k) = cmplx(0.2 + abs(sin(0.07 * n)), 0.)
                    fpls(i)%ctfsq_plane(h,k)    = real(fpls(i)%transfer_plane(h,k))**2
                end do
            end do
        end do
    end subroutine make_planes

    pure function rotation( a, b, c ) result( r )
        real, intent(in) :: a, b, c
        real :: r(3,3), za(3,3), yb(3,3), zc(3,3)
        za = reshape([cos(a), sin(a), 0., -sin(a), cos(a), 0., 0., 0., 1.], [3,3])
        yb = reshape([cos(b), 0., -sin(b), 0., 1., 0., sin(b), 0., cos(b)], [3,3])
        zc = reshape([cos(c), sin(c), 0., -sin(c), cos(c), 0., 0., 0., 1.], [3,3])
        r  = matmul(za, matmul(yb, zc))
    end function rotation

end module simple_flex_noise_var_tester
