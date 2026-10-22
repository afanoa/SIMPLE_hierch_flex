!@descr: unit tests of the fixed-background focused PCG solve (reconstructor_pcg%set_fixed_background)
module simple_pcg_focus_tester
use simple_defs,              only: dp, PI
use simple_ori,               only: ori
use simple_oris,              only: oris
use simple_image,             only: image
use simple_reconstructor_pcg, only: reconstructor_pcg, PCG_OP_KERNEL, PCG_LAMBDA
use simple_test_utils,        only: assert_true
implicit none
private
public :: run_all_pcg_focus_tests

integer, parameter :: BOX          = 32
real,    parameter :: SMPD         = 1.0
real,    parameter :: SIGMA        = 1.5
real,    parameter :: FOCUS_RAD    = 5.0
real,    parameter :: FOCUS_EDGE   = 2.0
integer, parameter :: NVIEWS       = 120
integer, parameter :: MAXITS       = 25
real,    parameter :: CORE(3,3)    = reshape([-6., -2., 1., -3., 4., -2., -5., 1., 5.], [3,3])
real,    parameter :: REGION(3,2)  = reshape([6., 1., 0., 7., -2., 2.], [3,2])
real,    parameter :: FOCUS_CTR(3) = [6.5, -0.5, 1.0]
real(dp), parameter :: GAIN_MIN    = 2.0_dp

contains

    subroutine run_all_pcg_focus_tests()
        write(*,'(A)') '**** running all pcg_focus tests ****'
        call test_zero_background()
        call test_background_outside_mask()
        call test_region_recovery()
    end subroutine run_all_pcg_focus_tests

    subroutine test_zero_background()
        class(reconstructor_pcg), allocatable :: hard, focus
        real, allocatable :: x0(:,:,:), x1(:,:,:)
        write(*,'(A)') 'test_zero_background'
        allocate(hard, focus)
        allocate(x0(BOX,BOX,BOX), x1(BOX,BOX,BOX), source=0.)
        call new_fixture(hard)
        call new_fixture(focus, x1)
        call hard%solve_accum(x0, maxits=MAXITS, rtol=0.)
        call focus%solve_accum(x1, maxits=MAXITS, rtol=0.)
        call assert_true(all(x1 == x0), 'a zero background reproduces the hard support bit for bit')
        call hard%kill
        call focus%kill
    end subroutine test_zero_background

    subroutine test_background_outside_mask()
        class(reconstructor_pcg), allocatable :: focus
        real, allocatable :: xfix(:,:,:), x(:,:,:), m(:,:,:)
        write(*,'(A)') 'test_background_outside_mask'
        allocate(focus)
        xfix = blob_map(CORE)
        m    = focus_mask()
        call new_fixture(focus, xfix)
        x = xfix
        call focus%solve_accum(x, maxits=MAXITS, rtol=0.)
        call assert_true(all(pack(x, m == 0.) == pack(xfix, m == 0.)), 'the map is the background outside the focus mask')
        call assert_true(any(pack(x, m > 0.) /= pack(xfix, m > 0.)), 'the focus region is solved')
        call focus%kill
    end subroutine test_background_outside_mask

    subroutine test_region_recovery()
        class(reconstructor_pcg), allocatable :: hard, focus
        real, allocatable :: core_map(:,:,:), truth(:,:,:), xh(:,:,:), xf(:,:,:), m(:,:,:)
        real(dp) :: err_hard, err_focus
        write(*,'(A)') 'test_region_recovery'
        allocate(hard, focus)
        core_map = blob_map(CORE)
        truth    = core_map + blob_map(REGION)
        m        = focus_mask()
        call new_fixture(hard)
        call new_fixture(focus, core_map)
        allocate(xh(BOX,BOX,BOX), source=0.)
        xf = core_map
        call hard%solve_accum(xh, maxits=MAXITS, rtol=0.)
        call focus%solve_accum(xf, maxits=MAXITS, rtol=0.)
        err_hard  = masked_err(xh, truth, m)
        err_focus = masked_err(xf, truth, m)
        write(*,'(A,2F9.4)') '  relative error in the focus mask, hard support and fixed background: ', err_hard, err_focus
        call assert_true(err_focus * GAIN_MIN < err_hard, 'a known background removes the core leakage into the focus region')
        call hard%kill
        call focus%kill
    end subroutine test_region_recovery

    subroutine new_fixture( pcgop, xfix )
        class(reconstructor_pcg), intent(inout) :: pcgop
        real, optional,           intent(in)    :: xfix(:,:,:)
        type(oris)  :: os
        type(ori)   :: o
        type(image) :: mskvol
        complex, allocatable :: planes(:,:,:)
        integer :: i, lims2(2,2)
        call os%new(NVIEWS, .false.)
        call os%spiral()
        call pcgop%new(BOX, SMPD, PCG_LAMBDA)
        lims2 = pcgop%get_lims2()
        allocate(planes(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2),NVIEWS))
        call pcgop%set_volume((blob_map(CORE) + blob_map(REGION)) * pcgop%get_invenv())
        call o%new(.false.)
        do i = 1, NVIEWS
            call os%get_ori(i, o)
            call pcgop%forward_plane(o, planes(:,:,i))
        end do
        call o%kill()
        call mskvol%new([BOX,BOX,BOX], SMPD)
        call mskvol%set_rmat(focus_mask(), .false.)
        call pcgop%set_window_volume(mskvol)
        call mskvol%kill()
        if( present(xfix) ) call pcgop%set_fixed_background(xfix)
        call pcgop%prep_particles(os, use_ctf=.false.)
        call pcgop%begin_accum
        call pcgop%accumulate_batch(planes, NVIEWS, 1)
        call pcgop%end_accum(.true.)
        call pcgop%set_op_mode(PCG_OP_KERNEL)
        call os%kill()
    end subroutine new_fixture

    function blob_map( ctrs ) result( x )
        real, intent(in)  :: ctrs(:,:)
        real, allocatable :: x(:,:,:)
        real    :: c
        integer :: i, j, k, b
        c = real(BOX/2 + 1)
        allocate(x(BOX,BOX,BOX), source=0.)
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    do b = 1, size(ctrs,2)
                        x(i,j,k) = x(i,j,k) + exp(-((real(i) - c - ctrs(1,b))**2 + (real(j) - c - ctrs(2,b))**2 + &
                            &(real(k) - c - ctrs(3,b))**2) / (2. * SIGMA**2))
                    end do
                end do
            end do
        end do
    end function blob_map

    function focus_mask() result( m )
        real, allocatable :: m(:,:,:)
        real    :: c, r
        integer :: i, j, k
        c = real(BOX/2 + 1)
        allocate(m(BOX,BOX,BOX), source=0.)
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    r = sqrt((real(i) - c - FOCUS_CTR(1))**2 + (real(j) - c - FOCUS_CTR(2))**2 + &
                        &(real(k) - c - FOCUS_CTR(3))**2)
                    if( r <= FOCUS_RAD )then
                        m(i,j,k) = 1.
                    else if( r < FOCUS_RAD + FOCUS_EDGE )then
                        m(i,j,k) = 0.5 * (1. + cos(PI * (r - FOCUS_RAD) / FOCUS_EDGE))
                    endif
                end do
            end do
        end do
    end function focus_mask

    real(dp) function masked_err( x, truth, m )
        real, intent(in) :: x(:,:,:), truth(:,:,:), m(:,:,:)
        masked_err = sqrt(sum(real((m * (x - truth))**2, dp)) / sum(real((m * truth)**2, dp)))
    end function masked_err

end module simple_pcg_focus_tester
