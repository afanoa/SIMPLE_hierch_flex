!@descr: unit tests of the per-site split's geometry and call (simple_site_split): the convex hull of a voxel set, and GAPw on synthetic two-lobe versus continuous-rod maps
! Noise-free identical halves make the first factor signal-dominated, so the call rests on the lobes and the density between them.
module simple_site_split_tester
use simple_test_utils
use simple_site_split, only: gap_call, gap_stats, hull_mask
implicit none
private
public :: run_all_site_split_tests

integer, parameter :: BOX = 40, XL = 14, XR = 26, YC = 20, ZC = 20
real,    parameter :: SMPD = 5., SIG = 2., RSITE = 12.

contains

    subroutine run_all_site_split_tests()
        write(*,'(A)') '**** running all site split tests ****'
        call test_hull()
        call test_gap_discrete_vs_continuum()
    end subroutine run_all_site_split_tests

    subroutine test_hull()
        logical, allocatable :: mk(:,:,:), inside(:,:,:)
        allocate(mk(BOX,BOX,BOX), source=.false.)
        mk(10:13,18:21,18:21) = .true.
        mk(27:30,18:21,18:21) = .true.
        call hull_mask(mk, inside)
        call assert_true(count(mk .and. .not. inside) == 0, 'the lobes lie inside their convex hull')
        call assert_true(inside(20,19,19), 'the segment between the lobes lies inside the hull')
        call assert_true(.not. inside(20,30,19), 'an off-axis voxel lies outside the hull')
        call assert_true(.not. inside(5,19,19), 'a voxel beyond the lobes lies outside the hull')
        call assert_int(21 * 4 * 4, count(inside), 'the hull of two aligned cubes is their bounding block')
        deallocate(mk)
        allocate(mk(BOX,BOX,BOX), source=.false.)
        mk(10:30,18:21,20) = .true.
        call hull_mask(mk, inside)
        call assert_true(.not. any(inside), 'a planar voxel set has no hull')
    end subroutine test_hull

    subroutine test_gap_discrete_vs_continuum()
        real,    allocatable :: bl(:,:,:), br(:,:,:), rod(:,:,:), e1(:,:,:), e2(:,:,:)
        logical, allocatable :: m(:,:,:)
        type(gap_stats) :: st
        integer :: i, j, k
        allocate(bl(BOX,BOX,BOX), br(BOX,BOX,BOX), rod(BOX,BOX,BOX), m(BOX,BOX,BOX))
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    bl(i,j,k)  = blob(i, j, k, XL)
                    br(i,j,k)  = blob(i, j, k, XR)
                    rod(i,j,k) = 0.
                    if( i >= XL .and. i <= XR ) rod(i,j,k) = exp(-real((j-YC)**2 + (k-ZC)**2) / (2. * SIG**2))
                    m(i,j,k) = real((i-20)**2 + (j-YC)**2 + (k-ZC)**2) <= RSITE**2
                end do
            end do
        end do
        ! two discrete states: the blob sits left or right, nothing between
        e1 = bl
        e2 = br
        call gap_call(e1, e2, e1, e2, 1000., 1000., m, SMPD, st)
        call assert_true(st%z1 > 3. .and. st%r1 > 0.5, 'identical halves: the first factor is signal-dominated')
        call assert_true(st%lsep > 10., 'the two blobs give lobes more than 10 A apart')
        call assert_true(st%ang < 1., 'identical halves: the lobe axes agree')
        call assert_true(st%nb >= 10, 'the hull between the blobs holds voxels beyond the dilated lobes')
        call assert_true(st%gapw < 0.28, 'two separated blobs: the density between the lobes is empty')
        call assert_char('DISCRETE', trim(st%call), 'two separated blobs are called discrete')
        ! a continuum: the same blobs joined by a rod of density
        e1 = bl + rod
        e2 = br + rod
        call gap_call(e1, e2, e1, e2, 1000., 1000., m, SMPD, st)
        call assert_true(st%nb >= 10, 'the rod case keeps a region between the lobes')
        call assert_true(st%gapw >= 0.28, 'a rod between the blobs fills the gap')
        call assert_char('not_discrete', trim(st%call), 'joined blobs are not called discrete')

        contains

            real function blob( i, j, k, x0 )
                integer, intent(in) :: i, j, k, x0
                blob = exp(-real((i-x0)**2 + (j-YC)**2 + (k-ZC)**2) / (2. * SIG**2))
            end function blob

    end subroutine test_gap_discrete_vs_continuum

end module simple_site_split_tester
