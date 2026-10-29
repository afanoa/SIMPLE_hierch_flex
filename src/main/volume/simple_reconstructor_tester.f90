!@descr: unit tests of weighted gridding and shared projection/insertion geometry
! Linearity: a plane inserted with weight w equals the plane inserted with weight one and the
! accumulators scaled by w, to round-off; weight one is the unweighted insertion bit for bit.
module simple_reconstructor_tester
!$ use omp_lib,             only: omp_get_thread_num, omp_get_dynamic, omp_set_dynamic
!$ use omp_lib,             only: omp_get_max_threads, omp_set_num_threads
use iso_fortran_env,        only: int32
use simple_defs,            only: OSMPL_PAD_FAC
use simple_image,           only: image
use simple_fplane_grid,     only: fplane_grid
use simple_ori,             only: ori
use simple_oris,            only: oris
use simple_parameters,      only: parameters
use simple_sp_project,      only: sp_project
use simple_sym,             only: sym
use simple_memoize_ft_maps, only: forget_ft_maps, memoize_ft_maps
use simple_type_defs,       only: CTFFLAG_NO, ctfparams, fplane_type
use simple_reconstructor,   only: reconstructor, exp_samples, insert_planes_targets, MAX_INSERT_TARGETS, MAX_INSERT_PLANES
use simple_test_utils
implicit none
private
public :: run_all_reconstructor_tests

integer, parameter :: BOX     = 16
integer, parameter :: NPLANES = 12
integer, parameter :: SEED    = 20261007
real,    parameter :: SMPD    = 2.0
real,    parameter :: RELTOL  = 1.e-6

contains

    subroutine run_all_reconstructor_tests()
        write(*,'(A)') '**** running all reconstructor tests ****'
        call test_weighted_insertion_linearity()
        call test_projection_geometry_reuse()
        call test_insertion_geometry_reuse()
    end subroutine run_all_reconstructor_tests

    subroutine test_insertion_geometry_reuse()
        integer, parameter :: NSLOTS = 6, MAXTEAM = 4, SLOTS(4) = [5,2,6,1], NINSERT = MAX_INSERT_PLANES+3
        real, parameter :: WW(6) = [0., 0.05, 0.37, 1., 1.25, -0.1]
        character(len=3), parameter :: SYMMETRIES(2) = ['c1 ', 'c11']
        class(parameters), target, allocatable :: params
        class(sp_project),         allocatable :: project
        class(reconstructor),      allocatable :: direct(:), shared(:)
        class(ori),                allocatable :: poses(:)
        class(sym),                allocatable :: symmetry
        class(image),              allocatable :: obs, obs_pad
        class(fplane_grid),        allocatable :: grid
        type(fplane_type) :: template, planes(MAX_INSERT_PLANES)
        type(ctfparams) :: ctfparms
        real :: img(BOX,BOX,1), weights(MAX_INSERT_TARGETS,MAX_INSERT_PLANES)
        complex, allocatable :: serial_cmat(:,:,:,:)
        real, allocatable :: serial_rho(:,:,:,:)
        integer :: lo(3), hi(3), seen(MAXTEAM), saved_team, wanted, nt, team, isym, pass, i, q, h, k, it, begin, ib, nr
        !$ logical :: saved_dynamic
        write(*,'(A)') 'test_insertion_geometry_reuse'
        allocate(params, project, direct(NSLOTS), shared(NSLOTS), poses(MAX_INSERT_PLANES), symmetry, obs, obs_pad, grid)
        params%box        = BOX
        params%box_crop   = BOX
        params%box_croppd = OSMPL_PAD_FAC * BOX
        params%smpd_crop  = SMPD
        params%nstates    = 1
        params%numlen     = 1
        params%oritype    = 'cls3D'
        do q = 1, NSLOTS
            call direct(q)%new_accumulator(params, project, expand=.true., wthreads=.false.)
            call shared(q)%new_accumulator(params, project, expand=.true., wthreads=.false.)
        end do
        lo = lbound(direct(1)%cmat_exp)
        hi = ubound(direct(1)%cmat_exp)
        allocate(serial_cmat(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),MAX_INSERT_TARGETS))
        allocate(serial_rho(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),MAX_INSERT_TARGETS))
        do ib = 1, MAX_INSERT_PLANES
            call poses(ib)%new(.true.)
        end do
        call obs%new([BOX,BOX,1], SMPD, wthreads=.false.)
        call obs_pad%new([OSMPL_PAD_FAC*BOX,OSMPL_PAD_FAC*BOX,1], SMPD, wthreads=.false.)
        call grid%new(OSMPL_PAD_FAC*BOX)
        call memoize_ft_maps([OSMPL_PAD_FAC*BOX,OSMPL_PAD_FAC*BOX,1], SMPD)
        do k = 1, BOX
            do h = 1, BOX
                img(h,k,1) = sin(real(3*h+7*k))
            end do
        end do
        call obs%set_rmat(img, .false.)
        call obs%pad(obs_pad, backgr=0., antialiasing=.false.)
        call obs_pad%fft()
        ctfparms%smpd    = SMPD
        ctfparms%ctfflag = CTFFLAG_NO
        call obs_pad%gen_fplane4rec(SMPD, ctfparms, [0.,0.], template, grid)
        saved_team = 1
        !$ saved_team = omp_get_max_threads()
        !$ saved_dynamic = omp_get_dynamic()
        !$ call omp_set_dynamic(.false.)
        do isym = 1, size(SYMMETRIES)
            call symmetry%new(trim(SYMMETRIES(isym)))
            do nt = 2, MAX_INSERT_TARGETS, 2
                do team = 1, MAXTEAM, MAXTEAM-1
                    !$ call omp_set_num_threads(team)
                    seen = 0
                    !$omp parallel default(shared) private(it)
                    it = 1
                    !$ it = omp_get_thread_num()+1
                    seen(it) = 1
                    !$omp end parallel
                    wanted = 1
                    !$ wanted = team
                    call assert_true(sum(seen) == wanted, 'insertion reuse exercises the requested team')
                    do q = 1, NSLOTS
                        call direct(q)%reset_exp()
                        call shared(q)%reset_exp()
                    end do
                    do pass = 1, 2
                        do begin = 1, NINSERT, MAX_INSERT_PLANES
                            nr = min(MAX_INSERT_PLANES,NINSERT-begin+1)
                            do ib = 1, nr
                                i = begin+ib-1
                                call poses(ib)%set_euler([real(19*i), real(11*i+pass), real(7*i)])
                                planes(ib) = template
                                planes(ib)%cmplx_plane = template%cmplx_plane * real(i+pass)/16.
                                do k = lbound(planes(ib)%ctfsq_plane,2), ubound(planes(ib)%ctfsq_plane,2)
                                    do h = lbound(planes(ib)%ctfsq_plane,1), ubound(planes(ib)%ctfsq_plane,1)
                                        planes(ib)%ctfsq_plane(h,k) = real(1+abs(h+k))/16.
                                        if( pass == 2 .and. h*h+k*k > 9 )then
                                            planes(ib)%cmplx_plane(h,k) = cmplx(0.,0.)
                                            planes(ib)%ctfsq_plane(h,k) = 0.
                                        endif
                                    end do
                                end do
                                do q = 1, nt
                                    weights(q,ib) = WW(mod(i+q+pass-3,size(WW))+1)
                                    if( weights(q,ib) <= 0. ) cycle
                                    call direct(SLOTS(q))%insert_plane_oversamp(symmetry, poses(ib), planes(ib), weights(q,ib))
                                end do
                            end do
                            call insert_planes_targets(shared, SLOTS(:nt), symmetry, poses, planes, weights(:nt,:nr), nr)
                        end do
                        do q = 1, NSLOTS
                            call assert_true(all(transfer(direct(q)%cmat_exp,[0_int32]) == &
                                &transfer(shared(q)%cmat_exp,[0_int32])), 'shared insertion data match bit for bit')
                            call assert_true(all(transfer(direct(q)%rho_exp,[0_int32]) == &
                                &transfer(shared(q)%rho_exp,[0_int32])), 'shared insertion density matches bit for bit')
                        end do
                    end do
                    do q = 1, nt
                        if( team == 1 )then
                            serial_cmat(:,:,:,q) = shared(SLOTS(q))%cmat_exp
                            serial_rho(:,:,:,q)  = shared(SLOTS(q))%rho_exp
                        else
                            call assert_true(all(transfer(serial_cmat(:,:,:,q),[0_int32]) == &
                                &transfer(shared(SLOTS(q))%cmat_exp,[0_int32])), 'insertion data match the serial team')
                            call assert_true(all(transfer(serial_rho(:,:,:,q),[0_int32]) == &
                                &transfer(shared(SLOTS(q))%rho_exp,[0_int32])), 'insertion density matches the serial team')
                        endif
                    end do
                end do
            end do
            call symmetry%kill()
        end do
        !$ call omp_set_num_threads(saved_team)
        !$ call omp_set_dynamic(saved_dynamic)
        do q = 1, NSLOTS
            call direct(q)%kill()
            call shared(q)%kill()
        end do
        do ib = 1, MAX_INSERT_PLANES
            call poses(ib)%kill()
        end do
        call obs%kill()
        call obs_pad%kill()
        call grid%kill()
        call forget_ft_maps()
    end subroutine test_insertion_geometry_reuse

    subroutine test_projection_geometry_reuse()
        integer, parameter :: NGROUPS = 3, MAXTEAM = 4, BANDS(3) = [4,6,3], MAXBAND = maxval(BANDS)
        class(parameters), target, allocatable :: params
        class(sp_project),         allocatable :: project
        class(reconstructor),      allocatable :: refs(:)
        class(ori),                allocatable :: poses(:)
        class(image),              allocatable :: obs, obs_pad
        class(fplane_grid),        allocatable :: grid
        type(exp_samples), allocatable :: samples(:)
        type(fplane_type), allocatable :: planes(:), output(:)
        type(fplane_type) :: template, direct
        type(ctfparams)   :: ctfparms
        complex, allocatable :: expected(:,:,:,:), values(:,:)
        logical :: failed(NPLANES,NGROUPS)
        !$ logical :: saved_dynamic
        integer :: seen(MAXTEAM), lo(3), hi(3), plo(2), phi(2), team, wanted, it, i, g, h, k, l, pass
        write(*,'(A)') 'test_projection_geometry_reuse'
        allocate(params, project, refs(NGROUPS), poses(NPLANES), obs, obs_pad, grid)
        allocate(samples(MAXTEAM), planes(NPLANES), output(MAXTEAM))
        params%box        = BOX
        params%box_crop   = BOX
        params%box_croppd = OSMPL_PAD_FAC * BOX
        params%smpd_crop  = SMPD
        params%nstates    = 1
        params%numlen     = 1
        params%oritype    = 'cls3D'
        do g = 1, NGROUPS
            call refs(g)%new_accumulator(params, project, expand=.true., wthreads=.false.)
        end do
        call obs%new([BOX,BOX,1], SMPD, wthreads=.false.)
        call obs_pad%new([OSMPL_PAD_FAC*BOX,OSMPL_PAD_FAC*BOX,1], SMPD, wthreads=.false.)
        call obs%zero_and_unflag_ft()
        call obs%pad(obs_pad, backgr=0., antialiasing=.false.)
        call obs_pad%fft()
        call grid%new(OSMPL_PAD_FAC*BOX)
        call memoize_ft_maps([OSMPL_PAD_FAC*BOX,OSMPL_PAD_FAC*BOX,1], SMPD)
        ctfparms%smpd    = SMPD
        ctfparms%ctfflag = CTFFLAG_NO
        call obs_pad%gen_fplane4rec(SMPD, ctfparms, [0.,0.], template, grid)
        do i = 1, NPLANES
            call poses(i)%new(.true.)
            call poses(i)%set_euler([real(19*i), real(11*i), real(7*i)])
            planes(i) = template
            planes(i)%nyq = OSMPL_PAD_FAC * BANDS(mod(i-1,size(BANDS))+1)
        end do
        plo = lbound(template%cmplx_plane)
        phi = ubound(template%cmplx_plane)
        allocate(expected(plo(1):phi(1),plo(2):phi(2),NPLANES,NGROUPS))
        allocate(values((2*MAXBAND+1)*(MAXBAND+1),MAXTEAM))
        lo = lbound(refs(1)%cmat_exp)
        hi = ubound(refs(1)%cmat_exp)
        !$ saved_dynamic = omp_get_dynamic()
        !$ call omp_set_dynamic(.false.)
        do pass = 1, 2
            do g = 1, NGROUPS
                do l = lo(3), hi(3)
                    do k = lo(2), hi(2)
                        do h = lo(1), hi(1)
                            refs(g)%cmat_exp(h,k,l) = cmplx(real(pass*g*h-k)/16., real(l+2*g)/32.)
                        end do
                    end do
                end do
            end do
            ! Fresh map values are gathered even when the same geometry is revisited.
            do i = 1, NPLANES
                do g = 1, NGROUPS
                    call refs(g)%project_fplane(poses(i), planes(i), direct, apply_ctf_amp=.false.)
                    expected(:,:,i,g) = direct%cmplx_plane
                end do
            end do
            do team = 1, MAXTEAM, MAXTEAM-1
                failed = .false.
                seen   = 0
                !$omp parallel do default(shared) private(i,it,g) num_threads(team) schedule(static)
                do i = 1, NPLANES
                    it = 1
                    !$ it = omp_get_thread_num() + 1
                    seen(it) = 1
                    call samples(it)%new_plane(refs(1), poses(i), planes(i))
                    do g = 1, NGROUPS
                        if( g == 1 .and. mod(i,3) == 0 ) cycle
                        call samples(it)%gather(refs(g), values(:,it))
                        call samples(it)%put_plane(values(:,it), planes(i), output(it), .false.)
                        failed(i,g) = any(transfer(output(it)%cmplx_plane, [0_int32]) /= &
                            &transfer(expected(:,:,i,g), [0_int32]))
                    end do
                end do
                !$omp end parallel do
                wanted = 1
                !$ wanted = team
                call assert_true(sum(seen) == wanted, 'projection reuse exercises the requested OpenMP team')
                call assert_true(.not. any(failed), 'reused projection geometry matches the original path bit for bit')
            end do
        end do
        !$ call omp_set_dynamic(saved_dynamic)
        do it = 1, MAXTEAM
            call samples(it)%kill()
        end do
        do g = 1, NGROUPS
            call refs(g)%kill()
        end do
        do i = 1, NPLANES
            call poses(i)%kill()
        end do
        call obs%kill()
        call obs_pad%kill()
        call grid%kill()
        call forget_ft_maps()
    end subroutine test_projection_geometry_reuse

    subroutine test_weighted_insertion_linearity()
        real, parameter :: WEIGHTS(3) = [0.37, 0.05, 0.81]
        type(parameters), target :: params
        type(sp_project)         :: project
        type(reconstructor)      :: rw, rs, rone, rplain
        type(oris)               :: poses
        type(ori)                :: o
        type(sym)                :: c1sym
        type(ctfparams)          :: ctfparms
        type(fplane_type)        :: fplane
        type(fplane_grid)        :: grid
        type(image)              :: obs, obs_pad
        complex, allocatable     :: cw(:,:,:), cs(:,:,:), cone(:,:,:), cplain(:,:,:)
        real,    allocatable     :: img(:,:,:), dw(:,:,:), ds(:,:,:), done(:,:,:), dplain(:,:,:)
        real    :: r, cerr, derr
        integer :: i, j, k, iw
        write(*,'(A)') 'test_weighted_insertion_linearity'
        params%box        = BOX
        params%box_crop   = BOX
        params%box_croppd = OSMPL_PAD_FAC * BOX
        params%smpd_crop  = SMPD
        params%nstates    = 1
        params%numlen     = 1
        params%oritype    = 'cls3D'
        call poses%new(NPLANES, .false.)
        call poses%spiral()
        call obs%new([BOX,BOX,1], SMPD, wthreads=.false.)
        call obs_pad%new([OSMPL_PAD_FAC*BOX,OSMPL_PAD_FAC*BOX,1], SMPD, wthreads=.false.)
        call grid%new(OSMPL_PAD_FAC*BOX)
        call o%new(.false.)
        call c1sym%new('c1')
        call memoize_ft_maps([OSMPL_PAD_FAC*BOX,OSMPL_PAD_FAC*BOX,1], SMPD)
        ctfparms%smpd    = SMPD
        ctfparms%ctfflag = CTFFLAG_NO
        call set_fixed_seed(SEED)
        allocate(img(BOX,BOX,1))
        do iw = 1, size(WEIGHTS)
            call rw%new_accumulator(params, project, expand=.true., wthreads=.false.)
            call rs%new_accumulator(params, project, expand=.true., wthreads=.false.)
            call rone%new_accumulator(params, project, expand=.true., wthreads=.false.)
            call rplain%new_accumulator(params, project, expand=.true., wthreads=.false.)
            do i = 1, NPLANES
                do k = 1, BOX
                    do j = 1, BOX
                        call random_number(r)
                        img(j,k,1) = r - 0.5
                    enddo
                enddo
                call poses%get_ori(i, o)
                call obs%set_rmat(img, .false.)
                call obs%pad(obs_pad, backgr=0., antialiasing=.false.)
                call obs_pad%fft()
                call obs_pad%gen_fplane4rec(SMPD, ctfparms, [0.,0.], fplane, grid)
                call rw%insert_plane_oversamp(c1sym, o, fplane, w=WEIGHTS(iw))
                call rs%insert_plane_oversamp(c1sym, o, fplane)
                call rone%insert_plane_oversamp(c1sym, o, fplane, w=1.0)
                call rplain%insert_plane_oversamp(c1sym, o, fplane)
            enddo
            ! one plane set per weight: scaling the unweighted sums by w must equal weighting each plane
            call rs%apply_weight(WEIGHTS(iw))
            call rw%compress_exp()
            call rs%compress_exp()
            call rone%compress_exp()
            call rplain%compress_exp()
            cw = rw%get_cmat()
            cs = rs%get_cmat()
            call rw%get_rho_copy(dw)
            call rs%get_rho_copy(ds)
            cerr = maxval(abs(cw - cs)) / max(maxval(abs(cs)), tiny(1.))
            derr = maxval(abs(dw - ds)) / max(maxval(abs(ds)), tiny(1.))
            call assert_true(cerr < RELTOL .and. derr < RELTOL, &
                &'insertion with weight w equals weight one scaled by w (data and density)')
            cone   = rone%get_cmat()
            cplain = rplain%get_cmat()
            call rone%get_rho_copy(done)
            call rplain%get_rho_copy(dplain)
            call assert_true(all(cone == cplain) .and. all(done == dplain), &
                &'insertion with weight one equals the unweighted insertion bit for bit')
            call rw%kill
            call rs%kill
            call rone%kill
            call rplain%kill
        enddo
        call forget_ft_maps()
        call c1sym%kill()
        call o%kill()
        call poses%kill()
        call obs%kill()
        call obs_pad%kill()
    end subroutine test_weighted_insertion_linearity

end module simple_reconstructor_tester
