!@descr: flex_pca sparse sites: greedy site pursuit with held-out stops, the support components and the site files
module simple_flex_sites_pursuit
!$ use omp_lib, only: omp_get_max_threads
use simple_core_module_api,          only: dp, fclose, fopen, int2str, int2str_pad, logfhandle, simple_exception, simple_mkdir, string
use simple_defs,                     only: CWD_GLOB
use simple_fileio,                   only: nlines, read_exit_code, simple_copy_file
use simple_syslib,                   only: del_file, dir_exists, file_exists, simple_abspath, simple_chdir, simple_getcwd, simple_rename
use simple_builder,                  only: builder
use simple_cmdline,                  only: cmdline
use simple_parameters,               only: parameters
use simple_qsys_env,                 only: qsys_env
use simple_image,                    only: image
use simple_rnd,                      only: gasdev, get_rnd_state, put_rnd_state
use simple_flex_pca_planes,          only: flex_plane_store
use simple_flex_pca_project_gateway, only: read_pind_list
use simple_flex_sites_fit,           only: flex_sites_ops, flex_sites_grid, flex_sites_state, sites_em_iter, sites_gains, sites_posterior_all, &
    &group_norm2, sites_diag_precision, sites_fit_nuisance, modes_ft, modes_real
implicit none
private
#include "simple_local_flags.inc"

public :: run_sparse_sites

integer,  parameter :: GSZ = 3, SEED_R = 5, MAXSEEDS = 60, NCTL = 12, ENV_DIL = 3, CTL_MARGIN = 4, SITE_R = 2
integer,  parameter :: SEED_ITS = 6, GROW_ITS = 6, FINAL_ITS = 8
real(dp), parameter :: ENV_THR = 0.1d0, KAPPA = 3.d0, PI_DP = 3.14159265358979323846d0
character(len=*), parameter :: SDIR = 'sites', CAP_DIR = 'sites_cap', CAP_FILE = 'sites_cap.txt', CAP_RES = 'sites_cap_gain.txt'
character(len=*), parameter :: CAP_EXIT = 'sites_cap.exit', CAP_PROJ = 'sites_cap.simple', CAP_STOP = 'sites_cap.stop'
integer,          parameter :: CAP_POLL = 5

contains

    subroutine run_sparse_sites( params, build, plane_store, pinds, nptcls, cline )
        class(parameters),       intent(in)    :: params
        class(builder),          intent(inout) :: build
        class(flex_plane_store), intent(inout) :: plane_store
        integer,                 intent(in)    :: pinds(:), nptcls
        class(cmdline),          intent(in)    :: cline
        type(flex_sites_ops), target  :: full, sub
        type(flex_sites_ops), pointer :: ops
        type(flex_sites_grid)    :: grd
        type(flex_sites_state)   :: t, r, tbest
        complex(dp), allocatable :: wn(:,:), wseed(:,:), wc(:,:), wcb(:,:)
        real(dp),    allocatable :: cent(:,:), x0(:,:), gn(:,:), dpm(:), gain(:), se(:), xnew(:,:), m(:,:), sdv(:,:), ll(:)
        real(dp),    allocatable :: lgain(:), lse(:), rec(:,:), en(:), mp(:,:), sdp(:,:), llp(:), recb(:,:)
        logical,     allocatable :: env(:), act(:,:), even(:), odd(:), allp(:), sv(:), dead(:), used(:), inball(:)
        integer,     allocatable :: sblk(:), comp(:), cb(:), rs0(:), rsb(:), krng(:), rsc(:), sblk1(:)
        complex(dp), allocatable :: wseed1(:,:), wn1(:,:)
        real(dp),    allocatable :: gn1(:,:)
        character(len=16), allocatable :: tnames(:)
        real(dp)   :: g, sd, thr, base, rnorm, cg, cse, cthr, gj, gbest
        integer    :: nseed, nctl, nb, j, k, rnd, best, gi, newb, ib, nrec, ntot, nrecb, kbj, kbest, jl, ncap, nrs0, nrsb
        integer    :: kcap, nrsc, nbc
        integer(8) :: tc0, tc1, tc2, trate
        logical    :: ok, okc, l_guard, l_capjob, l_subbest, l_r1fail
        if( file_exists(string(CAP_FILE)) .and. file_exists(string(CAP_STOP)) )then
            write(logfhandle,'(A)') '>>> SPARSE_SITES band ladder: the node''s own job ran this cap; nothing to do'
            return
        endif
        call full%new(params, build, plane_store, pinds, nptcls)
        call plane_store%release
        ops => full
        call make_grid(ops, grd)
        allocate(allp(ops%n), source=.true.)
        even = ops%eo == 0
        odd  = .not. even
        call seed_centres(ops, grd, cent, env, nseed)
        call control_centres(ops, grd, env, cent, nctl)
        nb = nseed + nctl
        allocate(act(nb,grd%ng), source=.false.)
        do j = 1, nb
            act(j,:) = ball_groups(ops, grd, cent(:,j))
        end do
        write(logfhandle,'(A,I0,A,I0,A,I0,A,I0,A)') '>>> SPARSE_SITES ', nseed, ' seeds and ', nctl, ' solvent controls, ', &
            &minval(count(act, 2)), '-', maxval(count(act, 2)), ' groups per ball'
        allocate(sv(ops%n3), source=.false.)
        sv(grd%sv) = .true.
        kbest     = 0
        gbest     = 0.d0
        ncap      = 0
        l_subbest = .false.
        l_r1fail  = .false.
        ! quarter-octave caps from the node's band down one octave, searched coarser until the pursued held-out gain drops; every
        ! cap draws what a run at its band alone would (cap 2 in its own job), and a seedless node band screens the coarser caps
        call get_rnd_state(rs0, nrs0)
        allocate(krng(0))
        do jl = 0, 4
            kbj = nint(real(full%kb, dp) * 2.d0**(-0.25d0 * real(jl, dp)))
            if( .not. any(krng == kbj) ) krng = [krng, kbj]
        end do
        l_capjob = file_exists(string(CAP_FILE))
        if( l_capjob )then
            call read_cap
        else if( size(krng) > 1 )then
            call launch_cap(krng(2))
        endif
        do jl = 1, size(krng)
            kbj = krng(jl)
            if( l_capjob .and. kbj /= kcap ) cycle
            call system_clock(tc0, trate)
            if( .not. l_capjob .and. jl == 2 )then
                if( .not. file_exists(string(CAP_DIR//'/'//CAP_EXIT//'.job')) ) call stop_cap
                if( wait_cap(gj, nbc) )then
                    call system_clock(tc2)
                    ncap = ncap + 1
                    call ladder_line(kbj, nbc, gj, real(tc2 - tc0) / real(trate), 0., ' (its own job)')
                    if( nbc > 0 .and. (kbest == 0 .or. gj > gbest) )then
                        kbest     = kbj
                        gbest     = gj
                        l_subbest = .true.
                    else if( kbest > 0 )then
                        exit
                    endif
                    cycle
                endif
            endif
            if( jl > 1 )then
                call sub%band(full, kbj)
                ops => sub
            endif
            if( .not. l_capjob .and. l_r1fail .and. kbest == 0 .and. jl > 2 )then
                if( .not. screen_cap() )then
                    call system_clock(tc2)
                    ncap = ncap + 1
                    call ladder_line(kbj, 0, 0.d0, real(tc2 - tc0) / real(trate), 0., ' (screened out)')
                    cycle
                endif
            endif
            if( l_capjob )then
                call put_rnd_state(rsc, nrsc)
            else
                call put_rnd_state(rs0, nrs0)
            endif
            call pursue(gj, tc1)
            call system_clock(tc2)
            ncap = ncap + 1
            call ladder_line(kbj, t%nb, gj, real(tc2 - tc0) / real(trate), real(tc1 - tc0) / real(trate), '')
            if( l_capjob ) call write_cap_gain(kbj, t%nb, gj)
            if( ncap == 1 )then
                recb  = rec(:,:nrec)
                nrecb = nrec
                if( .not. l_capjob .and. nrec > 0 ) l_r1fail = rec(9,1) < 0.5d0
                if( l_r1fail )then
                    wseed1 = wseed
                    sblk1  = sblk
                    wn1    = wn
                    gn1    = gn
                endif
            endif
            if( t%nb > 0 .and. (kbest == 0 .or. gj > gbest) )then
                kbest     = kbj
                gbest     = gj
                tbest     = t
                recb      = rec(:,:nrec)
                nrecb     = nrec
                wcb       = wseed(:,nseed+1:nb)
                l_subbest = .false.
                call get_rnd_state(rsb, nrsb)
            else if( kbest > 0 )then
                exit
            endif
        end do
        if( l_subbest )then
            if( adopt_cap() )then
                write(logfhandle,'(A,I0,A,F6.1,A,I0,A)') '>>> SPARSE_SITES band: ring ', kbest, ' (', real(ops%box) * ops%smpd / &
                    &real(kbest), ' A), the largest pursued held-out gain of ', ncap, ' caps; its own job''s sites are kept'
                write(logfhandle,'(A,I0,A,I0,A)') '>>> SPARSE_SITES done: ', nlines(string(SDIR//'/blocks.txt')), &
                    &' site blocks in ', nlines(string(SDIR//'/sites.txt')), ' sites'
                call flush(logfhandle)
                call full%kill
                call sub%kill
                return
            endif
            call sub%band(full, kbest)
            ops => sub
            call put_rnd_state(rs0, nrs0)
            call pursue(gj, tc1)
            tbest = t
            recb  = rec(:,:nrec)
            nrecb = nrec
            wcb   = wseed(:,nseed+1:nb)
            call get_rnd_state(rsb, nrsb)
        endif
        call simple_mkdir(SDIR)
        call write_rounds(recb, nrecb)
        if( kbest > 0 )then
            if( kbest /= ops%kb )then
                if( kbest == full%kb )then
                    call sub%kill
                    ops => full
                else
                    call sub%band(full, kbest)
                    ops => sub
                endif
            endif
            if( associated(ops, sub) )then
                call full%kill
            else
                deallocate(full%rs, full%rc, full%rt, full%ry)
            endif
            t = tbest
            wseed(:,nseed+1:nb) = wcb
            write(logfhandle,'(A,I0,A,F6.1,A,I0,A)') '>>> SPARSE_SITES band: ring ', kbest, ' (', real(ops%box) * ops%smpd / &
                &real(kbest), ' A), the largest pursued held-out gain of ', ncap, ' caps'
            call put_rnd_state(rsb, nrsb)
            call leave_one_out(t, lgain, lse)
        endif
        if( t%nb == 0 )then
            write(logfhandle,'(A)') '>>> SPARSE_SITES no site survives'
            call write_empty
            call full%kill
            call sub%kill
            return
        endif
        do j = 1, FINAL_ITS
            call sites_em_iter(ops, grd, t, allp, KAPPA, .true., 'final it '//int2str(j))
        end do
        call orient_blocks(t)
        comp = components(t)
        call sites_posterior_all(ops, t, m, sdv, ll, mp, sdp, llp)
        call write_sites(t, comp, lgain, lse)
        call write_coords(t, SDIR//'/coords.txt', m, sdv, ll)
        call write_coords(t, SDIR//'/coords_profiled.txt', mp, sdp, llp)
        write(logfhandle,'(A,I0,A,I0,A)') '>>> SPARSE_SITES done: ', t%nb, ' site blocks in ', maxval(comp), ' sites'
        call flush(logfhandle)
        call full%kill
        call sub%kill

      contains

        subroutine ladder_line( kk, nn, gg, secs, sph, tag )
            integer,          intent(in) :: kk, nn
            real(dp),         intent(in) :: gg
            real,             intent(in) :: secs, sph
            character(len=*), intent(in) :: tag
            write(logfhandle,'(A,I0,A,F6.1,A,I0,A,F12.1,A,F9.1,A,F9.1,A,A)') '>>> SPARSE_SITES band ladder: ring ', kk, ' (', &
                &real(ops%box) * ops%smpd / real(kk), ' A), ', nn, ' site blocks, pursued held-out gain ', gg, ' nats; ', &
                &secs, ' s, phase 1 ', sph, ' s', tag
            call flush(logfhandle)
        end subroutine ladder_line

        subroutine read_cap
            integer :: fu, ios, nsz
            call fopen(fu, string(CAP_FILE), status='old', action='read', iostat=ios)
            read(fu,*) kcap, nsz, nrsc
            allocate(rsc(nsz))
            read(fu,*) rsc
            call fclose(fu)
            write(logfhandle,'(A,I0,A)') '>>> SPARSE_SITES band ladder: this job runs ring ', kcap, ' of its node''s ladder'
            call flush(logfhandle)
        end subroutine read_cap

        subroutine launch_cap( kk )
            integer, intent(in) :: kk
            character(len=11), parameter :: PKEYS(5) = [character(len=11) :: 'vol1', 'rowmap', 'pindfile', 'mskfile', 'pcg_mskfile']
            type(qsys_env) :: q
            type(cmdline)  :: cl
            type(string)   :: root, here, ap, xp
            character(len=4096) :: exe
            real(dp)       :: mb
            integer        :: fu, ios, ik
            call simple_mkdir(CAP_DIR)
            call del_file(string(CAP_DIR//'/'//CAP_RES))
            call del_file(string(CAP_DIR//'/'//CAP_EXIT))
            call del_file(string(CAP_DIR//'/'//CAP_EXIT//'.job'))
            call del_file(string(CAP_DIR//'/'//CAP_STOP))
            call fopen(fu, string(CAP_DIR//'/'//CAP_FILE), status='replace', action='write', iostat=ios)
            write(fu,*) kk, size(rs0), nrs0
            write(fu,*) rs0
            call fclose(fu)
            ap = simple_abspath(params%projfile)
            call simple_copy_file(ap, string(CAP_DIR//'/'//CAP_PROJ))
            cl = cline
            call cl%set('projfile', CAP_PROJ)
            call cl%set('mkdir',    'no')
            do ik = 1, size(PKEYS)
                if( .not. cl%defined(trim(PKEYS(ik))) ) cycle
                ap = cl%get_carg(trim(PKEYS(ik)))
                ap = simple_abspath(ap)
                call cl%set(trim(PKEYS(ik)), ap%to_char())
            end do
            mb = (real(ops%nc, dp) * real(nb, dp) * real(nb + 1, dp) * 8.d0 + real(ops%ip_p(ops%n+1) - 1, dp) * 32.d0) / 1.d6 + 4000.d0
            if( cline%defined('qsys_partition') )then
                call q%new(params, 1, exec_bin=string('simple_exec'), qsys_nthr=params%nthr, &
                    &qsys_partition=cline%get_carg('qsys_partition'))
            else
                call q%new(params, 1, exec_bin=string('simple_exec'), qsys_nthr=params%nthr)
            endif
            call q%qdescr%set('job_memory_per_task', int2str(nint(mb)))
            call q%qdescr%set('job_cpus_per_task',   int2str(params%nthr))
            call get_command_argument(0, exe)
            xp = string('')
            if( index(exe, '/') > 0 ) xp = simple_abspath(trim(exe))
            call simple_getcwd(root)
            call simple_chdir(string(CAP_DIR))
            call simple_getcwd(here)
            CWD_GLOB = here%to_char()
            if( xp%strlen_trim() > 0 )then
                call q%exec_simple_prg_in_queue_async(cl, string('./sites_cap.sh'), string('sites_cap.log'), exec_bin=xp, &
                    &exit_code_fname=string(CAP_EXIT))
            else
                call q%exec_simple_prg_in_queue_async(cl, string('./sites_cap.sh'), string('sites_cap.log'), &
                    &exit_code_fname=string(CAP_EXIT))
            endif
            call simple_chdir(root)
            CWD_GLOB = root%to_char()
            call q%kill
            write(logfhandle,'(A,I0,A,A,A,I0,A)') '>>> SPARSE_SITES band ladder: ring ', kk, ' submitted as its own job in ', CAP_DIR, &
                &' (', nint(mb), ' MB)'
            call flush(logfhandle)
        end subroutine launch_cap

        subroutine stop_cap
            integer :: fu, ios
            call fopen(fu, string(CAP_DIR//'/'//CAP_STOP), status='replace', action='write', iostat=ios)
            write(fu,'(A)') 'the node''s own job runs this cap'
            call fclose(fu)
            write(logfhandle,'(A)') '>>> SPARSE_SITES band ladder: the cap job has not started; the cap runs here'
            call flush(logfhandle)
        end subroutine stop_cap

        logical function wait_cap( gg, nn )
            real(dp), intent(out) :: gg
            integer,  intent(out) :: nn
            integer :: fu, ios, kk, code
            logical :: err
            wait_cap = .false.
            gg = 0.d0
            nn = 0
            if( file_exists(string(CAP_DIR//'/'//CAP_STOP)) ) return
            do
                if( file_exists(string(CAP_DIR//'/'//CAP_RES)) )then
                    call fopen(fu, string(CAP_DIR//'/'//CAP_RES), status='old', action='read', iostat=ios)
                    read(fu,*,iostat=ios) kk, nn, gg
                    call fclose(fu)
                    wait_cap = ios == 0
                    return
                endif
                if( file_exists(string(CAP_DIR//'/'//CAP_EXIT)) )then
                    call read_exit_code(string(CAP_DIR//'/'//CAP_EXIT), code, err)
                    if( .not. err )then
                        write(logfhandle,'(A,I0,A)') '>>> SPARSE_SITES band ladder: the cap job ended (exit ', code, &
                            &') without its gain; the cap runs here'
                        return
                    endif
                endif
                call sleep(CAP_POLL)
            end do
        end function wait_cap

        subroutine write_cap_gain( kk, nn, gg )
            integer,  intent(in) :: kk, nn
            real(dp), intent(in) :: gg
            integer :: fu, ios
            call fopen(fu, string(CAP_RES//'.tmp'), status='replace', action='write', iostat=ios)
            write(fu,*) kk, nn, gg
            call fclose(fu)
            call simple_rename(string(CAP_RES//'.tmp'), string(CAP_RES))
        end subroutine write_cap_gain

        logical function adopt_cap()
            integer :: code
            logical :: err
            do while( .not. file_exists(string(CAP_DIR//'/'//CAP_EXIT)) )
                call sleep(CAP_POLL)
            end do
            call read_exit_code(string(CAP_DIR//'/'//CAP_EXIT), code, err)
            adopt_cap = .not. err .and. code == 0 .and. file_exists(string(CAP_DIR//'/'//SDIR//'/sites.txt')) .and. &
                &.not. dir_exists(string(SDIR))
            if( adopt_cap )then
                call simple_rename(string(CAP_DIR//'/'//SDIR), string(SDIR))
                adopt_cap = file_exists(string(SDIR//'/sites.txt'))
            endif
            if( .not. adopt_cap ) write(logfhandle,'(A)') '>>> SPARSE_SITES band ladder: the kept cap''s sites could not be '//&
                &'taken from its job; the cap runs here'
        end function adopt_cap

        logical function screen_cap()
            type(flex_sites_state) :: t0
            real(dp), allocatable  :: gs(:), ss(:)
            real(dp) :: bs
            integer  :: jb
            call init_state(t0, x0(:,:0), [integer ::], 0, wn1, gn1)
            allocate(gs(nb), ss(nb))
            call sites_gains(ops, t0, wseed1, sblk1, nb, odd, gs, ss, bs)
            jb = maxloc(gs(:nseed), 1)
            screen_cap = gs(jb) > 2.d0 * ss(jb) .and. gs(jb) > maxval(gs(nseed+1:nb)) .and. gs(jb) > 0.d0
        end function screen_cap

        subroutine pursue( gp, tph )
            real(dp),   intent(out) :: gp
            integer(8), intent(out) :: tph
            type(flex_sites_state) :: s1
            real(dp) :: base0
            if( allocated(x0) ) deallocate(x0, wc, wseed, rec, used, gain, se, gn)
            call templates(ops, wn, tnames)
            call sites_fit_nuisance(ops, wn, allp, gn)
            allocate(x0(ops%n3,nb), source=0.d0)
            do j = 1, nb
                do k = 1, grd%nsupp
                    if( act(j,grd%gid(k)) ) x0(grd%sv(k),j) = real(gasdev(), dp)
                end do
            end do
            call bandlimit(ops, x0)
            do j = 1, nb
                where( .not. sv ) x0(:,j) = 0.d0
            end do
            allocate(wc(ops%nc,nb))
            call modes_ft(ops, x0, wc)
            dpm = sites_diag_precision(ops, wc)
            do j = 1, nb
                x0(:,j) = x0(:,j) * sqrt(1.d0 / max(dpm(j), 1.d-300))
            end do
            call init_state(s1, x0, [(j, j=1,nb)], nb, wn, gn)
            do j = 1, SEED_ITS
                call sites_em_iter(ops, grd, s1, even, KAPPA, .false., 'seeds it '//int2str(j), bsupp=act)
            end do
            call system_clock(tph)
            sblk = s1%blk
            allocate(wseed(ops%nc,nb))
            call modes_ft(ops, s1%x, wseed)
            call init_state(t, s1%x(:,:0), [integer ::], 0, wn, s1%gn)
            t%alpha = s1%alpha; t%rho = s1%rho; t%rhor = s1%rhor
            allocate(rec(13,nseed), source=0.d0)
            allocate(used(nseed), source=.false.)
            nrec = 0
            l_guard = .true.
            allocate(gain(nb), se(nb))
            do rnd = 1, nseed
                call sites_gains(ops, t, wseed, sblk, nb, odd, gain, se, base)
                if( rnd == 1 ) base0 = base
                thr  = maxval(gain(nseed+1:nb))
                best = 0
                do j = 1, nseed
                    if( used(j) ) cycle
                    if( best == 0 )then
                        best = j
                    else if( gain(j) > gain(best) )then
                        best = j
                    endif
                end do
                if( best == 0 ) exit
                g  = gain(best); sd = se(best)
                ok = g > 2.d0 * sd .and. g > thr .and. g > 0.d0
                nrec = nrec + 1
                rec(:,nrec) = [real(rnd,dp), real(best,dp), cent(:,best), g, sd, thr, merge(1.d0, 0.d0, ok), 0.d0, 0.d0, 0.d0, -1.d0]
                write(logfhandle,'(A,I0,A,F10.5,A,I0,A,3F6.1,A,F10.1,A,F7.1,A,F7.1,A,A)') '>>> SPARSE_SITES pursuit round ', rnd, &
                    &' (held-out LL/ptcl of the sites ', base, '): best seed ', best, ' at', cent(:,best), ' gain ', g, ' nats (SE ', &
                    &sd, '); controls <= ', thr, ' -> ', merge('NEW SITE', 'stop    ', ok)
                call flush(logfhandle)
                if( .not. ok )then
                    l_guard = .false.
                    exit
                endif
                used(best) = .true.
                allocate(xnew(ops%n3,SITE_R), source=0.d0)
                xnew(:,1) = s1%x(:,best)
                do j = 2, SITE_R
                    do k = 1, grd%nsupp
                        if( act(best,grd%gid(k)) ) xnew(grd%sv(k),j) = real(gasdev(), dp)
                    end do
                end do
                call bandlimit(ops, xnew(:,2:))
                allocate(inball(ops%n3), source=.false.)
                do k = 1, grd%nsupp
                    inball(grd%sv(k)) = act(best,grd%gid(k))
                end do
                do j = 2, SITE_R
                    where( .not. inball ) xnew(:,j) = 0.d0
                    rnorm = sqrt(sum(xnew(:,j)**2))
                    xnew(:,j) = xnew(:,j) * 0.3d0 * sqrt(sum(xnew(:,1)**2)) / max(rnorm, 1.d-300)
                end do
                call add_block(t, xnew)
                deallocate(xnew, inball)
                newb = t%nb
                do gi = 1, GROW_ITS
                    call sites_em_iter(ops, grd, t, even, KAPPA, .true., 'round '//int2str(rnd)//' grow '//int2str(gi))
                end do
                en = block_energy(t)
                dead = en <= 0.d0
                if( dead(newb) )then
                    newb = 0
                else
                    newb = newb - count(dead(:newb-1))
                endif
                do ib = t%nb, 1, -1
                    if( dead(ib) ) call drop_block(t, ib)
                end do
                okc = newb > 0
                cg = 0.d0; cse = 0.d0; cthr = 0.d0
                if( okc )then
                    call block_gain(t, newb, cg, cse)
                    cthr = control_max(t)
                    okc  = cg > 2.d0 * cse .and. cg > cthr .and. cg >= g
                endif
                rec(10:13,nrec) = [cg, cse, cthr, merge(1.d0, 0.d0, okc)]
                write(logfhandle,'(A,I0,A,F10.1,A,F7.1,A,F7.1,A,F10.1,A,A)') '>>> SPARSE_SITES round ', rnd, ' check: the grown block adds ', &
                    &cg, ' nats held out given the other sites (SE ', cse, '); controls <= ', cthr, '; seed ', g, ' -> ', &
                    &merge('kept   ', 'dropped', okc)
                call flush(logfhandle)
                if( .not. okc .and. newb > 0 ) call drop_block(t, newb)
            end do
            if( l_guard )then
                write(logfhandle,'(A)') '>>> SPARSE_SITES runaway guard: every seed was used before the stop'
                call sites_gains(ops, t, wseed, sblk, nb, odd, gain, se, base)
            endif
            gp = 0.d0
            if( t%nb > 0 ) gp = (base - base0) * real(count(odd), dp)
        end subroutine pursue

        subroutine block_gain( st, bb, gg, ss )
            type(flex_sites_state), intent(in)  :: st
            integer,                intent(in)  :: bb
            real(dp),               intent(out) :: gg, ss
            type(flex_sites_state)   :: rr
            complex(dp), allocatable :: wb(:,:)
            integer,     allocatable :: jb(:)
            real(dp) :: g1(1), s1v(1), bl
            jb = pack([(j, j=1,st%qs)], st%blk == bb)
            allocate(wb(ops%nc,size(jb)))
            call modes_ft(ops, st%x(:,jb), wb)
            rr = st
            call drop_block(rr, bb)
            call sites_gains(ops, rr, wb, [(1, j=1,size(jb))], 1, odd, g1, s1v, bl)
            gg = g1(1); ss = s1v(1)
        end subroutine block_gain

        real(dp) function control_max( st )
            type(flex_sites_state), intent(in) :: st
            real(dp) :: gc(nctl), sc(nctl), bl
            call sites_gains(ops, st, wseed(:,nseed+1:nb), [(j, j=1,nctl)], nctl, odd, gc, sc, bl)
            control_max = maxval(gc)
        end function control_max

        subroutine leave_one_out( st, lg, ls )
            type(flex_sites_state), intent(inout) :: st
            real(dp), allocatable,  intent(inout) :: lg(:), ls(:)
            integer  :: w
            real(dp) :: cmx
            do while( st%nb >= 1 )
                if( allocated(lg) ) deallocate(lg, ls)
                allocate(lg(st%nb), ls(st%nb))
                do ib = 1, st%nb
                    call block_gain(st, ib, lg(ib), ls(ib))
                end do
                cmx = control_max(st)
                write(logfhandle,'(A,*(1X,I0,A,F0.1,A,F0.1,A))') '>>> SPARSE_SITES leave-one-block-out (held-out):', &
                    &(ib, ':', lg(ib), '(', ls(ib), ')', ib=1,st%nb)
                write(logfhandle,'(A,F0.1)') '>>> SPARSE_SITES    controls <= ', cmx
                w = minloc(lg, 1)
                if( lg(w) > 2.d0 * ls(w) .and. lg(w) > cmx ) exit
                write(logfhandle,'(A,I0,A,F0.1,A,F0.1,A)') '>>> SPARSE_SITES    dropped redundant site block ', w, ' (gain ', lg(w), &
                    &', SE ', ls(w), ')'
                call drop_block(st, w)
            end do
            call flush(logfhandle)
        end subroutine leave_one_out

        function block_energy( st ) result( e )
            type(flex_sites_state), intent(in) :: st
            real(dp), allocatable :: e(:)
            integer :: jj
            allocate(e(st%nb), source=0.d0)
            do jj = 1, st%qs
                e(st%blk(jj)) = e(st%blk(jj)) + sum(st%x(:,jj)**2)
            end do
        end function block_energy

        subroutine orient_blocks( st )
            type(flex_sites_state), intent(inout) :: st
            real(dp), allocatable :: c(:,:), ev(:), work(:), v(:,:)
            integer,  allocatable :: jb(:)
            integer :: nr, info
            do ib = 1, st%nb
                jb = pack([(j, j=1,st%qs)], st%blk == ib)
                nr = size(jb)
                if( nr < 2 ) cycle
                c = matmul(transpose(st%x(grd%sv,jb)), st%x(grd%sv,jb))
                allocate(ev(nr), work(3*nr))
                call dsyev('V', 'U', nr, c, nr, ev, work, 3*nr, info)
                v = c(:,nr:1:-1)
                st%x(:,jb) = matmul(st%x(:,jb), v)
                deallocate(ev, work)
            end do
        end subroutine orient_blocks

        function components( st ) result( lab )
            type(flex_sites_state), intent(in) :: st
            integer, allocatable :: lab(:), root(:)
            logical, allocatable :: sup(:,:)
            real(dp) :: kk, mm, mu, sdh, nv
            integer  :: a, b, inter, ra, rb, nsite
            allocate(sup(grd%nsupp,st%nb), source=.false.)
            do j = 1, st%qs
                sup(:,st%blk(j)) = sup(:,st%blk(j)) .or. st%x(grd%sv,j) /= 0.d0
            end do
            nv = real(grd%nsupp, dp)
            allocate(root(st%nb))
            root = [(a, a=1,st%nb)]
            do a = 1, st%nb
                do b = a + 1, st%nb
                    inter = count(sup(:,a) .and. sup(:,b))
                    kk  = real(count(sup(:,a)), dp)
                    mm  = real(count(sup(:,b)), dp)
                    mu  = kk * mm / nv
                    sdh = sqrt(mm * (kk / nv) * (1.d0 - kk / nv) * (nv - mm) / (nv - 1.d0))
                    if( real(inter, dp) > mu + 2.d0 * sdh )then
                        ra = a
                        do while( root(ra) /= ra )
                            ra = root(ra)
                        end do
                        rb = b
                        do while( root(rb) /= rb )
                            rb = root(rb)
                        end do
                        if( ra /= rb ) root(ra) = rb
                    endif
                end do
            end do
            allocate(lab(st%nb), source=0)
            nsite = 0
            do a = 1, st%nb
                ra = a
                do while( root(ra) /= ra )
                    ra = root(ra)
                end do
                if( lab(ra) == 0 )then
                    nsite   = nsite + 1
                    lab(ra) = nsite
                endif
                lab(a) = lab(ra)
            end do
        end function components

        subroutine write_rounds( rr, nn )
            real(dp), intent(in) :: rr(:,:)
            integer,  intent(in) :: nn
            integer :: fu, ios, kk
            call fopen(fu, string(SDIR//'/rounds.txt'), status='replace', action='write', iostat=ios)
            write(fu,'(A)') '# round seed x y z gain se controls_max accepted check_gain check_se check_controls check_ok'
            do kk = 1, nn
                write(fu,'(2I5,3F7.1,3F12.2,I3,3F12.2,I3)') nint(rr(1,kk)), nint(rr(2,kk)), rr(3:5,kk), rr(6:8,kk), nint(rr(9,kk)), &
                    &rr(10:12,kk), nint(rr(13,kk))
            end do
            call fclose(fu)
        end subroutine write_rounds

        subroutine write_empty
            integer :: fu, ios
            call fopen(fu, string(SDIR//'/sites.txt'), status='replace', action='write', iostat=ios)
            write(fu,'(A)') '# site blocks support_voxels centre_radius_A'
            call fclose(fu)
        end subroutine write_empty

        subroutine write_coords( st, fname, mm, sdd, lll )
            type(flex_sites_state), intent(in) :: st
            character(len=*),       intent(in) :: fname
            real(dp),               intent(in) :: mm(:,:), sdd(:,:), lll(:)
            integer, allocatable :: prow(:), keys(:)
            integer :: fu, ios, i, jj, nk
            allocate(prow(ops%n))
            prow = ops%pind
            if( params%rowmap%strlen() > 0 )then
                call read_pind_list(params%rowmap%to_char(), keys, nk)
                prow = keys(ops%pind)
            endif
            call fopen(fu, string(fname), status='replace', action='write', iostat=ios)
            write(fu,'(A,*(A))') '# row eo fit ll', (' s'//int2str(jj), jj=1,st%qs), (' '//trim(tnames(jj)), jj=1,size(tnames)), &
                &(' sd_s'//int2str(jj), jj=1,st%qs)
            do i = 1, ops%n
                write(fu,'(I0,1X,I0,A,ES14.6,*(1X,ES13.5))') prow(i), ops%eo(i), ' 1 ', lll(i), mm(:,i), sdd(:,i)
            end do
            call fclose(fu)
        end subroutine write_coords

        subroutine write_sites( st, lab, lg, ls )
            type(flex_sites_state), intent(in) :: st
            integer,                intent(in) :: lab(:)
            real(dp),               intent(in) :: lg(:), ls(:)
            real(dp), allocatable :: ev(:), sp(:)
            character(len=:), allocatable :: txt
            integer  :: fu, ios, i, s, nsite
            real(dp) :: cr, w
            call fopen(fu, string(SDIR//'/modes.txt'), status='replace', action='write', iostat=ios)
            write(fu,'(A)') '# mode block site'
            do j = 1, st%qs
                write(fu,'(3I6)') j, st%blk(j), lab(st%blk(j))
                call write_vol(st%x(:,j), SDIR//'/mode_'//int2str_pad(j,3)//'.mrc')
            end do
            call fclose(fu)
            allocate(ev(ops%n3), sp(ops%n3))
            call fopen(fu, string(SDIR//'/blocks.txt'), status='replace', action='write', iostat=ios)
            write(fu,'(A)') '# block site modes support_voxels energy loo_gain loo_se centroid_A'
            do ib = 1, st%nb
                ev = 0.d0
                do j = 1, st%qs
                    if( st%blk(j) == ib ) ev = ev + st%x(:,j)**2
                end do
                sp = merge(1.d0, 0.d0, ev > 0.d0)
                cr = 0.d0
                w  = max(sum(ev), 1.d-300)
                do i = 1, ops%n3
                    if( ev(i) > 0.d0 ) cr = cr + ev(i) * sqrt(real(sum(vox_xyz(ops, i)**2), dp))
                end do
                write(fu,'(I5,I5,I5,I8,ES12.4,2F12.1,F9.1)') ib, lab(ib), count(st%blk == ib), count(ev > 0.d0), sum(ev), lg(ib), ls(ib), &
                    &cr / w * ops%smpd
                call write_vol(ev, SDIR//'/block_'//int2str_pad(ib,2)//'_energy.mrc')
                call write_vol(sp, SDIR//'/block_'//int2str_pad(ib,2)//'_support.mrc')
            end do
            call fclose(fu)
            nsite = maxval(lab)
            call fopen(fu, string(SDIR//'/sites.txt'), status='replace', action='write', iostat=ios)
            write(fu,'(A)') '# site blocks support_voxels centre_radius_A'
            do s = 1, nsite
                sp = 0.d0
                ev = 0.d0
                do j = 1, st%qs
                    if( lab(st%blk(j)) == s ) ev = ev + st%x(:,j)**2
                end do
                sp = merge(1.d0, 0.d0, ev > 0.d0)
                cr = 0.d0
                do i = 1, ops%n3
                    if( ev(i) > 0.d0 ) cr = cr + ev(i) * sqrt(real(sum(vox_xyz(ops, i)**2), dp))
                end do
                txt = ''
                do ib = 1, st%nb
                    if( lab(ib) /= s ) cycle
                    if( len(txt) > 0 ) txt = txt//','
                    txt = txt//int2str(ib)
                end do
                write(fu,'(I5,1X,A,I8,F9.1)') s, txt, count(ev > 0.d0), cr / max(sum(ev), 1.d-300) * ops%smpd
                call write_site_mask(sp, SDIR//'/site_'//int2str_pad(s,2)//'.mrc')
            end do
            call fclose(fu)
            call write_vol(merge(1.d0, 0.d0, sv), SDIR//'/support.mrc')
        end subroutine write_sites

        subroutine write_vol( v, fname )
            real(dp),         intent(in) :: v(:)
            character(len=*), intent(in) :: fname
            type(image) :: img
            call img%new([ops%box,ops%box,ops%box], ops%smpd)
            call img%set_rmat(reshape(real(v), [ops%box,ops%box,ops%box]), .false.)
            call img%write(string(fname))
            call img%kill
        end subroutine write_vol

        subroutine write_site_mask( v, fname )
            real(dp),         intent(in) :: v(:)
            character(len=*), intent(in) :: fname
            type(image) :: img, big
            real, allocatable :: rm(:,:,:)
            call img%new([ops%box,ops%box,ops%box], ops%smpd)
            call img%set_rmat(reshape(real(v), [ops%box,ops%box,ops%box]), .false.)
            call img%fft()
            call big%new([params%box,params%box,params%box], params%smpd)
            call img%pad(big)
            call big%ifft()
            rm = big%get_rmat()
            rm = min(max(rm, 0.), 1.)
            call big%set_rmat(rm, .false.)
            call big%write(string(fname))
            call img%kill
            call big%kill
        end subroutine write_site_mask

    end subroutine run_sparse_sites

    subroutine init_state( st, x, blk, nb, wn, gn )
        type(flex_sites_state), intent(inout) :: st
        real(dp),               intent(in)    :: x(:,:), gn(:,:)
        integer,                intent(in)    :: blk(:), nb
        complex(dp),            intent(in)    :: wn(:,:)
        st%qs  = size(x,2)
        st%nb  = nb
        st%x   = x
        st%z   = x
        st%zr  = x
        st%u   = 0.d0 * x
        st%ur  = 0.d0 * x
        st%blk = blk
        st%wn  = wn
        st%gn  = gn
    end subroutine init_state

    subroutine add_block( st, xnew )
        type(flex_sites_state), intent(inout) :: st
        real(dp),               intent(in)    :: xnew(:,:)
        integer :: r, k
        r = size(xnew,2)
        st%x   = reshape([st%x,  xnew], [size(xnew,1), st%qs + r])
        st%z   = reshape([st%z,  xnew], [size(xnew,1), st%qs + r])
        st%zr  = reshape([st%zr, xnew], [size(xnew,1), st%qs + r])
        st%u   = reshape([st%u,  0.d0 * xnew], [size(xnew,1), st%qs + r])
        st%ur  = reshape([st%ur, 0.d0 * xnew], [size(xnew,1), st%qs + r])
        st%blk = [st%blk, (st%nb + 1, k = 1, r)]
        st%nb  = st%nb + 1
        st%qs  = st%qs + r
    end subroutine add_block

    subroutine drop_block( st, b )
        type(flex_sites_state), intent(inout) :: st
        integer,                intent(in)    :: b
        integer, allocatable :: idx(:)
        integer :: k
        idx    = pack([(k, k=1,st%qs)], st%blk /= b)
        st%x   = st%x(:,idx)
        st%z   = st%z(:,idx)
        st%zr  = st%zr(:,idx)
        st%u   = st%u(:,idx)
        st%ur  = st%ur(:,idx)
        st%blk = st%blk(idx)
        where( st%blk > b ) st%blk = st%blk - 1
        st%qs  = size(idx)
        st%nb  = st%nb - 1
    end subroutine drop_block

    function vox_xyz( ops, v ) result( xyz )
        type(flex_sites_ops), intent(in) :: ops
        integer,              intent(in) :: v
        integer :: xyz(3), c
        c = ops%box / 2 + 1
        xyz(1) = mod(v - 1, ops%box) + 1 - c
        xyz(2) = mod((v - 1) / ops%box, ops%box) + 1 - c
        xyz(3) = (v - 1) / (ops%box * ops%box) + 1 - c
    end function vox_xyz

    subroutine make_grid( ops, grd )
        type(flex_sites_ops),  intent(in)    :: ops
        type(flex_sites_grid), intent(inout) :: grd
        integer, allocatable :: full(:), map(:)
        integer :: v, xyz(3), ngx, o, s
        ngx = (ops%box + GSZ - 1) / GSZ
        o   = ops%box / 2
        allocate(full(ops%n3), source=0)
        grd%nsupp = 0
        do v = 1, ops%n3
            xyz = vox_xyz(ops, v)
            if( real(sum(xyz**2)) > ops%msk**2 ) cycle
            grd%nsupp = grd%nsupp + 1
            full(v) = (((xyz(3) + o) / GSZ) * ngx + (xyz(2) + o) / GSZ) * ngx + (xyz(1) + o) / GSZ + 1
        end do
        allocate(map(ngx**3), source=0)
        do v = 1, ops%n3
            if( full(v) > 0 ) map(full(v)) = 1
        end do
        grd%ng = 0
        do s = 1, size(map)
            if( map(s) > 0 )then
                grd%ng = grd%ng + 1
                map(s) = grd%ng
            endif
        end do
        allocate(grd%sv(grd%nsupp), grd%gid(grd%nsupp), grd%gvox(grd%ng))
        grd%gvox = 0
        s = 0
        do v = 1, ops%n3
            if( full(v) == 0 ) cycle
            s = s + 1
            grd%sv(s)  = v
            grd%gid(s) = map(full(v))
            grd%gvox(grd%gid(s)) = grd%gvox(grd%gid(s)) + 1
        end do
        grd%kappa = KAPPA
        write(logfhandle,'(A,I0,A,F6.3,A,I0,A,I0,A)') '>>> SPARSE_SITES support: ', grd%nsupp, ' voxels (', &
            &real(grd%nsupp) / real(ops%n3), ' of the box), ', grd%ng, ' groups of ', GSZ, '^3 voxels'
    end subroutine make_grid

    subroutine bandlimit( ops, x )
        type(flex_sites_ops), intent(inout) :: ops
        real(dp),             intent(inout) :: x(:,:)
        complex(dp), allocatable :: w(:,:)
        integer :: j
        allocate(w(ops%nc,size(x,2)))
        call modes_ft(ops, x, w)
        do j = 1, size(x,2)
            where( ops%cell_shell > ops%kb ) w(:,j) = cmplx(0.d0, 0.d0, kind=dp)
        end do
        call modes_real(ops, w, x)
    end subroutine bandlimit

    subroutine templates( ops, wn, names )
        type(flex_sites_ops),           intent(inout) :: ops
        complex(dp),       allocatable, intent(out)   :: wn(:,:)
        character(len=16), allocatable, intent(out)   :: names(:)
        complex(dp), allocatable :: dk(:,:)
        real(dp),    allocatable :: d(:,:), t(:,:), xyz(:,:)
        real(dp) :: f(3), ed(5)
        integer  :: c, v, a, j
        allocate(dk(ops%nc,3))
        do c = 1, ops%nc
            do a = 1, 3
                f(a) = real(ops%cell_lg(a,c), dp) / real(ops%box, dp)
                if( abs(ops%cell_lg(a,c)) == ops%nyq ) f(a) = 0.d0
            end do
            dk(c,:) = cmplx(0.d0, 2.d0 * PI_DP, kind=dp) * f * ops%mu(c)
        end do
        allocate(d(ops%n3,3), t(ops%n3,12), xyz(ops%n3,3))
        call modes_real(ops, dk, d)
        do v = 1, ops%n3
            xyz(v,:) = real(vox_xyz(ops, v), dp)
        end do
        t(:,1) = -xyz(:,3) * d(:,2) + xyz(:,2) * d(:,3)
        t(:,2) =  xyz(:,3) * d(:,1) - xyz(:,1) * d(:,3)
        t(:,3) = -xyz(:,2) * d(:,1) + xyz(:,1) * d(:,2)
        t(:,4:6) = d
        t(:,7) = xyz(:,1) * d(:,1) + xyz(:,2) * d(:,2) + xyz(:,3) * d(:,3)
        t(:,8) = merge(1.d0, 0.d0, sum(xyz**2, 2) <= real(ops%msk, dp)**2)
        ed = [0.d0, ops%kb / 4.d0, ops%kb / 2.d0, 3.d0 * ops%kb / 4.d0, ops%kb + 1.d0]
        deallocate(dk)
        allocate(dk(ops%nc,4))
        do j = 1, 4
            where( real(ops%cell_shell, dp) >= ed(j) .and. real(ops%cell_shell, dp) < ed(j+1) )
                dk(:,j) = ops%mu
            elsewhere
                dk(:,j) = cmplx(0.d0, 0.d0, kind=dp)
            end where
        end do
        call modes_real(ops, dk, t(:,9:12))
        call bandlimit(ops, t)
        do j = 1, 12
            t(:,j) = t(:,j) / max(sqrt(sum(t(:,j)**2)), 1.d-300)
        end do
        allocate(wn(ops%nc,12))
        call modes_ft(ops, t, wn)
        names = [character(len=16) :: 'rotx', 'roty', 'rotz', 'trx', 'try', 'trz', 'dil', 'bg', 'mu_b1', 'mu_b2', 'mu_b3', 'mu_b4']
    end subroutine templates

    subroutine seed_centres( ops, grd, cent, env, nseed )
        type(flex_sites_ops),  intent(in)  :: ops
        type(flex_sites_grid), intent(in)  :: grd
        real(dp), allocatable, intent(out) :: cent(:,:)
        logical,  allocatable, intent(out) :: env(:)
        integer,               intent(out) :: nseed
        real(dp), allocatable :: p(:,:), dmin(:)
        logical,  allocatable :: insp(:)
        real(dp) :: c0(3)
        integer  :: v, np, k
        env = ops%mu_real > ENV_THR * maxval(ops%mu_real)
        env = dilate(ops, env, ENV_DIL)
        allocate(insp(ops%n3), source=.false.)
        insp(grd%sv) = .true.
        env = env .and. insp
        np  = count(env)
        allocate(p(3,np))
        k = 0
        do v = 1, ops%n3
            if( .not. env(v) ) cycle
            k = k + 1
            p(:,k) = real(vox_xyz(ops, v), dp)
        end do
        c0 = sum(p, 2) / real(np, dp)
        allocate(cent(3,MAXSEEDS))
        cent(:,1) = p(:,maxloc(sum((p - spread(c0,2,np))**2, 1), 1))
        dmin  = sqrt(sum((p - spread(cent(:,1),2,np))**2, 1))
        nseed = 1
        do while( nseed < MAXSEEDS .and. maxval(dmin) >= real(SEED_R, dp) )
            nseed = nseed + 1
            cent(:,nseed) = p(:,maxloc(dmin, 1))
            dmin = min(dmin, sqrt(sum((p - spread(cent(:,nseed),2,np))**2, 1)))
        end do
        cent = cent(:,:nseed)
        write(logfhandle,'(A,I0,A,I0,A,F5.1,A)') '>>> SPARSE_SITES seeds: ', nseed, ' centres over the envelope (', np, &
            &' voxels); every envelope voxel within ', maxval(dmin), ' voxels of a centre'
    end subroutine seed_centres

    subroutine control_centres( ops, grd, env, cent, nctl )
        type(flex_sites_ops),  intent(in)    :: ops
        type(flex_sites_grid), intent(in)    :: grd
        logical,               intent(in)    :: env(:)
        real(dp), allocatable, intent(inout) :: cent(:,:)
        integer,               intent(out)   :: nctl
        real(dp), allocatable :: p(:,:), dmin(:), cc(:,:)
        logical,  allocatable :: far(:), sol(:)
        integer :: v, np, k
        far = dilate(ops, env, CTL_MARGIN)
        allocate(sol(ops%n3), source=.false.)
        sol(grd%sv) = .true.
        sol = sol .and. .not. far
        np  = count(sol)
        nctl = 0
        if( np == 0 ) return
        allocate(p(3,np))
        k = 0
        do v = 1, ops%n3
            if( .not. sol(v) ) cycle
            k = k + 1
            p(:,k) = real(vox_xyz(ops, v), dp)
        end do
        allocate(cc(3,NCTL))
        cc(:,1) = p(:,1)
        dmin = sqrt(sum((p - spread(cc(:,1),2,np))**2, 1))
        nctl = 1
        do while( nctl < NCTL .and. maxval(dmin) >= real(SEED_R, dp) )
            nctl = nctl + 1
            cc(:,nctl) = p(:,maxloc(dmin, 1))
            dmin = min(dmin, sqrt(sum((p - spread(cc(:,nctl),2,np))**2, 1)))
        end do
        cent = reshape([cent, cc(:,:nctl)], [3, size(cent,2) + nctl])
    end subroutine control_centres

    function dilate( ops, m, nit ) result( d )
        type(flex_sites_ops), intent(in) :: ops
        logical,              intent(in) :: m(:)
        integer,              intent(in) :: nit
        logical, allocatable :: d(:), e(:,:,:), f(:,:,:)
        integer :: it, n
        n = ops%box
        e = reshape(m, [n,n,n])
        do it = 1, nit
            f = e
            f(2:n,:,:) = f(2:n,:,:) .or. e(1:n-1,:,:)
            f(1:n-1,:,:) = f(1:n-1,:,:) .or. e(2:n,:,:)
            f(:,2:n,:) = f(:,2:n,:) .or. e(:,1:n-1,:)
            f(:,1:n-1,:) = f(:,1:n-1,:) .or. e(:,2:n,:)
            f(:,:,2:n) = f(:,:,2:n) .or. e(:,:,1:n-1)
            f(:,:,1:n-1) = f(:,:,1:n-1) .or. e(:,:,2:n)
            e = f
        end do
        d = reshape(e, [n**3])
    end function dilate

    function ball_groups( ops, grd, c ) result( g )
        type(flex_sites_ops),  intent(in) :: ops
        type(flex_sites_grid), intent(in) :: grd
        real(dp),              intent(in) :: c(3)
        logical, allocatable :: g(:)
        integer :: s
        allocate(g(grd%ng), source=.false.)
        do s = 1, grd%nsupp
            if( sum((real(vox_xyz(ops, grd%sv(s)), dp) - c)**2) <= real(SEED_R, dp)**2 ) g(grd%gid(s)) = .true.
        end do
    end function ball_groups

end module simple_flex_sites_pursuit
