!@descr: distributed kernel-latent EM board: the master (prg=klem_board with nparts > 1) holds no particles; NPARTS
!! persistent workers (prg=klem_board_part) each load their share of the particles once (stack order) and answer
!! file requests round by round: consensus sums, calibration sums, M-step partial sums (with their own kernel-smoothed
!! responsibilities), the E-step (responsibilities, per-part balance, optional pose step). The master merges the
!! partial sums, forms the Wiener-shrunk cells, applies focus and deviation weights, and writes the standard board
!! outputs (klem_latent.txt, klem_states.txt, klem_resp.bin, klem_halves.bin, cells/, klem_poses.txt).
!! Exchange files live in DIR/distr: req_NNNN.bin + req_NNNN.ready (master -> workers), rep_PP_NNNN.bin + .done
!! (workers -> master); Fourier cubes travel as the half sphere (h >= 0, |k| <= kmax) and are unpacked by Friedel.
module simple_klem_distr
use simple_core_module_api
use simple_sp_project,  only: sp_project
use simple_parameters,  only: parameters
use simple_cmdline,     only: cmdline
use simple_image,       only: image
use simple_string,      only: string
use simple_qsys_env,    only: qsys_env
use simple_stack_io,    only: stack_io
use simple_srch_sort_loc, only: hpsort
use simple_klem_core
!$ use omp_lib
implicit none
public :: klem_cfg, klem_read_cfg, klem_board_master, klem_board_worker, klem_read_init_positions, apply_focus_d
private
#include "simple_local_flags.inc"

integer, parameter :: REQ_CONS = 1, REQ_SCALE = 2, REQ_NOISE = 3, REQ_MSTEP = 4, REQ_ESTEP = 5, REQ_NEWBOARD = 7, REQ_FLUSH = 8, REQ_FINISH = 9
integer, parameter :: WAIT_START = 21600, WAIT_MAX = 3600   ! the first worker reply may wait for queued workers (6 h); later waits 1 h    ! seconds a worker/master waits for the other side before giving up (a dead master frees its workers within the hour)

type klem_cfg
    integer :: n1 = 7, n2 = 7, nfree = 2, nsub = 0, seed = 1, conv = 0, balance_until = 0, write_halves = 0, pose_every = 0, pnoise = 0, pscale = 0, ctfabs = 1, nthr_part = 0, kmax_load = 0, gridding = 0
    real    :: stop_frac = 0.03, shrink = 0.02, eps0f = 0.05, balance = 0.5, sshift = -1., fscale = 1., csign = 1., pose_step = 3., pose_prior = 4., devw = 0., hp = 0.
    real,    allocatable :: lps(:), sigmas(:)
    character(len=512) :: focus = 'none', subset = 'none', init = 'none', server_list = 'none', solvent = 'none'
end type klem_cfg

contains

    !> the board config ("key value" lines), the same keys as the sequential board plus nthr_part
    subroutine klem_read_cfg( fname, cfg )
        character(len=*), intent(in)    :: fname
        type(klem_cfg),   intent(inout) :: cfg
        character(len=512) :: line, tok(64)
        integer :: u, io, ntok, j
        if( allocated(cfg%lps) ) deallocate(cfg%lps, cfg%sigmas)
        allocate(cfg%lps(9), cfg%sigmas(9)); cfg%lps = [12.,12.,10.,10.,8.,8.,8.,8.,8.]; cfg%sigmas = [1.2,1.2,1.0,1.0,0.8,0.8,0.7,0.7,0.7]
        open(newunit=u, file=trim(fname), status='old', action='read')
        do
            read(u,'(A)',iostat=io) line
            if( io /= 0 ) exit
            if( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
            call cfg_tokenize(line, tok, ntok)
            if( ntok < 2 ) cycle
            select case(trim(tok(1)))
                case('n1');            read(tok(2),*) cfg%n1
                case('n2');            read(tok(2),*) cfg%n2
                case('nfree');         read(tok(2),*) cfg%nfree
                case('shrink');        read(tok(2),*) cfg%shrink
                case('eps0');          read(tok(2),*) cfg%eps0f
                case('balance');       read(tok(2),*) cfg%balance
                case('balance_until'); read(tok(2),*) cfg%balance_until
                case('focus_mask');    cfg%focus = tok(2)
                case('solvent_mask');  cfg%solvent = tok(2)
                case('subset');        cfg%subset = tok(2)
                case('nsub');          read(tok(2),*) cfg%nsub
                case('seed');          read(tok(2),*) cfg%seed
                case('ss');            read(tok(2),*) cfg%sshift
                case('fs');            read(tok(2),*) cfg%fscale
                case('conv');          read(tok(2),*) cfg%conv
                case('sign');          read(tok(2),*) cfg%csign
                case('ctfabs');        read(tok(2),*) cfg%ctfabs
                case('write_halves');  read(tok(2),*) cfg%write_halves
                case('pose_every');    read(tok(2),*) cfg%pose_every
                case('pose_step');     read(tok(2),*) cfg%pose_step
                case('pose_prior');    read(tok(2),*) cfg%pose_prior
                case('devw');          read(tok(2),*) cfg%devw
                case('pnoise');        read(tok(2),*) cfg%pnoise
                case('hp');            read(tok(2),*) cfg%hp
                case('pscale');        read(tok(2),*) cfg%pscale
                case('stop_frac');     read(tok(2),*) cfg%stop_frac
                case('kmax_load');     read(tok(2),*) cfg%kmax_load
                case('gridding');      read(tok(2),*) cfg%gridding
                case('server_list');   cfg%server_list = trim(tok(2))
                case('nthr_part');     read(tok(2),*) cfg%nthr_part
                case('init_file');     cfg%init = tok(2)
                case('lps')
                    deallocate(cfg%lps); allocate(cfg%lps(ntok-1))
                    do j = 2,ntok
                        read(tok(j),*) cfg%lps(j-1)
                    enddo
                case('sigmas')
                    deallocate(cfg%sigmas); allocate(cfg%sigmas(ntok-1))
                    do j = 2,ntok
                        read(tok(j),*) cfg%sigmas(j-1)
                    enddo
            end select
        enddo
        close(u)
        if( size(cfg%lps) /= size(cfg%sigmas) ) THROW_HARD('lps and sigmas must have the same length')
    end subroutine klem_read_cfg

    !> lattice start positions "pind g1 g2" for the given pinds (random where absent)
    subroutine klem_read_init_positions( initfile, n1_, n2_, pinds, nall, g1o, g2o )
        character(len=*), intent(in)  :: initfile
        integer,          intent(in)  :: n1_, n2_, pinds(:), nall
        real,             intent(out) :: g1o(:), g2o(:)
        real, allocatable :: a1(:), a2(:)
        integer :: u_, io_, i_, c_, n
        real    :: t1, t2
        n = size(pinds); call random_number(g1o); call random_number(g2o); g1o = g1o * real(n1_ - 1); g2o = g2o * real(n2_ - 1)
        if( trim(initfile) == 'none' ) return
        allocate(a1(nall), a2(nall), source=-1.)
        open(newunit=u_, file=trim(initfile), status='old', action='read')
        do
            read(u_,*,iostat=io_) i_, t1, t2
            if( io_ /= 0 ) exit
            if( i_ >= 1 .and. i_ <= nall )then
                a1(i_) = t1; a2(i_) = t2
            endif
        enddo
        close(u_); c_ = 0
        do i_ = 1,n
            if( a1(pinds(i_)) >= 0. )then
                g1o(i_) = min(max(a1(pinds(i_)), 0.), real(n1_ - 1)); g2o(i_) = min(max(a2(pinds(i_)), 0.), real(n2_ - 1)); c_ = c_ + 1
            endif
        enddo
        write(logfhandle,'(A,I0,A,I0,A)') '>>> init (served board): positions from '//trim(initfile)//' for ', c_, ' of ', n, ' particles'
        deallocate(a1, a2)
    end subroutine klem_read_init_positions

    subroutine cfg_tokenize( str, toks, nt )
        character(len=*), intent(in)  :: str
        character(len=*), intent(out) :: toks(:)
        integer,          intent(out) :: nt
        integer :: i, l, s
        logical :: intok
        nt = 0; intok = .false.; l = len_trim(str); s = 1
        do i = 1,l+1
            if( i > l )then
                if( intok )then
                    nt = nt + 1; toks(nt) = str(s:l)
                endif
            else if( str(i:i) == ' ' .or. str(i:i) == char(9) )then
                if( intok )then
                    nt = nt + 1; toks(nt) = str(s:i-1); intok = .false.
                endif
            else if( .not. intok )then
                intok = .true.; s = i
            endif
        enddo
    end subroutine cfg_tokenize

    ! ------------------------------------------------------------------ half-sphere packing of Fourier cubes
    subroutine half_sphere_index( kmax, hh, kk, ll, nhs )
        integer,              intent(in)  :: kmax
        integer, allocatable, intent(out) :: hh(:), kk(:), ll(:)
        integer,              intent(out) :: nhs
        integer :: h, k, l, c
        c = 0
        do l = -kmax,kmax
            do k = -kmax,kmax
                do h = 0,kmax
                    if( h*h + k*k + l*l <= kmax*kmax ) c = c + 1
                enddo
            enddo
        enddo
        nhs = c; allocate(hh(nhs), kk(nhs), ll(nhs)); c = 0
        do l = -kmax,kmax
            do k = -kmax,kmax
                do h = 0,kmax
                    if( h*h + k*k + l*l <= kmax*kmax )then
                        c = c + 1; hh(c) = h; kk(c) = k; ll(c) = l
                    endif
                enddo
            enddo
        enddo
    end subroutine half_sphere_index

    ! ------------------------------------------------------------------ file waiting
    subroutine wait_for_file( fname, what, wmax )
        character(len=*), intent(in) :: fname, what
        integer, optional, intent(in) :: wmax
        integer :: wlim
        integer :: waited
        waited = 0
        wlim = WAIT_MAX; if( present(wmax) ) wlim = wmax
        do
            if( file_exists(trim(fname)) ) return
            call sleep(1); waited = waited + 1
            if( waited > wlim ) THROW_HARD('timed out waiting for '//trim(what)//': '//trim(fname))
        enddo
    end subroutine wait_for_file

    subroutine touch_file( fname )
        character(len=*), intent(in) :: fname
        integer :: u
        open(newunit=u, file=trim(fname), status='replace', action='write'); write(u,'(A)') 'ok'; close(u)
    end subroutine touch_file

    ! ================================================================== MASTER
    !> pinds/eos/g1/g2: the particle rows, halves and lattice start positions chosen by the sequential setup
    subroutine klem_board_master( params, spproj, cfg, pinds, eos, g1, g2, odir, box, smpd, cfgs_more, odirs_more, g1_more, g2_more )
        class(parameters), intent(inout) :: params
        class(sp_project), intent(inout) :: spproj
        type(klem_cfg),    intent(in)    :: cfg
        integer,           intent(in)    :: pinds(:), eos(:), box
        real,              intent(in)    :: g1(:), g2(:), smpd
        character(len=*),  intent(in)    :: odir
        type(qsys_env)  :: qenv
        type(klem_band) :: bd
        type(image)     :: mimg
        type(string)    :: scr, outf, fn, eb
        type(cmdline), allocatable :: clines(:)
        integer, allocatable :: ord(:), key(:), stk_of(:), ind_of(:), hh(:), kk(:), lh(:), i1(:), i2(:), pstart(:), pend(:), am(:), pind_all(:)
        real,    allocatable :: sig2(:), scl(:), sig2_base(:), wsh(:), sc_sh(:), occ(:), fm(:,:,:), z1(:), z2(:), maxr(:), gmass(:), crep(:), Ksm(:,:), rrow(:)
        real(dp), allocatable :: num_sh(:), den_sh(:), sig2_sh(:), cnt_sh(:), dpos(:), fnum(:), fd0(:), fd1(:)
        complex(dp), allocatable :: num(:,:,:,:,:), numc(:,:,:,:)
        real(dp),    allocatable :: den(:,:,:,:,:), denc(:,:,:,:)
        complex, allocatable :: V(:,:,:,:,:), Gc(:,:,:,:), Vprev(:,:,:,:,:)
        real(dp) :: dnum, dden
        real     :: vchg
        complex(dp), allocatable :: pnum(:)
        real(dp),    allocatable :: pden(:)
        complex,     allocatable :: pnum4(:)
        real,        allocatable :: pden4(:)
        character(len=512) :: ddir, fname, nb_infile, nb_ddir, odirb, line, list_file
        logical :: l_server, l_stop
        integer :: nserved, nl, nall
        type(klem_cfg), optional, intent(in) :: cfgs_more(:)
        character(len=512), optional, intent(in) :: odirs_more(:)
        real, optional, intent(in) :: g1_more(:,:), g2_more(:,:)
        type(klem_cfg) :: cfgb
        real, allocatable :: g1b(:), g2b(:)
        integer :: ib, nboards
        character(len=:), allocatable :: dirc
        real(dp) :: eps0med, acc, tsum, tsum2, t0, tm1, tm2
        real     :: lp, sigma, eps, eps0, cur_lp, fsc_s, tmp, pgain, pgain_all
        integer  :: n, nparts, ngrid, ncl, nit, it, ireq, ip, i, j, q, u, ih, ic, h, k, l, kmax, kmaxf, n3, nhs, ish, ncal, nmoved, nmoved_all, c2, nthr_w, io, cnt, nchg, nchg_all
        logical  :: l_focus, l_solv
        real,    allocatable :: dm(:,:,:)
        complex, allocatable :: G0d(:,:,:,:)
        real,    allocatable :: sm(:,:,:)
        complex, allocatable :: G0(:,:,:,:)
        cfgb = cfg; odirb = odir; g1b = g1; g2b = g2; nboards = 1; if( present(cfgs_more) ) nboards = 1 + size(cfgs_more)
        l_server = trim(cfg%server_list) /= 'none'; list_file = trim(cfg%server_list); nserved = 0; nall = spproj%os_ptcl3D%get_noris()
        n = size(pinds); nparts = params%nparts; ngrid = cfgb%n1 * cfgb%n2; ncl = ngrid + cfgb%nfree; nit = size(cfgb%lps); l_focus = trim(cfgb%focus) /= 'none'
        ddir = trim(odirb)//'/distr'; call simple_mkdir(trim(ddir)); dirc = trim(odirb)
        nthr_w = cfgb%nthr_part; if( nthr_w <= 0 ) nthr_w = params%nthr
        write(logfhandle,'(A,I0,A,I0,A,I0,A)') '>>> KLEM BOARD (distributed): ', n, ' particles over ', nparts, ' workers x ', nthr_w, ' threads'
        ! ---- the particles in stack order, split into contiguous parts
        allocate(stk_of(n), ind_of(n), ord(n), key(n))
        do i = 1,n
            call spproj%map_ptcl_ind2stk_ind('ptcl3D', pinds(i), stk_of(i), ind_of(i)); key(i) = stk_of(i) * 10000000 + ind_of(i); ord(i) = i
        enddo
        call hpsort(key, ord); deallocate(key, stk_of, ind_of)
        allocate(pstart(nparts), pend(nparts))
        do ip = 1,nparts
            pstart(ip) = (ip - 1) * n / nparts + 1; pend(ip) = ip * n / nparts
            open(newunit=u, file=trim(ddir)//'/part_'//int2str_pad(ip,2)//'.txt', status='replace', action='write')
            do j = pstart(ip),pend(ip)
                i = ord(j); write(u,'(I0,1X,I0,1X,F8.4,1X,F8.4)') pinds(i), eos(i), g1b(i), g2b(i)
            enddo
            close(u)
        enddo
        allocate(pind_all(n)); pind_all = pinds(ord)
        ! ---- launch the persistent workers
        ! the queue layer submits "<cwd>/<script>"; workers run through simple_private_exec (the public executable appends a job
        ! record to the project at exit, and concurrent worker exits corrupt it)
        call qenv%new(params, nparts, qsys_nthr=nthr_w)
        allocate(clines(nparts))
        do ip = 1,nparts
            call clines(ip)%set('prg', 'klem_board_part'); call clines(ip)%set('projfile', params%projfile); call clines(ip)%set('infile', params%infile)
            call clines(ip)%set('dir', dirc); call clines(ip)%set('part', ip); call clines(ip)%set('nparts', nparts); call clines(ip)%set('nthr', nthr_w); call clines(ip)%set('mkdir', 'no')
        enddo
        call simple_mkdir('stderrout')
        do ip = 1,nparts   ! one job per worker (a multi-command script would run them one after the other)
            scr = 'klem_worker_'//int2str_pad(ip,2); outf = 'klem_worker_'//int2str_pad(ip,2)//'.log'
            call qenv%exec_simple_prg_in_queue_async(clines(ip), scr, outf)
        enddo
        write(logfhandle,'(A,I0,A)') '>>> ', nparts, ' workers submitted (scripts klem_worker_NN in the working directory, logs klem_worker_NN.log)'; flush(logfhandle)
        ! ---- band object at the finest band (coefficient bookkeeping only)
        kmaxf = floor(real(box) * smpd / minval(cfgb%lps))
        if( present(cfgs_more) )then   ! load once at the finest band of ALL boards the workers will serve
            do ib = 1,size(cfgs_more)
                kmaxf = max(kmaxf, floor(real(box) * smpd / minval(cfgs_more(ib)%lps)))
            enddo
        endif
        if( l_server .and. cfg%kmax_load > 0 ) kmaxf = max(kmaxf, cfg%kmax_load)
        if( present(cfgs_more) .or. l_server )then
            open(newunit=u, file=trim(ddir)//'/kmax_load.txt', status='replace', action='write'); write(u,'(I0)') kmaxf; close(u)
        endif
        call klem_band_new(bd, kmaxf)
        allocate(Ksm(ngrid,ngrid), i1(ngrid), i2(ngrid)); call klem_lattice_kernel(cfgb%n1, cfgb%n2, 1., Ksm, i1, i2)
        allocate(sig2(bd%nband), scl(bd%nband), sig2_base(bd%nband), occ(ncl)); sig2 = 1.; scl = 1.
        ireq = 0; ib = 1
        boards: do
        write(logfhandle,'(A,I0,A,I0,A)') '>>> board ', ib, ' of ', nboards, ': '//trim(odirb); flush(logfhandle)
        nit = size(cfgb%lps); l_focus = trim(cfgb%focus) /= 'none'; l_solv = trim(cfgb%solvent) /= 'none'; cur_lp = -1.; it = 0
        call klem_set_gridding(cfgb%gridding)
        if( klem_gridding_on() )then
            allocate(dm(box,box,box)); call klem_deapod_envelope(box, 4., dm); write(logfhandle,'(A)') '>>> gridding on (1 = Kaiser-Bessel 3x3x3, 2 = trilinear 2x2x2): stencil insertion/interpolation + real-space deapodisation (gain capped at 4)'
        endif
        if( allocated(Vprev) ) deallocate(Vprev)
        do
            it = it + 1; if( it > nit ) exit
            lp = cfgb%lps(it); sigma = cfgb%sigmas(it); t0 = wtime()
            if( abs(lp - cur_lp) > 1.e-6 )then
                kmax = floor(real(box) * smpd / lp); n3 = 2 * kmax + 1
                if( allocated(hh) ) deallocate(hh, kk, lh)
                call half_sphere_index(kmax, hh, kk, lh, nhs)
                write(logfhandle,'(A,F5.1,A,I0,A,I0,A,I0,A)') '>>> band ', lp, ' A: kmax ', kmax, ', grid ', n3, '^3, ', nhs, ' half-sphere coefficients'; flush(logfhandle)
                ! consensus of each half from all workers
                if( allocated(numc) ) deallocate(numc, denc)
                if( allocated(Gc) ) deallocate(Gc, V, num, den)
                if( allocated(pnum) ) deallocate(pnum, pden, pnum4, pden4)
                allocate(numc(-kmax:kmax,-kmax:kmax,-kmax:kmax,2), source=cmplx(0.d0,0.d0,dp)); allocate(denc(-kmax:kmax,-kmax:kmax,-kmax:kmax,2), source=0.d0)
                allocate(pnum(nhs), pden(nhs), pnum4(nhs), pden4(nhs))
                ireq = ireq + 1; call write_req(REQ_CONS, it, kmax)
                do ip = 1,nparts
                    call wait_rep(ip, ireq); open(newunit=u, file=rep_name(ip, ireq), access='stream', form='unformatted', status='old', action='read')
                    do ih = 1,2
                        read(u) pnum4; pnum = pnum4; read(u) pden4; pden = pden4; call unpack_add(pnum, pden, numc(:,:,:,ih), denc(:,:,:,ih))
                    enddo
                    close(u)
                enddo
                call rm_reps(ireq)
                allocate(dpos(2 * n3**3)); cnt = 0
                do ih = 1,2
                    do l = -kmax,kmax
                        do k = -kmax,kmax
                            do h = -kmax,kmax
                                if( denc(h,k,l,ih) > 0.d0 )then
                                    cnt = cnt + 1; dpos(cnt) = denc(h,k,l,ih)
                                endif
                            enddo
                        enddo
                    enddo
                enddo
                eps0med = median_dp(dpos(1:cnt)); deallocate(dpos); eps0 = real(cfgb%eps0f * eps0med); eps = real(cfgb%shrink * eps0med)
                allocate(Gc(-kmax:kmax,-kmax:kmax,-kmax:kmax,2))
                do ih = 1,2
                    Gc(:,:,:,ih) = cmplx(numc(:,:,:,ih) / (denc(:,:,:,ih) + real(eps0,dp)))
                enddo
                deallocate(numc, denc)
                allocate(num(-kmax:kmax,-kmax:kmax,-kmax:kmax,ncl,2), source=cmplx(0.d0,0.d0,dp)); allocate(den(-kmax:kmax,-kmax:kmax,-kmax:kmax,ncl,2), source=0.d0)
                allocate(V(-kmax:kmax,-kmax:kmax,-kmax:kmax,ncl,2), source=cmplx(0.,0.))
                ! calibration: reference scale per shell, then the noise per shell (workers sum over a subsample)
                if( allocated(num_sh) ) deallocate(num_sh, den_sh, sig2_sh, cnt_sh, sc_sh, wsh)
                allocate(num_sh(0:kmax), den_sh(0:kmax), sig2_sh(0:kmax), cnt_sh(0:kmax), sc_sh(0:kmax), wsh(0:kmax)); num_sh = 0.d0; den_sh = 0.d0
                ireq = ireq + 1; call write_req(REQ_SCALE, it, kmax)
                do ip = 1,nparts
                    call wait_rep(ip, ireq); open(newunit=u, file=rep_name(ip, ireq), access='stream', form='unformatted', status='old', action='read')
                    read(u) pnum4(1:kmax+1); pnum(1:kmax+1) = pnum4(1:kmax+1); read(u) pden4(1:kmax+1); pden(1:kmax+1) = pden4(1:kmax+1); close(u)
                    num_sh = num_sh + real(pnum(1:kmax+1), dp); den_sh = den_sh + pden(1:kmax+1)
                enddo
                call rm_reps(ireq)
                do ish = 0,kmax
                    sc_sh(ish) = merge(real(num_sh(ish) / max(den_sh(ish), 1.d-30)), 0., den_sh(ish) > 0.d0)
                enddo
                num_sh = sc_sh
                do ish = 0,kmax
                    sc_sh(ish) = real(num_sh(max(ish-1,0)) + num_sh(ish) + num_sh(min(ish+1,kmax))) / 3.
                enddo
                tmp = 0.05 * max(sum(sc_sh(1:min(4,kmax))) / real(min(4,kmax)), 1.e-30); sc_sh = max(sc_sh, tmp)
                do q = 1,bd%nband
                    scl(q) = 1.
                    if( bd%shell(q) <= kmax ) scl(q) = sc_sh(bd%shell(q))
                enddo
                sig2_sh = 0.d0; cnt_sh = 0.d0
                ireq = ireq + 1; call write_req(REQ_NOISE, it, kmax)
                do ip = 1,nparts
                    call wait_rep(ip, ireq); open(newunit=u, file=rep_name(ip, ireq), access='stream', form='unformatted', status='old', action='read')
                    read(u) pden4(1:kmax+1); pden(1:kmax+1) = pden4(1:kmax+1); sig2_sh = sig2_sh + pden(1:kmax+1); read(u) pden4(1:kmax+1); pden(1:kmax+1) = pden4(1:kmax+1); cnt_sh = cnt_sh + pden(1:kmax+1); close(u)
                enddo
                call rm_reps(ireq)
                do ish = 0,kmax
                    sig2_sh(ish) = max(merge(sig2_sh(ish) / cnt_sh(ish), 1.d0, cnt_sh(ish) > 0.d0), 1.d-12)
                enddo
                do q = 1,bd%nband
                    sig2(q) = 1.
                    if( bd%shell(q) <= kmax ) sig2(q) = real(sig2_sh(bd%shell(q)))
                    if( cfg%hp > 0. .and. real(bd%shell(q)) < real(box) * smpd / cfg%hp ) sig2(q) = 1.e30
                enddo
                sig2_base = sig2; wsh = 1.
                write(logfhandle,'(A,F5.1,A,ES9.2)') '>>> band ', lp, ' A: shrinkage eps ', eps
                write(logfhandle,'(A,*(F6.2,1X))') '>>> reference scale relative to shell 1 (every few shells): ', (sc_sh(min(ish,kmax)) / max(sc_sh(1),1e-30), ish=1,kmax,max(1,kmax/6))
                if( l_focus )then
                    if( allocated(fm) ) deallocate(fm)
                    fn = trim(cfgb%focus); call mimg%new([box,box,box], smpd); call mimg%read(fn); allocate(fm(box,box,box)); call mimg%get_rmat_sub(fm); call mimg%kill
                    fm = min(max(fm, 0.), 1.)
                endif
                if( l_solv )then   ! solvent mask on the full cell maps (applied after the focus; the consensus is not masked)
                    if( allocated(sm) ) deallocate(sm)
                    fn = trim(cfgb%solvent); call mimg%new([box,box,box], smpd); call mimg%read(fn); allocate(sm(box,box,box)); call mimg%get_rmat_sub(sm); call mimg%kill
                    sm = min(max(sm, 0.), 1.)
                    write(logfhandle,'(A,F6.2,A)') '>>> solvent mask: ', 100. * real(count(sm > 0.5)) / real(box)**3, ' % of the box above 0.5'
                    write(logfhandle,'(A,F6.2,A)') '>>> focal mask: ', 100. * real(count(fm > 0.5)) / real(box)**3, ' % of the box above 0.5 (references differ only inside it)'
                endif
                cur_lp = lp
            endif
            ! ---- M-step: workers accumulate with their kernel-smoothed responsibilities at this round's sigma
            num = cmplx(0.d0,0.d0,dp); den = 0.d0
            ireq = ireq + 1; call write_req(REQ_MSTEP, it, kmax)
            do ip = 1,nparts
                call wait_rep(ip, ireq); open(newunit=u, file=rep_name(ip, ireq), access='stream', form='unformatted', status='old', action='read')
                do ih = 1,2
                    do ic = 1,ncl
                        read(u) pnum4; pnum = pnum4; read(u) pden4; pden = pden4; call unpack_add(pnum, pden, num(:,:,:,ic,ih), den(:,:,:,ic,ih))
                    enddo
                enddo
                close(u)
            enddo
            call rm_reps(ireq); tm1 = wtime()
            do ih = 1,2
                do ic = 1,ncl
                    V(:,:,:,ic,ih) = cmplx((num(:,:,:,ic,ih) + real(eps,dp) * Gc(:,:,:,ih)) / (den(:,:,:,ic,ih) + real(eps,dp)))
                enddo
            enddo
            if( klem_gridding_on() )then   ! deapodise the gridded cells before the focus
                allocate(G0d(-kmax:kmax,-kmax:kmax,-kmax:kmax,2), source=cmplx(0.,0.)); call apply_focus_d(V, G0d, dm, kmax, ncl, box, smpd); deallocate(G0d)
            endif
            if( l_focus ) call apply_focus_d(V, Gc, fm, kmax, ncl, box, smpd)
            if( l_solv )then
                allocate(G0(-kmax:kmax,-kmax:kmax,-kmax:kmax,2), source=cmplx(0.,0.)); call apply_focus_d(V, G0, sm, kmax, ncl, box, smpd); deallocate(G0)
            endif
            tm2 = wtime()
            if( allocated(Vprev) )then   ! model change: relative change of the grid-cell deviations since the last iteration (1 at a band change)
                vchg = 1.
                if( all(shape(Vprev) == shape(V)) )then
                    dnum = 0.d0; dden = 0.d0
                    do ih = 1,2
                        do ic = 1,ngrid
                            dnum = dnum + sum(real(abs(V(:,:,:,ic,ih) - Vprev(:,:,:,ic,ih)),dp)**2); dden = dden + sum(real(abs(V(:,:,:,ic,ih) - Gc(:,:,:,ih)),dp)**2)
                        enddo
                    enddo
                    vchg = real(sqrt(dnum / max(dden, 1.d-30)))
                endif
                deallocate(Vprev)
            else
                vchg = 1.
            endif
            allocate(Vprev, source=V)
            if( cfgb%devw > 0. .and. it > 1 )then
                allocate(fnum(0:kmax), fd0(0:kmax), fd1(0:kmax)); fnum = 0.d0; fd0 = 0.d0; fd1 = 0.d0
                do ic = 1,ngrid
                    do l = -kmax,kmax
                        do k = -kmax,kmax
                            do h = -kmax,kmax
                                ish = nint(sqrt(real(h*h + k*k + l*l)))
                                if( ish > kmax ) cycle
                                fnum(ish) = fnum(ish) + real(conjg(V(h,k,l,ic,1) - Gc(h,k,l,1)) * (V(h,k,l,ic,2) - Gc(h,k,l,2)), dp)
                                fd0(ish) = fd0(ish) + real(abs(V(h,k,l,ic,1) - Gc(h,k,l,1)),dp)**2; fd1(ish) = fd1(ish) + real(abs(V(h,k,l,ic,2) - Gc(h,k,l,2)),dp)**2
                            enddo
                        enddo
                    enddo
                enddo
                do ish = 0,kmax
                    fsc_s = real(fnum(ish) / max(sqrt(fd0(ish) * fd1(ish)), 1.d-30)); wsh(ish) = min(max(2. * fsc_s / (1. + max(fsc_s, 1.e-6)), cfgb%devw), 1.)
                enddo
                wsh(0:min(1,kmax)) = 1.
                do q = 1,bd%nband
                    if( bd%shell(q) <= kmax ) sig2(q) = sig2_base(q) / wsh(bd%shell(q))
                enddo
                write(logfhandle,'(A,*(F5.2,1X))') '>>> deviation-FSC likelihood weights (every few shells): ', (wsh(min(ish,kmax)), ish=1,kmax,max(1,kmax/6))
                deallocate(fnum, fd0, fd1)
            endif
            ! ---- E-step: workers update their responsibilities against V (cross-half), balance per part, pose step
            ireq = ireq + 1; call write_req(REQ_ESTEP, it, kmax)
            occ = 0.; nmoved_all = 0; pgain_all = 0.; nchg_all = 0
            do ip = 1,nparts
                call wait_rep(ip, ireq); open(newunit=u, file=rep_name(ip, ireq), access='stream', form='unformatted', status='old', action='read')
                read(u) pden4(1:ncl); pden(1:ncl) = pden4(1:ncl); read(u) nmoved; read(u) pgain; read(u) nchg; close(u)
                occ = occ + real(pden(1:ncl)); nmoved_all = nmoved_all + nmoved; nchg_all = nchg_all + nchg; pgain_all = pgain_all + pgain * real(pend(ip) - pstart(ip) + 1)
            enddo
            call rm_reps(ireq); occ = occ / real(n)
            write(logfhandle,'(A,I0,A,F7.1,A,F7.1,A,F7.1,A)') '>>> it ', it, ' timing: M-step round ', tm1 - t0, ' s, solve + focus ', tm2 - tm1, ' s, E-step round ', wtime() - tm2, ' s'
            if( cfgb%pose_every > 0 )then
                if( mod(it, cfgb%pose_every) == 0 ) write(logfhandle,'(A,I0,A,F5.1,A,F6.3,A,F8.2)') '>>> pose step (it ', it, ', step ', cfgb%pose_step, ' deg): fraction moved ', real(nmoved_all) / real(n), ', mean log-likelihood gain ', pgain_all / real(n)
            endif
            write(logfhandle,'(A,I0,A,F5.1,A,F4.1,A,3F7.3,A,F6.3,A,F6.3,A,F6.3,A,F7.1,A)') '>>> it ', it, ' (lp ', lp, ', sigma ', sigma, '): grid occupancy min/median/max ', minval(occ(1:ngrid)), median_r(occ(1:ngrid)), maxval(occ(1:ngrid)), &
                &', free mass ', sum(occ(ngrid+1:ncl)), ', argmax changed ', real(nchg_all) / real(n), ', model change ', vchg, ' (', wtime() - t0, ' s)'
            flush(logfhandle)
            if( vchg < cfgb%stop_frac )then   ! converged at this band: jump to its last iteration (final sigma), skipping the middle ones
                j = it
                do while( j < nit )
                    if( abs(cfgb%lps(j+1) - lp) > 1.e-6 ) exit
                    j = j + 1
                enddo
                if( j - 1 > it )then
                    write(logfhandle,'(A,I0,A,I0)') '>>> early stop at this band: skipping iterations ', it + 1, ' to ', j - 1; it = j - 1
                endif
            endif
        enddo
        ! ---- finish: workers write their latent/resp/poses parts; the master merges and writes the board outputs
        ireq = ireq + 1; call write_req(REQ_FLUSH, nit, kmax)
        do ip = 1,nparts
            call wait_rep(ip, ireq)
        enddo
        allocate(am(n), z1(n), z2(n), maxr(n), gmass(n), rrow(ncl)); j = 0
        open(newunit=u, file=trim(odirb)//'/klem_latent.txt', status='replace', action='write')
        write(u,'(A)') '# pind z1 z2 argmax maxr gridmass layer2mass'
        do ip = 1,nparts
            open(newunit=io, file=trim(ddir)//'/part_'//int2str_pad(ip,2)//'_latent.txt', status='old', action='read')
            do
                read(io,*,iostat=cnt) i, tmp, fsc_s, ic, pgain, sigma
                if( cnt /= 0 ) exit
                j = j + 1; z1(j) = tmp; z2(j) = fsc_s; am(j) = ic; maxr(j) = pgain; gmass(j) = sigma
                write(u,'(I0,1X,F8.4,1X,F8.4,1X,I0,1X,F7.4,1X,F7.4,1X,F7.4)') i, tmp, fsc_s, ic, pgain, sigma, 0.
            enddo
            close(io)
        enddo
        close(u)
        if( j /= n ) write(logfhandle,'(A,I0,A,I0)') 'WARNING: merged latent rows ', j, ' of ', n
        open(newunit=u, file=trim(odirb)//'/klem_resp.bin', access='stream', form='unformatted', status='replace', action='write')
        write(u) n, ncl, ngrid, cfgb%n1, cfgb%n2
        write(u) pind_all
        write(u) i1, i2
        do ip = 1,nparts
            open(newunit=io, file=trim(ddir)//'/part_'//int2str_pad(ip,2)//'_resp.bin', access='stream', form='unformatted', status='old', action='read')
            do j = pstart(ip),pend(ip)
                read(io) rrow; write(u) rrow
            enddo
            close(io)
        enddo
        close(u)
        allocate(crep(ncl)); c2 = n3 / 2
        do ic = 1,ncl
            acc = 0.d0; tsum = 0.d0; tsum2 = 0.d0
            do l = -kmax,kmax
                do k = -kmax,kmax
                    do h = -kmax,kmax
                        if( h*h + k*k + l*l < 1 .or. h*h + k*k + l*l > c2*c2 ) cycle
                        acc = acc + real(conjg(V(h,k,l,ic,1) - Gc(h,k,l,1)) * (V(h,k,l,ic,2) - Gc(h,k,l,2)), dp)
                        tsum = tsum + real(abs(V(h,k,l,ic,1) - Gc(h,k,l,1)),dp)**2; tsum2 = tsum2 + real(abs(V(h,k,l,ic,2) - Gc(h,k,l,2)),dp)**2
                    enddo
                enddo
            enddo
            crep(ic) = real(acc / max(sqrt(tsum * tsum2), 1.d-30))
        enddo
        write(logfhandle,'(A,F5.2,A,I0,A,I0,A)') '>>> states: cross-half reproducibility of the deviation median ', median_r(crep(1:ngrid)), ', ', count(crep(1:ngrid) > 0.5), ' of ', ngrid, ' above 0.5'
        open(newunit=u, file=trim(odirb)//'/klem_states.txt', status='replace', action='write')
        write(u,'(A)') '# cell population reproducibility isgrid'
        do ic = 1,ncl
            write(u,'(I0,1X,I0,1X,F8.4,1X,I0)') ic, count(am == ic), crep(ic), merge(1, 0, ic <= ngrid)
        enddo
        close(u)
        open(newunit=u, file=trim(odirb)//'/klem_halves.bin', access='stream', form='unformatted', status='replace', action='write')
        write(u) ncl, ngrid, kmax, n3
        write(u) (real(count(am == ic)), ic=1,ngrid)
        do ih = 1,2
            do ic = 1,ngrid
                do l = -kmax,kmax
                    do k = -kmax,kmax
                        do h = -kmax,kmax
                            write(u) V(h,k,l,ic,ih) - Gc(h,k,l,ih)
                        enddo
                    enddo
                enddo
            enddo
        enddo
        close(u)
        call simple_mkdir(trim(odirb)//'/cells'); if( cfgb%write_halves /= 0 ) call simple_mkdir(trim(odirb)//'/cells_half')
        do ic = 1,ncl
            if( ic <= ngrid )then
                fname = trim(odirb)//'/cells/cell_'//int2str_pad(i1(ic),2)//'_'//int2str_pad(i2(ic),2)//'.mrc'
            else
                fname = trim(odirb)//'/cells/free_'//int2str_pad(ic - ngrid, 1)//'.mrc'
            endif
            call klem_write_volume(0.5 * (V(:,:,:,ic,1) + V(:,:,:,ic,2)), kmax, box, smpd, trim(fname))
            if( cfgb%write_halves /= 0 )then
                call klem_write_volume(V(:,:,:,ic,1), kmax, box, smpd, trim(odirb)//'/cells_half/cell_'//int2str_pad(ic,2)//'_even.mrc')
                call klem_write_volume(V(:,:,:,ic,2), kmax, box, smpd, trim(odirb)//'/cells_half/cell_'//int2str_pad(ic,2)//'_odd.mrc')
            endif
        enddo
        if( cfgb%pose_every > 0 )then
            open(newunit=u, file=trim(odirb)//'/klem_poses.txt', status='replace', action='write')
            write(u,'(A)') '# pind e1 e2 e3'
            do ip = 1,nparts
                open(newunit=io, file=trim(ddir)//'/part_'//int2str_pad(ip,2)//'_poses.txt', status='old', action='read')
                do
                    read(io,'(A)',iostat=cnt) fname
                    if( cnt /= 0 ) exit
                    if( fname(1:1) == '#' ) cycle
                    write(u,'(A)') trim(fname)
                enddo
                close(io)
            enddo
            close(u)
        endif
        call rm_reps(ireq)
        write(logfhandle,'(A)') '>>> KLEM BOARD (distributed) written to '//trim(odirb)
        open(newunit=u, file=trim(odirb)//'/BOARD_DONE', status='replace', action='write'); close(u)
        deallocate(am, z1, z2, maxr, gmass, rrow, crep)
        l_stop = .true.
        if( ib < nboards )then   ! next board on the same workers: part files (init positions) in its distr dir, then the switch request in the current one
            ib = ib + 1; cfgb = cfgs_more(ib - 1); odirb = odirs_more(ib - 1); g1b = g1_more(:,ib - 1); g2b = g2_more(:,ib - 1); l_stop = .false.
        else if( l_server )then   ! served mode: the next board dir is the next unread line of the list file; the line STOP ends the service
            do
                nl = 0; line = ''
                open(newunit=u, file=trim(list_file), status='old', action='read', iostat=io)
                if( io == 0 )then
                    do
                        read(u,'(A)',iostat=io) line
                        if( io /= 0 )then
                            line = ''; exit
                        endif
                        if( len_trim(line) == 0 ) cycle
                        nl = nl + 1
                        if( nl > nserved ) exit
                        line = ''
                    enddo
                    close(u)
                endif
                if( len_trim(line) > 0 ) exit
                call sleep(5)
            enddo
            nserved = nserved + 1
            if( trim(adjustl(line)) /= 'STOP' )then
                ib = ib + 1; odirb = trim(adjustl(line)); call klem_read_cfg(trim(odirb)//'/board.txt', cfgb)
                call klem_read_init_positions(trim(cfgb%init), cfgb%n1, cfgb%n2, pinds, nall, g1b, g2b); l_stop = .false.
                write(logfhandle,'(A)') '>>> board server: next board '//trim(odirb); flush(logfhandle)
            endif
        endif
        if( .not. l_stop )then
            if( cfgb%n1 * cfgb%n2 + cfgb%nfree /= ncl ) THROW_HARD('board server: every board must have the same lattice and free components')
            nb_ddir = trim(odirb)//'/distr'; nb_infile = trim(odirb)//'/board.txt'; call simple_mkdir(trim(nb_ddir))
            do ip = 1,nparts
                open(newunit=u, file=trim(nb_ddir)//'/part_'//int2str_pad(ip,2)//'.txt', status='replace', action='write')
                do j = pstart(ip),pend(ip)
                    i = ord(j); write(u,'(I0,1X,I0,1X,F8.4,1X,F8.4)') pinds(i), eos(i), g1b(i), g2b(i)
                enddo
                close(u)
            enddo
            ireq = ireq + 1; call write_req(REQ_NEWBOARD, 0, kmax); ddir = nb_ddir
            cycle boards
        endif
        ireq = ireq + 1; call write_req(REQ_FINISH, nit, kmax)
        exit boards
        enddo boards
        call qenv%kill

      contains

        function rep_name( ip_, ir_ ) result( s_ )
            integer, intent(in) :: ip_, ir_
            character(len=:), allocatable :: s_
            s_ = trim(ddir)//'/rep_'//int2str_pad(ip_,2)//'_'//int2str_pad(ir_,4)//'.bin'
        end function rep_name

        subroutine wait_rep( ip_, ir_ )
            integer, intent(in) :: ip_, ir_
            call wait_for_file(trim(ddir)//'/rep_'//int2str_pad(ip_,2)//'_'//int2str_pad(ir_,4)//'.done', 'worker '//int2str(ip_)//' reply', merge(WAIT_START, WAIT_MAX, ir_ == 1))
        end subroutine wait_rep

        subroutine rm_reps( ir_ )
            integer, intent(in) :: ir_
            integer :: ip_
            do ip_ = 1,nparts
                call del_file(trim(ddir)//'/rep_'//int2str_pad(ip_,2)//'_'//int2str_pad(ir_,4)//'.bin')
                call del_file(trim(ddir)//'/rep_'//int2str_pad(ip_,2)//'_'//int2str_pad(ir_,4)//'.done')
            enddo
            if( ir_ > 1 )then
                call del_file(trim(ddir)//'/req_'//int2str_pad(ir_-1,4)//'.bin'); call del_file(trim(ddir)//'/req_'//int2str_pad(ir_-1,4)//'.ready')
            endif
        end subroutine rm_reps

        !> request header + payload: CONS/SCALE/NOISE/MSTEP carry what the workers need at this stage; ESTEP carries V
        subroutine write_req( typ, it_, kmax_ )
            integer, intent(in) :: typ, it_, kmax_
            integer :: u_, ih_, ic_, q_, nb_
            open(newunit=u_, file=trim(ddir)//'/req_'//int2str_pad(ireq,4)//'.bin', access='stream', form='unformatted', status='replace', action='write')
            nb_ = bd%nband
            write(u_) typ, it_, kmax_, ncl, ngrid, nb_
            select case(typ)
                case(REQ_SCALE)
                    do ih_ = 1,2
                        call pack_c(Gc(:,:,:,ih_)); pnum4 = cmplx(pnum); write(u_) pnum4
                    enddo
                case(REQ_NOISE)
                    do ih_ = 1,2
                        call pack_c(Gc(:,:,:,ih_)); pnum4 = cmplx(pnum); write(u_) pnum4
                    enddo
                    write(u_) scl
                case(REQ_NEWBOARD)
                    write(u_) nb_infile, nb_ddir
                case(REQ_MSTEP)
                    write(u_) sigma
                case(REQ_ESTEP)
                    do ih_ = 1,2
                        do ic_ = 1,ncl
                            call pack_c(V(:,:,:,ic_,ih_)); pnum4 = cmplx(pnum); write(u_) pnum4
                        enddo
                    enddo
                    write(u_) sig2; write(u_) scl
                    write(u_) merge(1, 0, cfgb%pose_every > 0 .and. mod(it_, max(cfgb%pose_every,1)) == 0)
            end select
            close(u_)
            call touch_file(trim(ddir)//'/req_'//int2str_pad(ireq,4)//'.ready')
        end subroutine write_req

        subroutine pack_c( C_ )
            complex, intent(in) :: C_(-kmax:kmax,-kmax:kmax,-kmax:kmax)
            integer :: m_
            do m_ = 1,nhs
                pnum(m_) = cmplx(C_(hh(m_),kk(m_),lh(m_)), kind=dp)
            enddo
        end subroutine pack_c

        subroutine unpack_add( pn_, pd_, A_, D_ )
            complex(dp), intent(in)    :: pn_(nhs)
            real(dp),    intent(in)    :: pd_(nhs)
            complex(dp), intent(inout) :: A_(-kmax:kmax,-kmax:kmax,-kmax:kmax)
            real(dp),    intent(inout) :: D_(-kmax:kmax,-kmax:kmax,-kmax:kmax)
            integer :: m_
            do m_ = 1,nhs
                A_(hh(m_),kk(m_),lh(m_)) = A_(hh(m_),kk(m_),lh(m_)) + pn_(m_); D_(hh(m_),kk(m_),lh(m_)) = D_(hh(m_),kk(m_),lh(m_)) + pd_(m_)
                if( hh(m_) > 0 )then
                    A_(-hh(m_),-kk(m_),-lh(m_)) = A_(-hh(m_),-kk(m_),-lh(m_)) + conjg(pn_(m_)); D_(-hh(m_),-kk(m_),-lh(m_)) = D_(-hh(m_),-kk(m_),-lh(m_)) + pd_(m_)
                endif
            enddo
        end subroutine unpack_add

    end subroutine klem_board_master

    !> cells = consensus + mask x (cell - consensus), the mask applied in real space at the full even box
    subroutine apply_focus_d( V, Gc, fm, kmax, ncl, box, smpd )
        integer, intent(in)    :: kmax, ncl, box
        complex, intent(inout) :: V(-kmax:kmax,-kmax:kmax,-kmax:kmax,ncl,2)
        complex, intent(in)    :: Gc(-kmax:kmax,-kmax:kmax,-kmax:kmax,2)
        real,    intent(in)    :: fm(box,box,box), smpd
        type(image), allocatable :: dimgs(:)
        real, pointer :: pr(:,:,:) => null()
        integer :: ic, ih, h, k, l, phys(3), nthr, ithr, j
        nthr = 1
        !$ nthr = omp_get_max_threads()
        allocate(dimgs(nthr))
        do j = 1,nthr   ! one image per thread (plan creation is not thread-safe, execution is)
            call dimgs(j)%new([box,box,box], smpd)
        enddo
        !$omp parallel do default(shared) schedule(dynamic,1) collapse(2) private(ih,ic,h,k,l,phys,ithr,pr)
        do ih = 1,2
            do ic = 1,ncl
                ithr = 1
                !$ ithr = omp_get_thread_num() + 1
                call dimgs(ithr)%zero_and_flag_ft()
                do l = -kmax,kmax
                    do k = -kmax,kmax
                        do h = 0,kmax
                            phys = dimgs(ithr)%comp_addr_phys(h,k,l)
                            call dimgs(ithr)%set_fcomp([h,k,l], phys, (V(h,k,l,ic,ih) - Gc(h,k,l,ih)) / real(box)**3)
                        enddo
                    enddo
                enddo
                call dimgs(ithr)%ifft(); call dimgs(ithr)%get_rmat_ptr(pr); pr(1:box,1:box,1:box) = pr(1:box,1:box,1:box) * fm; call dimgs(ithr)%fft()
                do l = -kmax,kmax
                    do k = -kmax,kmax
                        do h = -kmax,kmax
                            phys = dimgs(ithr)%comp_addr_phys(h,k,l)
                            V(h,k,l,ic,ih) = Gc(h,k,l,ih) + dimgs(ithr)%get_fcomp([h,k,l], phys) * real(box)**3
                        enddo
                    enddo
                enddo
            enddo
        enddo
        !$omp end parallel do
        do j = 1,nthr
            call dimgs(j)%kill
        enddo
        deallocate(dimgs)
    end subroutine apply_focus_d

    ! ================================================================== WORKER
    subroutine klem_board_worker( params )
        class(parameters), intent(inout) :: params
        type(sp_project) :: spproj
        type(klem_cfg)   :: cfg
        type(klem_band)  :: bd
        type(klem_data)  :: dat
        integer, allocatable :: pinds(:), eos(:), am(:), hh(:), kk(:), lh(:), i1(:), i2(:)
        real,    allocatable :: g1(:), g2(:), r(:,:), ll(:,:), W(:,:), W1(:,:), Ksm(:,:), sig2(:), scl(:), pfac(:), pres(:), occ(:), psc(:), pscn(:)
        complex, allocatable :: V(:,:,:,:,:), Gc(:,:,:,:), P(:)
        complex(dp), allocatable :: num(:,:,:,:,:), num1(:,:,:,:,:), pnum(:)
        real(dp),    allocatable :: den(:,:,:,:,:), den1(:,:,:,:,:), pden(:), num_sh(:), den_sh(:), sig2_sh(:), cnt_sh(:)
        complex,     allocatable :: pnum4(:)
        real,        allocatable :: pden4(:)
        integer,     allocatable :: am_prev(:)
        character(len=512) :: ddir, nb_infile, nb_ddir
        real(dp) :: acc, accn, accd
        real     :: sigma, mx, s, tmp, pmed, pgain, box_r
        integer  :: ip, nparts, n, i, m, j, q, u, io, ngrid, ncl, ireq, typ, it, kmax, kmaxf, nhs, ncl_r, ngrid_r, nb_r, ih, ic, ish, nmoved, box, ncal, dopose, cnt, nchg
        real     :: smpd
        logical  :: ctfabs
        ip = params%part; nparts = params%nparts; ddir = trim(params%dir%to_char())//'/distr'; box = params%box; smpd = params%smpd
        call klem_read_cfg(params%infile%to_char(), cfg); ctfabs = cfg%ctfabs /= 0; ngrid = cfg%n1 * cfg%n2; ncl = ngrid + cfg%nfree
        call klem_set_gridding(cfg%gridding)
        call spproj%read(params%projfile)
        ! ---- the part's particles: rows, halves, lattice start positions
        n = 0
        open(newunit=u, file=trim(ddir)//'/part_'//int2str_pad(ip,2)//'.txt', status='old', action='read')
        do
            read(u,*,iostat=io) i, j, tmp, s
            if( io /= 0 ) exit
            n = n + 1
        enddo
        rewind(u); allocate(pinds(n), eos(n), g1(n), g2(n))
        do i = 1,n
            read(u,*) pinds(i), eos(i), g1(i), g2(i)
        enddo
        close(u)
        write(logfhandle,'(A,I0,A,I0,A,I0,A)') '>>> KLEM WORKER ', ip, ' of ', nparts, ': ', n, ' particles'; flush(logfhandle)
        allocate(Ksm(ngrid,ngrid), i1(ngrid), i2(ngrid)); call klem_lattice_kernel(cfg%n1, cfg%n2, 1., Ksm, i1, i2)
        allocate(r(ncl,n), source=0.); call init_resp()
        ! ---- load once at the finest band of every board this worker will serve
        kmaxf = floor(real(box) * smpd / minval(cfg%lps))
        if( file_exists(trim(ddir)//'/kmax_load.txt') )then
            open(newunit=u, file=trim(ddir)//'/kmax_load.txt', status='old', action='read'); read(u,*) m; close(u); kmaxf = max(kmaxf, m)
        endif
        call klem_band_new(bd, kmaxf)
        call klem_load_particles(spproj, bd, pinds, eos, box, smpd, cfg%sshift, cfg%fscale, cfg%csign, ctfabs, dat)
        call spproj%kill
        write(logfhandle,'(A,I0,A,I0,A)') '>>> worker ', ip, ': loaded at kmax ', kmaxf, ', waiting for requests'; flush(logfhandle)
        allocate(ll(ncl,n), W(ncl,n), am(n), occ(ncl), sig2(bd%nband), scl(bd%nband), pfac(n), pres(n), P(bd%nband), psc(n), pscn(n), am_prev(n)); pfac = 1.; psc = 1.; am_prev = 0
        ireq = 0; kmax = -1
        do
            ireq = ireq + 1
            call wait_for_file(trim(ddir)//'/req_'//int2str_pad(ireq,4)//'.ready', 'master request')
            open(newunit=u, file=trim(ddir)//'/req_'//int2str_pad(ireq,4)//'.bin', access='stream', form='unformatted', status='old', action='read')
            read(u) typ, it, kmax, ncl_r, ngrid_r, nb_r
            if( ncl_r /= ncl .or. nb_r /= bd%nband ) THROW_HARD('worker/master mismatch in ncl or nband')
            if( .not. allocated(Gc) )then
                allocate(Gc(-kmax:kmax,-kmax:kmax,-kmax:kmax,2), source=cmplx(0.,0.)); allocate(V(-kmax:kmax,-kmax:kmax,-kmax:kmax,ncl,2), source=cmplx(0.,0.))
                allocate(num(-kmax:kmax,-kmax:kmax,-kmax:kmax,ncl,2), source=cmplx(0.d0,0.d0,dp)); allocate(den(-kmax:kmax,-kmax:kmax,-kmax:kmax,ncl,2), source=0.d0)
                call half_sphere_index(kmax, hh, kk, lh, nhs); allocate(pnum(nhs), pden(nhs), pnum4(nhs), pden4(nhs))
            else if( ubound(Gc,1) /= kmax )then
                deallocate(Gc, V, num, den, hh, kk, lh, pnum, pden, pnum4, pden4)
                allocate(Gc(-kmax:kmax,-kmax:kmax,-kmax:kmax,2), source=cmplx(0.,0.)); allocate(V(-kmax:kmax,-kmax:kmax,-kmax:kmax,ncl,2), source=cmplx(0.,0.))
                allocate(num(-kmax:kmax,-kmax:kmax,-kmax:kmax,ncl,2), source=cmplx(0.d0,0.d0,dp)); allocate(den(-kmax:kmax,-kmax:kmax,-kmax:kmax,ncl,2), source=0.d0)
                call half_sphere_index(kmax, hh, kk, lh, nhs); allocate(pnum(nhs), pden(nhs), pnum4(nhs), pden4(nhs))
            endif
            select case(typ)
                case(REQ_CONS)
                    close(u)
                    allocate(W1(1,n), source=1.); allocate(num1(-kmax:kmax,-kmax:kmax,-kmax:kmax,1,2), source=cmplx(0.d0,0.d0,dp)); allocate(den1(-kmax:kmax,-kmax:kmax,-kmax:kmax,1,2), source=0.d0)
                    call klem_mstep_accum(bd, kmax, dat, W1, 1, cfg%conv, num1, den1)
                    call open_rep(u)
                    do ih = 1,2
                        call pack_cd(num1(:,:,:,1,ih), den1(:,:,:,1,ih)); pnum4 = cmplx(pnum); write(u) pnum4; pden4 = real(pden); write(u) pden4
                    enddo
                    close(u); call done_rep(); deallocate(W1, num1, den1)
                case(REQ_SCALE)
                    do ih = 1,2
                        read(u) pnum4; pnum = pnum4; call unpack_c(pnum, Gc(:,:,:,ih))
                    enddo
                    close(u)
                    if( allocated(num_sh) ) deallocate(num_sh, den_sh, sig2_sh, cnt_sh)
                    allocate(num_sh(0:kmax), den_sh(0:kmax), sig2_sh(0:kmax), cnt_sh(0:kmax)); num_sh = 0.d0; den_sh = 0.d0
                    ncal = max(1, n * nparts / 2000)
                    do i = 1,n,ncal
                        ih = merge(2, 1, dat%eo(i) == 0); call klem_project_cell(bd, kmax, Gc(:,:,:,ih), dat%R(:,:,i), cfg%conv, P)
                        do q = 1,bd%nband
                            if( bd%shell(q) > kmax ) cycle
                            tmp = dat%C(q,i); ish = bd%shell(q)
                            num_sh(ish) = num_sh(ish) + real(bd%wq(q) * real(conjg(tmp * P(q)) * dat%F(q,i)), dp); den_sh(ish) = den_sh(ish) + real(bd%wq(q) * abs(tmp * P(q))**2, dp)
                        enddo
                    enddo
                    call open_rep(u); pnum(1:kmax+1) = cmplx(num_sh, 0.d0, dp); pden(1:kmax+1) = den_sh; pnum4(1:kmax+1) = cmplx(pnum(1:kmax+1)); write(u) pnum4(1:kmax+1); pden4(1:kmax+1) = real(pden(1:kmax+1)); write(u) pden4(1:kmax+1); close(u); call done_rep()
                case(REQ_NOISE)
                    do ih = 1,2
                        read(u) pnum4; pnum = pnum4; call unpack_c(pnum, Gc(:,:,:,ih))
                    enddo
                    read(u) scl; close(u)
                    sig2_sh = 0.d0; cnt_sh = 0.d0; ncal = max(1, n * nparts / 2000)
                    do i = 1,n,ncal
                        ih = merge(2, 1, dat%eo(i) == 0); call klem_project_cell(bd, kmax, Gc(:,:,:,ih), dat%R(:,:,i), cfg%conv, P)
                        do q = 1,bd%nband
                            if( bd%shell(q) > kmax ) cycle
                            ish = bd%shell(q); sig2_sh(ish) = sig2_sh(ish) + real(abs(dat%F(q,i) - dat%C(q,i) * scl(q) * P(q))**2, dp); cnt_sh(ish) = cnt_sh(ish) + 1.d0
                        enddo
                    enddo
                    call open_rep(u); pden(1:kmax+1) = sig2_sh; pden4(1:kmax+1) = real(pden(1:kmax+1)); write(u) pden4(1:kmax+1); pden(1:kmax+1) = cnt_sh; pden4(1:kmax+1) = real(pden(1:kmax+1)); write(u) pden4(1:kmax+1); close(u); call done_rep()
                    if( cfg%pnoise /= 0 .or. cfg%pscale /= 0 )then
                        do ish = 0,kmax
                            sig2_sh(ish) = max(merge(sig2_sh(ish) / cnt_sh(ish), 1.d0, cnt_sh(ish) > 0.d0), 1.d-12)
                        enddo
                        do q = 1,bd%nband
                            sig2(q) = 1.
                            if( bd%shell(q) <= kmax ) sig2(q) = real(sig2_sh(bd%shell(q)))
                            if( cfg%hp > 0. .and. real(bd%shell(q)) < real(box) * smpd / cfg%hp ) sig2(q) = 1.e30
                        enddo
                    endif
                    if( cfg%pscale /= 0 )then   ! per-particle contrast scale against the other half's consensus, folded into the CTF weights (ratio: C persists across bands)
                        !$omp parallel do default(shared) private(i,ih,q,accn,accd,P)
                        do i = 1,n
                            ih = merge(2, 1, dat%eo(i) == 0); call klem_project_cell(bd, kmax, Gc(:,:,:,ih), dat%R(:,:,i), cfg%conv, P); accn = 0.d0; accd = 0.d0
                            do q = 1,bd%nband
                                if( bd%shell(q) > kmax ) cycle
                                accn = accn + real(bd%wq(q),dp) * real(real(conjg(dat%C(q,i) * scl(q) * P(q)) * dat%F(q,i)),dp) / real(sig2(q),dp)
                                accd = accd + real(bd%wq(q),dp) * real(abs(dat%C(q,i) * scl(q) * P(q)),dp)**2 / real(sig2(q),dp)
                            enddo
                            pscn(i) = real(min(max(accn / max(accd, 1.d-30), 0.3d0), 3.d0)) * psc(i)
                        enddo
                        !$omp end parallel do
                        do i = 1,n
                            dat%C(:,i) = dat%C(:,i) * (pscn(i) / psc(i)); psc(i) = pscn(i)
                        enddo
                        write(logfhandle,'(A,I0,A,F5.2,A,F5.2,A,F5.2)') '>>> worker ', ip, ': per-particle contrast scale median ', median_r(psc), ', min-max ', minval(psc), ' - ', maxval(psc)
                    endif
                    if( cfg%pnoise /= 0 )then   ! per-particle noise factor against the other half's consensus, this part's median as the reference
                        !$omp parallel do default(shared) private(i,ih,q,acc,P)
                        do i = 1,n
                            ih = merge(2, 1, dat%eo(i) == 0); call klem_project_cell(bd, kmax, Gc(:,:,:,ih), dat%R(:,:,i), cfg%conv, P); acc = 0.d0
                            do q = 1,bd%nband
                                if( bd%shell(q) > kmax ) cycle
                                acc = acc + real(bd%wq(q),dp) * real(abs(dat%F(q,i) - dat%C(q,i) * scl(q) * P(q)),dp)**2 / real(sig2(q),dp)
                            enddo
                            pres(i) = real(acc)
                        enddo
                        !$omp end parallel do
                        pmed = median_r(pres)
                        do i = 1,n
                            pfac(i) = min(max(pres(i) / max(pmed, 1.e-30), 0.3), 3.)
                        enddo
                    endif
                case(REQ_MSTEP)
                    read(u) sigma; close(u)
                    call klem_lattice_kernel(cfg%n1, cfg%n2, sigma, Ksm, i1, i2)
                    !$omp parallel do default(shared) private(i,m,j)
                    do i = 1,n
                        do m = 1,ngrid
                            W(m,i) = 0.
                            do j = 1,ngrid
                                W(m,i) = W(m,i) + r(j,i) * Ksm(j,m)
                            enddo
                        enddo
                        do m = ngrid+1,ncl
                            W(m,i) = r(m,i)
                        enddo
                    enddo
                    !$omp end parallel do
                    num = cmplx(0.d0,0.d0,dp); den = 0.d0
                    call klem_mstep_accum(bd, kmax, dat, W, ncl, cfg%conv, num, den)
                    call open_rep(u)
                    do ih = 1,2
                        do ic = 1,ncl
                            call pack_cd(num(:,:,:,ic,ih), den(:,:,:,ic,ih)); pnum4 = cmplx(pnum); write(u) pnum4; pden4 = real(pden); write(u) pden4
                        enddo
                    enddo
                    close(u); call done_rep()
                case(REQ_ESTEP)
                    do ih = 1,2
                        do ic = 1,ncl
                            read(u) pnum4; pnum = pnum4; call unpack_c(pnum, V(:,:,:,ic,ih))
                        enddo
                    enddo
                    read(u) sig2; read(u) scl; read(u) dopose; close(u)
                    call klem_estep_ll(bd, kmax, dat, V, ncl, cfg%conv, sig2, scl, ll)
                    if( cfg%pnoise /= 0 )then
                        !$omp parallel do default(shared) private(i)
                        do i = 1,n
                            ll(:,i) = ll(:,i) / pfac(i)
                        enddo
                        !$omp end parallel do
                    endif
                    !$omp parallel do default(shared) private(i,ic,mx,s)
                    do i = 1,n
                        mx = maxval(ll(:,i)); s = 0.
                        do ic = 1,ncl
                            r(ic,i) = exp(ll(ic,i) - mx); s = s + r(ic,i)
                        enddo
                        r(:,i) = r(:,i) / s; ll(:,i) = ll(:,i) - mx
                    enddo
                    !$omp end parallel do
                    if( cfg%balance > 0. .and. (cfg%balance_until <= 0 .or. it <= cfg%balance_until) ) call klem_sinkhorn(ll, ncl, ngrid, n, cfg%balance, r)
                    nmoved = 0; pgain = 0.
                    do i = 1,n
                        am(i) = maxloc(r(:,i), 1)
                    enddo
                    nchg = count(am /= am_prev); am_prev = am
                    if( dopose /= 0 ) call klem_pose_step(bd, kmax, dat, V, ncl, cfg%conv, sig2, scl, am, cfg%pose_step, cfg%pose_prior, nmoved, pgain)
                    occ = sum(r, 2)
                    call open_rep(u); pden(1:ncl) = real(occ, dp); pden4(1:ncl) = real(pden(1:ncl)); write(u) pden4(1:ncl); write(u) nmoved; write(u) pgain; write(u) nchg; close(u); call done_rep()
                case(REQ_NEWBOARD)   ! same particles, a new board: new config, new distr dir, new start positions
                    read(u) nb_infile, nb_ddir; close(u)
                    call klem_read_cfg(trim(nb_infile), cfg)
                    if( cfg%n1 * cfg%n2 + cfg%nfree /= ncl ) THROW_HARD('worker: the new board has a different lattice or free components')
                    ddir = nb_ddir
                    open(newunit=u, file=trim(ddir)//'/part_'//int2str_pad(ip,2)//'.txt', status='old', action='read')
                    do i = 1,n
                        read(u,*) m, j, g1(i), g2(i)
                        if( m /= pinds(i) ) THROW_HARD('worker: particle order differs between boards')
                    enddo
                    close(u)
                    call init_resp(); pfac = 1.; psc = 1.; am_prev = 0
                    write(logfhandle,'(A,I0,A)') '>>> worker ', ip, ': new board '//trim(nb_infile); flush(logfhandle)
                    cycle
                case(REQ_FINISH)
                    close(u); write(logfhandle,'(A,I0,A)') '>>> worker ', ip, ': finished'; flush(logfhandle)
                    exit
                case(REQ_FLUSH)
                    close(u)
                    open(newunit=io, file=trim(ddir)//'/part_'//int2str_pad(ip,2)//'_latent.txt', status='replace', action='write')
                    do i = 1,n
                        s = sum(r(1:ngrid,i)); tmp = 0.; mx = 0.
                        do m = 1,ngrid
                            tmp = tmp + r(m,i) * real(i1(m)) / max(s,1.e-12); mx = mx + r(m,i) * real(i2(m)) / max(s,1.e-12)
                        enddo
                        write(io,'(I0,1X,F8.4,1X,F8.4,1X,I0,1X,F7.4,1X,F7.4)') pinds(i), tmp, mx, maxloc(r(:,i), 1), maxval(r(:,i)), s
                    enddo
                    close(io)
                    open(newunit=io, file=trim(ddir)//'/part_'//int2str_pad(ip,2)//'_resp.bin', access='stream', form='unformatted', status='replace', action='write')
                    do i = 1,n
                        write(io) r(:,i)
                    enddo
                    close(io)
                    if( cfg%pose_every > 0 )then
                        open(newunit=io, file=trim(ddir)//'/part_'//int2str_pad(ip,2)//'_poses.txt', status='replace', action='write')
                        do i = 1,n
                            write(io,'(I0,1X,F9.4,1X,F9.4,1X,F9.4)') pinds(i), klem_rotmat2euler(dat%R(:,:,i))
                        enddo
                        close(io)
                    endif
                    call done_rep()
                    write(logfhandle,'(A,I0,A)') '>>> worker ', ip, ': board outputs written'; flush(logfhandle)
                case default
                    close(u); THROW_HARD('unknown request type')
            end select
        enddo

      contains

        subroutine init_resp()   ! responsibilities from the lattice start positions (g1, g2), a small share on the free components
            integer :: i_, m_
            real    :: s_
            do i_ = 1,n
                s_ = 0.
                do m_ = 1,ngrid
                    r(m_,i_) = exp(-0.5 * ((g1(i_) - real(i1(m_)))**2 + (g2(i_) - real(i2(m_)))**2) / 0.8**2); s_ = s_ + r(m_,i_)
                enddo
                do m_ = ngrid+1,ncl
                    r(m_,i_) = 0.02 * s_ / real(max(cfg%nfree,1))
                enddo
                r(:,i_) = r(:,i_) / sum(r(:,i_))
            enddo
        end subroutine init_resp

        subroutine open_rep( u_ )
            integer, intent(out) :: u_
            open(newunit=u_, file=trim(ddir)//'/rep_'//int2str_pad(ip,2)//'_'//int2str_pad(ireq,4)//'.bin', access='stream', form='unformatted', status='replace', action='write')
        end subroutine open_rep

        subroutine done_rep()
            call touch_file(trim(ddir)//'/rep_'//int2str_pad(ip,2)//'_'//int2str_pad(ireq,4)//'.done')
        end subroutine done_rep

        subroutine pack_cd( A_, D_ )
            complex(dp), intent(in) :: A_(-kmax:kmax,-kmax:kmax,-kmax:kmax)
            real(dp),    intent(in) :: D_(-kmax:kmax,-kmax:kmax,-kmax:kmax)
            integer :: m_
            do m_ = 1,nhs
                pnum(m_) = A_(hh(m_),kk(m_),lh(m_)); pden(m_) = D_(hh(m_),kk(m_),lh(m_))
            enddo
        end subroutine pack_cd

        subroutine unpack_c( pn_, C_ )
            complex(dp), intent(in)  :: pn_(nhs)
            complex,     intent(out) :: C_(-kmax:kmax,-kmax:kmax,-kmax:kmax)
            integer :: m_
            do m_ = 1,nhs
                C_(hh(m_),kk(m_),lh(m_)) = cmplx(pn_(m_))
                if( hh(m_) > 0 ) C_(-hh(m_),-kk(m_),-lh(m_)) = conjg(cmplx(pn_(m_)))
            enddo
        end subroutine unpack_c

    end subroutine klem_board_worker

    ! ------------------------------------------------------------------ small helpers (as in the sequential board)
    function median_dp( v ) result( med )
        real(dp), intent(in) :: v(:)
        real(dp) :: med
        real(dp), allocatable :: w_(:)
        integer :: n_
        allocate(w_(size(v))); w_ = v; call hpsort_dp(w_); n_ = size(w_)
        if( mod(n_,2) == 1 )then
            med = w_((n_+1)/2)
        else
            med = 0.5d0 * (w_(n_/2) + w_(n_/2+1))
        endif
        deallocate(w_)
    end function median_dp

    function median_r( v ) result( med )
        real, intent(in) :: v(:)
        real :: med
        real(dp), allocatable :: w_(:)
        allocate(w_(size(v))); w_ = real(v, dp); med = real(median_dp(w_)); deallocate(w_)
    end function median_r

    subroutine hpsort_dp( ra )   ! heapsort, Numerical Recipes form
        real(dp), intent(inout) :: ra(:)
        integer  :: n_, i_, ir_, j_, l_
        real(dp) :: rra
        n_ = size(ra)
        if( n_ < 2 ) return
        l_ = n_/2 + 1; ir_ = n_
        do
            if( l_ > 1 )then
                l_ = l_ - 1; rra = ra(l_)
            else
                rra = ra(ir_); ra(ir_) = ra(1); ir_ = ir_ - 1
                if( ir_ == 1 )then
                    ra(1) = rra; return
                endif
            endif
            i_ = l_; j_ = l_ + l_
            do while( j_ <= ir_ )
                if( j_ < ir_ )then
                    if( ra(j_) < ra(j_+1) ) j_ = j_ + 1
                endif
                if( rra < ra(j_) )then
                    ra(i_) = ra(j_); i_ = j_; j_ = j_ + j_
                else
                    j_ = ir_ + 1
                endif
            enddo
            ra(i_) = rra
        enddo
    end subroutine hpsort_dp

    function wtime() result( t )
        real(dp) :: t
        integer(8) :: cnt_, rate_
        call system_clock(cnt_, rate_); t = real(cnt_,dp) / real(max(rate_,1_8),dp)
    end function wtime

end module simple_klem_distr
