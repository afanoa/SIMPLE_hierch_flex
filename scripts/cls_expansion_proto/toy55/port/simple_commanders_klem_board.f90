!@descr: kernel-latent EM board (klem): a lattice of free cell volumes plus free junk components fitted to a particle
!! set at fixed poses by cross-half EM over a band schedule; the Fortran twin of klem.py's board loop (no pose
!! marginalisation, no sub-class modes). Driven by a config file (infile, "key value" lines):
!!   n1 7 | n2 7 | nfree 2 | lps 12 12 10 10 8 8 8 8 8 | sigmas 1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7 | shrink 0.02 |
!!   eps0 0.05 | balance 0.5 | balance_until 0 | focus_mask none|FILE | subset none|FILE (1-based project rows) |
!!   nsub 0 | seed 1 | ss -1 | fs 1 | conv 0 | sign 1 | ctfabs 1 | write_halves 0 | pose_every 0 | pose_step 3 | pose_prior 4 |
!!   devw 0 (0.05 = deviation-FSC shell weights) | pnoise 0|1 (per-particle noise factor) | init_file none|FILE (pind g1 g2)
!! Outputs in dir: klem_latent.txt, klem_states.txt, klem_resp.bin, klem_halves.bin, cells/cell_XX_YY.mrc (+ cells_half/).
module simple_commanders_klem_board
use simple_core_module_api
use simple_commander_base, only: commander_base
use simple_cmdline,        only: cmdline
use simple_parameters,     only: parameters
use simple_sp_project,     only: sp_project
use simple_image,          only: image
use simple_string,         only: string
use simple_klem_core
use simple_klem_distr,     only: klem_cfg, klem_read_cfg, klem_board_master, klem_board_worker, apply_focus_d
!$ use omp_lib
implicit none
public :: commander_klem_board, commander_klem_board_part
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_klem_board
  contains
    procedure :: execute => exec_klem_board
end type commander_klem_board

type, extends(commander_base) :: commander_klem_board_part
  contains
    procedure :: execute => exec_klem_board_part
end type commander_klem_board_part

contains

    subroutine exec_klem_board( self, cline )
        class(commander_klem_board), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(parameters)   :: params
        type(sp_project)   :: spproj
        type(klem_band)    :: bd
        type(klem_data)    :: dat
        type(klem_cfg)     :: cfg
        type(image)        :: mimg, dimg
        type(string)       :: fn
        character(len=512) :: line, key, cfg_focus, cfg_solv, cfg_subset, odir, fname
        character(len=512) :: tok(64)
        real,    allocatable :: lps(:), sigmas(:), r(:,:), ll(:,:), W(:,:), Ksm(:,:), sig2(:), scl(:), fm(:,:,:), rsub(:)
        real,    allocatable :: sc_sh(:), sig2_sh(:), num_sh(:), den_sh(:), occ(:), z1(:), z2(:), crep(:), g1(:), g2(:)
        complex, allocatable :: V(:,:,:,:,:), Gc(:,:,:,:), P(:), Vprev(:,:,:,:,:)
        real(dp) :: dnum, dden
        real     :: vchg
        complex(dp), allocatable :: num(:,:,:,:,:)
        real(dp),    allocatable :: den(:,:,:,:,:), dpos(:)
        integer, allocatable :: pinds(:), eos(:), i1(:), i2(:), am(:), stsel(:), subrows(:), idx(:)
        real(dp) :: tsum, t0, eps0med, acc, accn, accd
        real     :: shrink, eps0f, balance, sshift, fscale, csign, cur_lp, lp, sigma, eps, eps0, dmax, smpd, mx, s, tmp
        integer  :: n1, n2, nfree, ngrid, ncl, nit, seed, nsub, conv, balance_until, write_halves, ntok, io, u, i, j, m, ic, it, pose_every, nmoved
        real     :: pose_step, pose_prior, pgain, devw, fsc_s, pmed
        integer  :: pnoise, pscale, nchg
        real     :: stop_frac
        integer, allocatable :: am_prev(:)
        character(len=512) :: cfg_init, cfg_boards
        type(klem_cfg), allocatable :: cfgs_more(:)
        character(len=512), allocatable :: odirs_more(:)
        real, allocatable :: g1_more(:,:), g2_more(:,:)
        integer :: nmore, kb
        real,    allocatable :: sig2_base(:), wsh(:), pfac(:), pres(:), gi1(:), gi2(:), psc(:), pscn(:)
        real(dp), allocatable :: fnum(:), fd0(:), fd1(:)
        complex  :: d0c, d1c
        integer  :: box, kmax, n, nall, q, ih, h, k, l, ncal, nsh, ish, cnt, n3, ib, nbm, c2, ctfabs_i
        logical  :: ctfabs, l_focus, l_solv, l_exists
        real,    allocatable :: sm(:,:,:)
        complex, allocatable :: G0(:,:,:,:)
        ! ---- defaults and config
        n1 = 7; n2 = 7; nfree = 2; shrink = 0.02; eps0f = 0.05; balance = 0.5; balance_until = 0; cfg_focus = 'none'; cfg_solv = 'none'; cfg_subset = 'none'
        nsub = 0; seed = 1; sshift = -1.; fscale = 1.; conv = 0; csign = 1.; ctfabs_i = 1; write_halves = 0; pose_every = 0; pose_step = 3.; pose_prior = 4.; devw = 0.; pnoise = 0; pscale = 0; stop_frac = 0.03; cfg_init = 'none'; cfg_boards = 'none'
        allocate(lps(9), sigmas(9)); lps = [12.,12.,10.,10.,8.,8.,8.,8.,8.]; sigmas = [1.2,1.2,1.0,1.0,0.8,0.8,0.7,0.7,0.7]
        call params%new(cline)
        open(newunit=u, file=params%infile%to_char(), status='old', action='read')
        do
            read(u,'(A)',iostat=io) line
            if( io /= 0 ) exit
            if( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
            call tokenize(line, tok, ntok)
            if( ntok < 2 ) cycle
            key = tok(1)
            select case(trim(key))
                case('n1');            read(tok(2),*) n1
                case('n2');            read(tok(2),*) n2
                case('nfree');         read(tok(2),*) nfree
                case('shrink');        read(tok(2),*) shrink
                case('eps0');          read(tok(2),*) eps0f
                case('balance');       read(tok(2),*) balance
                case('balance_until'); read(tok(2),*) balance_until
                case('focus_mask');    cfg_focus = tok(2)
                case('solvent_mask');  cfg_solv = tok(2)
                case('subset');        cfg_subset = tok(2)
                case('nsub');          read(tok(2),*) nsub
                case('seed');          read(tok(2),*) seed
                case('ss');            read(tok(2),*) sshift
                case('fs');            read(tok(2),*) fscale
                case('conv');          read(tok(2),*) conv
                case('sign');          read(tok(2),*) csign
                case('ctfabs');        read(tok(2),*) ctfabs_i
                case('write_halves');  read(tok(2),*) write_halves
                case('pose_every');    read(tok(2),*) pose_every
                case('pose_step');     read(tok(2),*) pose_step
                case('pose_prior');    read(tok(2),*) pose_prior
                case('devw');          read(tok(2),*) devw
                case('pnoise');        read(tok(2),*) pnoise
                case('pscale');        read(tok(2),*) pscale
                case('stop_frac');     read(tok(2),*) stop_frac
                case('init_file');     cfg_init = tok(2)
                case('boards_more');   cfg_boards = tok(2)
                case('lps')
                    deallocate(lps); allocate(lps(ntok-1))
                    do j = 2,ntok
                        read(tok(j),*) lps(j-1)
                    enddo
                case('sigmas')
                    deallocate(sigmas); allocate(sigmas(ntok-1))
                    do j = 2,ntok
                        read(tok(j),*) sigmas(j-1)
                    enddo
            end select
        enddo
        close(u)
        if( size(lps) /= size(sigmas) ) THROW_HARD('lps and sigmas must have the same length')
        nit = size(lps); ctfabs = ctfabs_i /= 0; l_focus = trim(cfg_focus) /= 'none'; l_solv = trim(cfg_solv) /= 'none'; odir = trim(params%dir%to_char())
        call simple_mkdir(trim(odir)); call simple_mkdir(trim(odir)//'/cells')
        if( write_halves /= 0 ) call simple_mkdir(trim(odir)//'/cells_half')
        ! ---- project, particle set, halves
        call spproj%read(params%projfile)
        box = params%box; smpd = params%smpd; nall = spproj%os_ptcl3D%get_noris()
        allocate(stsel(nall), source=0)
        do i = 1,nall
            if( spproj%os_ptcl3D%get_state(i) > 0 ) stsel(i) = 1
        enddo
        if( trim(cfg_subset) /= 'none' )then
            allocate(subrows(nall), source=0)
            open(newunit=u, file=trim(cfg_subset), status='old', action='read')
            do
                read(u,*,iostat=io) i
                if( io /= 0 ) exit
                if( i >= 1 .and. i <= nall ) subrows(i) = 1
            enddo
            close(u)
            stsel = stsel * subrows; deallocate(subrows)
        endif
        n = sum(stsel); allocate(pinds(n)); j = 0
        do i = 1,nall
            if( stsel(i) == 1 )then
                j = j + 1; pinds(j) = i
            endif
        enddo
        call seed_rnd_simple(seed)
        if( nsub > 0 .and. nsub < n )then
            allocate(rsub(n)); call random_number(rsub); allocate(idx(n))
            call argsort_real(rsub, idx); allocate(subrows(nsub)); subrows = pinds(idx(1:nsub)); call sort_int(subrows)
            deallocate(pinds); allocate(pinds(nsub)); pinds = subrows; n = nsub; deallocate(rsub, idx, subrows)
        endif
        allocate(eos(n)); cnt = 0
        do i = 1,n
            eos(i) = spproj%os_ptcl3D%get_eo(pinds(i))
            if( eos(i) < 0 .or. eos(i) > 1 ) eos(i) = mod(pinds(i), 2)
            cnt = cnt + eos(i)
        enddo
        if( cnt == 0 .or. cnt == n )then   ! no usable eo labels: split by row parity
            do i = 1,n
                eos(i) = mod(pinds(i), 2)
            enddo
        endif
        ngrid = n1 * n2; ncl = ngrid + nfree
        write(logfhandle,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') '>>> KLEM BOARD: ', n, ' particles, grid ', n1, ' x ', n2, ' + ', nfree, ' free components, ', nit, ' rounds'
        write(logfhandle,'(A,9F6.1)') '>>> bands (A): ', lps
        write(logfhandle,'(A,9F6.2)') '>>> kernel widths (cells): ', sigmas
        ! ---- initial responsibilities: random lattice positions, Gaussian width 0.8, 2 % on the free components
        allocate(Ksm(ngrid,ngrid), i1(ngrid), i2(ngrid)); call klem_lattice_kernel(n1, n2, 1., Ksm, i1, i2)
        allocate(r(ncl,n), source=0.); allocate(g1(n), g2(n)); call random_number(g1); call random_number(g2); g1 = g1 * real(n1 - 1); g2 = g2 * real(n2 - 1)
        if( trim(cfg_init) /= 'none' )then   ! lattice positions from a file "pind g1 g2" (the prototype's npz/latent init, made by klem_init_positions.py)
            allocate(gi1(nall), gi2(nall), source=-1.)
            open(newunit=u, file=trim(cfg_init), status='old', action='read')
            do
                read(u,*,iostat=io) i, tmp, s
                if( io /= 0 ) exit
                if( i >= 1 .and. i <= nall )then
                    gi1(i) = tmp; gi2(i) = s
                endif
            enddo
            close(u); cnt = 0
            do i = 1,n
                if( gi1(pinds(i)) >= 0. )then
                    g1(i) = min(max(gi1(pinds(i)), 0.), real(n1 - 1)); g2(i) = min(max(gi2(pinds(i)), 0.), real(n2 - 1)); cnt = cnt + 1
                endif
            enddo
            write(logfhandle,'(A,I0,A,I0,A)') '>>> init: lattice positions from '//trim(cfg_init)//' for ', cnt, ' of ', n, ' particles (the rest random)'
            deallocate(gi1, gi2)
        endif
        if( params%nparts > 1 )then   ! distributed board: persistent workers hold the particles, this process merges
            call klem_read_cfg(params%infile%to_char(), cfg)
            if( trim(cfg_boards) /= 'none' )then   ! board server: the listed board dirs (each with board.txt + its init file) follow on the same workers
                nmore = 0
                open(newunit=u, file=trim(cfg_boards), status='old', action='read')
                do
                    read(u,'(A)',iostat=io) line
                    if( io /= 0 ) exit
                    if( len_trim(line) > 0 ) nmore = nmore + 1
                enddo
                rewind(u); allocate(cfgs_more(nmore), odirs_more(nmore), g1_more(n,nmore), g2_more(n,nmore)); kb = 0
                do
                    read(u,'(A)',iostat=io) line
                    if( io /= 0 ) exit
                    if( len_trim(line) == 0 ) cycle
                    kb = kb + 1; odirs_more(kb) = trim(adjustl(line)); call klem_read_cfg(trim(odirs_more(kb))//'/board.txt', cfgs_more(kb))
                    call read_init_positions(trim(cfgs_more(kb)%init), cfgs_more(kb)%n1, cfgs_more(kb)%n2, g1_more(:,kb), g2_more(:,kb))
                enddo
                close(u)
                write(logfhandle,'(A,I0,A)') '>>> board server: ', nmore, ' further boards on the same workers'
                call klem_board_master(params, spproj, cfg, pinds, eos, g1, g2, trim(odir), box, smpd, cfgs_more, odirs_more, g1_more, g2_more)
            else
                call klem_board_master(params, spproj, cfg, pinds, eos, g1, g2, trim(odir), box, smpd)
            endif
            call spproj%kill
            call simple_end('**** SIMPLE_KLEM_BOARD NORMAL STOP ****')
            return
        endif
        do i = 1,n
            s = 0.
            do m = 1,ngrid
                r(m,i) = exp(-0.5 * ((g1(i) - real(i1(m)))**2 + (g2(i) - real(i2(m)))**2) / 0.8**2); s = s + r(m,i)
            enddo
            do m = ngrid+1,ncl
                r(m,i) = 0.02 * s / real(max(nfree,1))
            enddo
            r(:,i) = r(:,i) / sum(r(:,i))
        enddo
        deallocate(g1, g2)
        ! ---- EM over the schedule
        cur_lp = -1.; allocate(occ(ncl), am(n), z1(n), z2(n), am_prev(n)); am_prev = 0; it = 0
        do
            it = it + 1; if( it > nit ) exit
            lp = lps(it); sigma = sigmas(it); t0 = wtime()
            if( abs(lp - cur_lp) > 1.e-6 )then
                kmax = floor(real(box) * smpd / lp); n3 = 2 * kmax + 1
                if( .not. allocated(dat%F) )then   ! load once at the finest band of the schedule; coarser bands use the coefficients within their kmax
                    call klem_band_new(bd, floor(real(box) * smpd / minval(lps)))
                    write(logfhandle,'(A,F5.1,A,I0,A,I0,A)') '>>> loading the particles once at the finest band ', minval(lps), ' A: kmax ', bd%kmax, ', ', bd%nband, ' coefficients'
                    flush(logfhandle)
                    call klem_load_particles(spproj, bd, pinds, eos, box, smpd, sshift, fscale, csign, ctfabs, dat)
                endif
                write(logfhandle,'(A,F5.1,A,I0,A,I0,A)') '>>> band ', lp, ' A: kmax ', kmax, ', grid ', n3, '^3'
                flush(logfhandle)
                ! consensus of each half (all particles, weight 1) -> shrinkage target
                if( allocated(num) ) deallocate(num, den, Gc, V, ll, W, sig2, scl)
                allocate(num(-kmax:kmax,-kmax:kmax,-kmax:kmax,1,2), source=cmplx(0.d0,0.d0,dp)); allocate(den(-kmax:kmax,-kmax:kmax,-kmax:kmax,1,2), source=0.d0)
                allocate(W(1,n), source=1.); call klem_mstep_accum(bd, kmax, dat, W, 1, conv, num, den); deallocate(W)
                allocate(dpos(2 * n3**3)); cnt = 0
                do ih = 1,2
                    do l = -kmax,kmax
                        do k = -kmax,kmax
                            do h = -kmax,kmax
                                if( den(h,k,l,1,ih) > 0.d0 )then
                                    cnt = cnt + 1; dpos(cnt) = den(h,k,l,1,ih)
                                endif
                            enddo
                        enddo
                    enddo
                enddo
                eps0med = median_dp(dpos(1:cnt)); deallocate(dpos); eps0 = real(eps0f * eps0med); eps = real(shrink * eps0med)
                allocate(Gc(-kmax:kmax,-kmax:kmax,-kmax:kmax,2))
                do ih = 1,2
                    Gc(:,:,:,ih) = cmplx(num(:,:,:,1,ih) / (den(:,:,:,1,ih) + real(eps0,dp)))
                enddo
                deallocate(num, den)
                allocate(num(-kmax:kmax,-kmax:kmax,-kmax:kmax,ncl,2), source=cmplx(0.d0,0.d0,dp)); allocate(den(-kmax:kmax,-kmax:kmax,-kmax:kmax,ncl,2), source=0.d0)
                allocate(V(-kmax:kmax,-kmax:kmax,-kmax:kmax,ncl,2), source=cmplx(0.,0.)); allocate(ll(ncl,n), source=0.); allocate(W(ncl,n), source=0.)
                ! calibration on up to 2000 particles: reference scale and noise per shell against the other half's consensus
                nsh = kmax + 1; allocate(sig2(bd%nband), scl(bd%nband), sc_sh(0:kmax), sig2_sh(0:kmax), num_sh(0:kmax), den_sh(0:kmax), P(bd%nband))
                num_sh = 0.; den_sh = 0.; ncal = min(2000, n)
                do i = 1,n,max(1, n / ncal)
                    ih = merge(2, 1, dat%eo(i) == 0); call klem_project_cell(bd, kmax, Gc(:,:,:,ih), dat%R(:,:,i), conv, P)
                    do q = 1,bd%nband
                        if( bd%shell(q) > kmax ) cycle
                        tmp = dat%C(q,i); ish = bd%shell(q)
                        num_sh(ish) = num_sh(ish) + bd%wq(q) * real(conjg(tmp * P(q)) * dat%F(q,i)); den_sh(ish) = den_sh(ish) + bd%wq(q) * abs(tmp * P(q))**2
                    enddo
                enddo
                do ish = 0,kmax
                    sc_sh(ish) = merge(num_sh(ish) / max(den_sh(ish), 1.e-30), 0., den_sh(ish) > 0.)
                enddo
                ! 3-tap smoothing with edge padding, then the floor at 5 % of the mean of shells 1-4
                num_sh = sc_sh
                do ish = 0,kmax
                    sc_sh(ish) = (num_sh(max(ish-1,0)) + num_sh(ish) + num_sh(min(ish+1,kmax))) / 3.
                enddo
                tmp = 0.05 * max(sum(sc_sh(1:min(4,kmax))) / real(min(4,kmax)), 1.e-30); sc_sh = max(sc_sh, tmp)
                do q = 1,bd%nband
                    scl(q) = 1.
                    if( bd%shell(q) <= kmax ) scl(q) = sc_sh(bd%shell(q))
                enddo
                sig2_sh = 0.; num_sh = 0.
                do i = 1,n,max(1, n / ncal)
                    ih = merge(2, 1, dat%eo(i) == 0); call klem_project_cell(bd, kmax, Gc(:,:,:,ih), dat%R(:,:,i), conv, P)
                    do q = 1,bd%nband
                        if( bd%shell(q) > kmax ) cycle
                        ish = bd%shell(q); sig2_sh(ish) = sig2_sh(ish) + abs(dat%F(q,i) - dat%C(q,i) * scl(q) * P(q))**2; num_sh(ish) = num_sh(ish) + 1.
                    enddo
                enddo
                do ish = 0,kmax
                    sig2_sh(ish) = max(merge(sig2_sh(ish) / num_sh(ish), 1., num_sh(ish) > 0.), 1.e-12)
                enddo
                do q = 1,bd%nband
                    sig2(q) = 1.
                    if( bd%shell(q) <= kmax ) sig2(q) = sig2_sh(bd%shell(q))
                enddo
                if( allocated(sig2_base) ) deallocate(sig2_base, wsh)
                allocate(sig2_base(bd%nband), wsh(0:kmax)); sig2_base = sig2; wsh = 1.
                if( pscale /= 0 )then   ! per-particle contrast scale against the other half's consensus, folded into the CTF weights (ratio to the previous scale: C persists)
                    if( .not. allocated(psc) )then
                        allocate(psc(n), pscn(n)); psc = 1.
                    endif
                    !$omp parallel do default(shared) private(i,ih,q,accn,accd,P)
                    do i = 1,n
                        ih = merge(2, 1, dat%eo(i) == 0); call klem_project_cell(bd, kmax, Gc(:,:,:,ih), dat%R(:,:,i), conv, P); accn = 0.d0; accd = 0.d0
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
                    write(logfhandle,'(A,F5.2,A,F5.2,A,F5.2)') '>>> per-particle contrast scale: median ', median_r(psc), ', 5-95 pct ', percentile_r(psc, 5.), ' - ', percentile_r(psc, 95.)
                endif
                if( pnoise /= 0 )then   ! per-particle noise factor: residual against the other half's consensus relative to the median
                    if( .not. allocated(pfac) ) allocate(pfac(n), pres(n))
                    !$omp parallel do default(shared) private(i,ih,q,acc,P)
                    do i = 1,n   ! P is private and already allocated (private copies inherit the allocation)
                        ih = merge(2, 1, dat%eo(i) == 0); call klem_project_cell(bd, kmax, Gc(:,:,:,ih), dat%R(:,:,i), conv, P); acc = 0.d0
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
                    write(logfhandle,'(A,F5.2,A,F5.2)') '>>> per-particle noise factor: 5-95 pct ', percentile_r(pfac, 5.), ' - ', percentile_r(pfac, 95.)
                endif
                write(logfhandle,'(A,F5.1,A,ES9.2)') '>>> band ', lp, ' A: shrinkage eps ', eps
                write(logfhandle,'(A,*(F6.2,1X))') '>>> reference scale relative to shell 1 (every few shells): ', (sc_sh(min(ish,kmax)) / max(sc_sh(1),1e-30), ish=1,kmax,max(1,kmax/6))
                deallocate(sc_sh, sig2_sh, num_sh, den_sh, P)
                if( l_focus )then   ! the focal mask on the full box (the deviation is masked in real space at the full even box)
                    if( allocated(fm) ) deallocate(fm)
                    fn = trim(cfg_focus); call mimg%new([box,box,box], smpd); call mimg%read(fn); allocate(fm(box,box,box)); call mimg%get_rmat_sub(fm); call mimg%kill
                    fm = min(max(fm, 0.), 1.)
                    write(logfhandle,'(A,F6.2,A)') '>>> focal mask: ', 100. * real(count(fm > 0.5)) / real(box)**3, ' % of the box above 0.5 (references differ only inside it)'
                endif
                if( l_solv )then
                    if( allocated(sm) ) deallocate(sm)
                    fn = trim(cfg_solv); call mimg%new([box,box,box], smpd); call mimg%read(fn); allocate(sm(box,box,box)); call mimg%get_rmat_sub(sm); call mimg%kill
                    sm = min(max(sm, 0.), 1.)
                endif
                cur_lp = lp
            endif
            ! ---- M-step: kernel-smoothed weights, Wiener-shrunk cells, free components unsmoothed
            call klem_lattice_kernel(n1, n2, sigma, Ksm, i1, i2)
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
            call klem_mstep_accum(bd, kmax, dat, W, ncl, conv, num, den)
            do ih = 1,2
                do ic = 1,ncl
                    V(:,:,:,ic,ih) = cmplx((num(:,:,:,ic,ih) + real(eps,dp) * Gc(:,:,:,ih)) / (den(:,:,:,ic,ih) + real(eps,dp)))
                enddo
            enddo
            if( l_focus ) call apply_focus()
            if( l_solv )then
                allocate(G0(-kmax:kmax,-kmax:kmax,-kmax:kmax,2), source=cmplx(0.,0.)); call apply_focus_d(V, G0, sm, kmax, ncl, box, smpd); deallocate(G0)
            endif
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
            if( devw > 0. .and. it > 1 )then   ! pooled cross-half FSC of the grid-cell deviations per shell -> likelihood shell weights 2F/(1+F)
                allocate(fnum(0:kmax), fd0(0:kmax), fd1(0:kmax)); fnum = 0.d0; fd0 = 0.d0; fd1 = 0.d0
                do ic = 1,ngrid
                    do l = -kmax,kmax
                        do k = -kmax,kmax
                            do h = -kmax,kmax
                                ish = nint(sqrt(real(h*h + k*k + l*l)))
                                if( ish > kmax ) cycle
                                d0c = V(h,k,l,ic,1) - Gc(h,k,l,1); d1c = V(h,k,l,ic,2) - Gc(h,k,l,2)
                                fnum(ish) = fnum(ish) + real(conjg(d0c) * d1c, dp); fd0(ish) = fd0(ish) + real(abs(d0c),dp)**2; fd1(ish) = fd1(ish) + real(abs(d1c),dp)**2
                            enddo
                        enddo
                    enddo
                enddo
                do ish = 0,kmax
                    fsc_s = real(fnum(ish) / max(sqrt(fd0(ish) * fd1(ish)), 1.d-30)); wsh(ish) = min(max(2. * fsc_s / (1. + max(fsc_s, 1.e-6)), devw), 1.)
                enddo
                wsh(0:min(1,kmax)) = 1.
                do q = 1,bd%nband
                    if( bd%shell(q) <= kmax ) sig2(q) = sig2_base(q) / wsh(bd%shell(q))
                enddo
                write(logfhandle,'(A,*(F5.2,1X))') '>>> deviation-FSC likelihood weights (every few shells): ', (wsh(min(ish,kmax)), ish=1,kmax,max(1,kmax/6))
                deallocate(fnum, fd0, fd1)
            endif
            ! ---- E-step (cross-half), responsibilities, balance
            call klem_estep_ll(bd, kmax, dat, V, ncl, conv, sig2, scl, ll)
            if( pnoise /= 0 )then
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
            if( balance > 0. .and. (balance_until <= 0 .or. it <= balance_until) ) call klem_sinkhorn(ll, ncl, ngrid, n, balance, r)
            if( pose_every > 0 )then
                if( mod(it, pose_every) == 0 )then
                    do i = 1,n
                        am(i) = maxloc(r(:,i), 1)
                    enddo
                    call klem_pose_step(bd, kmax, dat, V, ncl, conv, sig2, scl, am, pose_step, pose_prior, nmoved, pgain)
                    write(logfhandle,'(A,I0,A,F5.1,A,F6.3,A,F8.2)') '>>> pose step (it ', it, ', step ', pose_step, ' deg): fraction moved ', real(nmoved) / real(n), ', mean log-likelihood gain ', pgain
                endif
            endif
            occ = sum(r, 2) / real(n)
            do i = 1,n
                am(i) = maxloc(r(:,i), 1)
            enddo
            nchg = count(am /= am_prev); am_prev = am
            write(logfhandle,'(A,I0,A,F5.1,A,F4.1,A,3F7.3,A,F6.3,A,F6.3,A,F6.3,A,F7.1,A)') '>>> it ', it, ' (lp ', lp, ', sigma ', sigma, '): grid occupancy min/median/max ', minval(occ(1:ngrid)), median_r(occ(1:ngrid)), maxval(occ(1:ngrid)), &
                &', free mass ', sum(occ(ngrid+1:ncl)), ', argmax changed ', real(nchg) / real(n), ', model change ', vchg, ' (', wtime() - t0, ' s)'
            flush(logfhandle)
            if( vchg < stop_frac )then   ! converged at this band: jump to its last iteration (final sigma), skipping the middle ones
                j = it
                do while( j < nit )
                    if( abs(lps(j+1) - lp) > 1.e-6 ) exit
                    j = j + 1
                enddo
                if( j - 1 > it )then
                    write(logfhandle,'(A,I0,A,I0)') '>>> early stop at this band: skipping iterations ', it + 1, ' to ', j - 1; it = j - 1
                endif
            endif
        enddo
        ! ---- readouts: posterior means on the lattice, argmax, reproducibility of every cell's deviation
        do i = 1,n
            s = sum(r(1:ngrid,i)); z1(i) = 0.; z2(i) = 0.
            do m = 1,ngrid
                z1(i) = z1(i) + r(m,i) * real(i1(m)) / max(s,1.e-12); z2(i) = z2(i) + r(m,i) * real(i2(m)) / max(s,1.e-12)
            enddo
            am(i) = maxloc(r(:,i), 1)
        enddo
        allocate(crep(ncl)); c2 = n3 / 2
        do ic = 1,ncl
            acc = 0.d0; tsum = 0.d0; eps0med = 0.d0
            do l = -kmax,kmax
                do k = -kmax,kmax
                    do h = -kmax,kmax
                        if( h*h + k*k + l*l < 1 .or. h*h + k*k + l*l > c2*c2 ) cycle
                        acc    = acc    + real(conjg(V(h,k,l,ic,1) - Gc(h,k,l,1)) * (V(h,k,l,ic,2) - Gc(h,k,l,2)), dp)
                        tsum   = tsum   + real(abs(V(h,k,l,ic,1) - Gc(h,k,l,1)),dp)**2
                        eps0med = eps0med + real(abs(V(h,k,l,ic,2) - Gc(h,k,l,2)),dp)**2
                    enddo
                enddo
            enddo
            crep(ic) = real(acc / max(sqrt(tsum * eps0med), 1.d-30))
        enddo
        write(logfhandle,'(A,F5.2,A,I0,A,I0,A)') '>>> states: cross-half reproducibility of the deviation median ', median_r(crep(1:ngrid)), ', ', count(crep(1:ngrid) > 0.5), ' of ', ngrid, ' above 0.5'
        ! ---- outputs
        open(newunit=u, file=trim(odir)//'/klem_latent.txt', status='replace', action='write')
        write(u,'(A)') '# pind z1 z2 argmax maxr gridmass layer2mass'
        do i = 1,n
            write(u,'(I0,1X,F8.4,1X,F8.4,1X,I0,1X,F7.4,1X,F7.4,1X,F7.4)') pinds(i), z1(i), z2(i), am(i), maxval(r(:,i)), sum(r(1:ngrid,i)), 0.
        enddo
        close(u)
        open(newunit=u, file=trim(odir)//'/klem_states.txt', status='replace', action='write')
        write(u,'(A)') '# cell population reproducibility isgrid'
        do ic = 1,ncl
            write(u,'(I0,1X,I0,1X,F8.4,1X,I0)') ic, count(am == ic), crep(ic), merge(1, 0, ic <= ngrid)
        enddo
        close(u)
        open(newunit=u, file=trim(odir)//'/klem_resp.bin', access='stream', form='unformatted', status='replace', action='write')
        write(u) n, ncl, ngrid, n1, n2
        write(u) pinds
        write(u) i1, i2
        do i = 1,n
            write(u) r(:,i)
        enddo
        close(u)
        open(newunit=u, file=trim(odir)//'/klem_halves.bin', access='stream', form='unformatted', status='replace', action='write')
        write(u) ncl, ngrid, kmax, n3
        write(u) sum(r(1:ngrid,:), 2)
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
        do ic = 1,ncl
            if( ic <= ngrid )then
                fname = trim(odir)//'/cells/cell_'//int2str_pad(i1(ic),2)//'_'//int2str_pad(i2(ic),2)//'.mrc'
            else
                fname = trim(odir)//'/cells/free_'//int2str_pad(ic - ngrid, 1)//'.mrc'
            endif
            call klem_write_volume(0.5 * (V(:,:,:,ic,1) + V(:,:,:,ic,2)), kmax, box, smpd, trim(fname))
            if( write_halves /= 0 )then
                call klem_write_volume(V(:,:,:,ic,1), kmax, box, smpd, trim(odir)//'/cells_half/cell_'//int2str_pad(ic,2)//'_even.mrc')
                call klem_write_volume(V(:,:,:,ic,2), kmax, box, smpd, trim(odir)//'/cells_half/cell_'//int2str_pad(ic,2)//'_odd.mrc')
            endif
        enddo
        if( pose_every > 0 )then   ! the moved poses, SIMPLE Euler convention, one row per particle
            open(newunit=u, file=trim(odir)//'/klem_poses.txt', status='replace', action='write')
            write(u,'(A)') '# pind e1 e2 e3'
            do i = 1,n
                write(u,'(I0,1X,F9.4,1X,F9.4,1X,F9.4)') pinds(i), klem_rotmat2euler(dat%R(:,:,i))
            enddo
            close(u)
        endif
        write(logfhandle,'(A)') '>>> KLEM BOARD written to '//trim(odir)
        call spproj%kill
        call simple_end('**** SIMPLE_KLEM_BOARD NORMAL STOP ****')

      contains

        !> cells = consensus + mask x (cell - consensus), the mask applied in real space at the full even box
        subroutine apply_focus()   ! the distributed board's parallel version, one image per thread
            call apply_focus_d(V, Gc, fm, kmax, ncl, box, smpd)
        end subroutine apply_focus

        subroutine tokenize( str, toks, nt )
            character(len=*), intent(in)  :: str
            character(len=*), intent(out) :: toks(:)
            integer,          intent(out) :: nt
            integer :: p, q, lenr
            nt = 0; p = 1; lenr = len_trim(str)
            do while( p <= lenr )
                do while( p <= lenr .and. (str(p:p) == ' ' .or. str(p:p) == char(9)) )
                    p = p + 1
                enddo
                if( p > lenr ) exit
                q = p
                do while( q <= lenr .and. str(q:q) /= ' ' .and. str(q:q) /= char(9) )
                    q = q + 1
                enddo
                nt = nt + 1
                if( nt > size(toks) ) exit
                toks(nt) = str(p:q-1); p = q
            enddo
        end subroutine tokenize

        subroutine seed_rnd_simple( sd )
            integer, intent(in) :: sd
            integer, allocatable :: sarr(:)
            integer :: ns, j_
            call random_seed(size=ns); allocate(sarr(ns))
            do j_ = 1,ns
                sarr(j_) = sd * 7919 + j_ * 104729
            enddo
            call random_seed(put=sarr); deallocate(sarr)
        end subroutine seed_rnd_simple

        subroutine argsort_real( v, ord )
            real,    intent(in)  :: v(:)
            integer, intent(out) :: ord(:)
            integer :: a_, b_, t_
            do a_ = 1,size(v)
                ord(a_) = a_
            enddo
            do a_ = 2,size(v)   ! insertion sort (nsub draws only)
                t_ = ord(a_); b_ = a_ - 1
                do while( b_ >= 1 )
                    if( v(ord(b_)) <= v(t_) ) exit
                    ord(b_+1) = ord(b_); b_ = b_ - 1
                enddo
                ord(b_+1) = t_
            enddo
        end subroutine argsort_real

        subroutine sort_int( v )
            integer, intent(inout) :: v(:)
            integer :: a_, b_, t_
            do a_ = 2,size(v)
                t_ = v(a_); b_ = a_ - 1
                do while( b_ >= 1 )
                    if( v(b_) <= t_ ) exit
                    v(b_+1) = v(b_); b_ = b_ - 1
                enddo
                v(b_+1) = t_
            enddo
        end subroutine sort_int

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

        function percentile_r( v, pc ) result( val )
            real, intent(in) :: v(:), pc
            real :: val
            real(dp), allocatable :: w_(:)
            integer :: k_
            allocate(w_(size(v))); w_ = real(v, dp); call hpsort_dp(w_); k_ = max(1, min(size(v), nint(pc / 100. * real(size(v))))); val = real(w_(k_)); deallocate(w_)
        end function percentile_r

        function median_r( v ) result( med )
            real, intent(in) :: v(:)
            real :: med
            real(dp), allocatable :: w_(:)
            allocate(w_(size(v))); w_ = real(v, dp); med = real(median_dp(w_)); deallocate(w_)
        end function median_r

        subroutine read_init_positions( initfile, n1_, n2_, g1o, g2o )   ! lattice start positions "pind g1 g2" for the current pinds (random where absent)
            character(len=*), intent(in)  :: initfile
            integer,          intent(in)  :: n1_, n2_
            real,             intent(out) :: g1o(:), g2o(:)
            real, allocatable :: a1(:), a2(:)
            integer :: u_, io_, i_, c_
            real    :: t1, t2
            call random_number(g1o); call random_number(g2o); g1o = g1o * real(n1_ - 1); g2o = g2o * real(n2_ - 1)
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
            write(logfhandle,'(A,I0,A,I0,A)') '>>> init (server board): positions from '//trim(initfile)//' for ', c_, ' of ', n, ' particles'
            deallocate(a1, a2)
        end subroutine read_init_positions

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

    end subroutine exec_klem_board

    !> one persistent worker of the distributed board: loads its part of the particles once and answers the master's
    !! file requests (prg=klem_board_part projfile infile dir part nparts nthr mkdir=no)
    subroutine exec_klem_board_part( self, cline )
        class(commander_klem_board_part), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        type(parameters) :: params
        call params%new(cline)
        call klem_board_worker(params)
        call simple_end('**** SIMPLE_KLEM_BOARD_PART NORMAL STOP ****')
    end subroutine exec_klem_board_part

end module simple_commanders_klem_board
