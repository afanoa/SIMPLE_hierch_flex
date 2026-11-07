!@descr: cross-half split test of the states of a node (SSNR junk screen, docked equal-N pair maps against a stratified shuffle null, cross-half EM, re-test beyond the EM band) and the merge-up of states along a state tree
module simple_state_mergeup
!$ use omp_lib
use simple_defs,          only: dp, logfhandle, OSMPL_PAD_FAC, MAXIMGBATCHSZ, PI, nthr_glob
use simple_type_defs,     only: fplane_type
use simple_error,         only: simple_exception
use simple_string,        only: string
use simple_string_utils,  only: int2str, int2str_pad
use simple_fileio,        only: fopen, fclose, fileiochk, nlines
use simple_ori,           only: ori
use simple_sym,           only: sym
use simple_image,         only: image
use simple_reconstructor, only: reconstructor, exp_samples
use simple_parameters,    only: parameters
use simple_builder,       only: builder
use simple_state_tree,    only: state_tree
use simple_stat,          only: avg_ranks
use simple_linalg,        only: matinv
implicit none

public :: state_mergeup, merge_poles, merge_poles_subspace, gram_corr, gram_resid, fisher_z, ivec, EM_LIM
public :: HALVES_KSIG, mask_neff, seed_poles, split_zcrit, locality_guard, effect_floor, write_split_stats, read_split_stats
public :: SITE_GRID
private
#include "simple_local_flags.inc"

integer, parameter :: NSHUF = 20, NSHUF_BONF = 50, NVIEW = 40, NEM_ITERS = 3, MAXJOBS = 1000, RND_SEED = 0
real(dp), parameter :: ALPHA_SPLIT = 1._dp / 21._dp
integer,  parameter :: NSHUF_STEP = 10
real,     parameter :: NSHUF_KSIG = 4.
integer, parameter :: NDOCK_ITS = 8, RB_MARGIN = 2, RB_MARGIN_FOCUS = 4, SITE_ITS = 50, SITE_ITS_MIN = 10, XS_ITS = 30
real,    parameter :: SITE_CONV = 0.005, HALVES_KSIG = 3.
integer, parameter :: SITE_GRID(12) = [500, 1000, 1500, 2500, 4000, 6000, 10000, 15000, 25000, 40000, 60000, 90000]
real,    parameter :: SITE_PFLOOR = 1.e-3
real,    parameter :: EPS_REG = 0.05, EM_LIM = 20., JUNK_SSNR = 0.35, JUNK_PROJ = 0.75, DOCK_MAX = 3.
real,    parameter :: JUNK_LIM = 12., JUNK_GAP = 2., JUNK_MINW = 0.005, MAD_SD = 1.4826
integer, parameter :: JUNK_ITS = 500, NQ_CON = 81, LOC_BOX = 88

type :: ivec
    integer, allocatable :: v(:)
end type ivec

type :: shuffle_null
    integer, allocatable :: sjidx(:), soff(:)
    real,    allocatable :: nul2(:)
    integer :: ng = 0, neq = 0, nsh = 0, ndone = 0, nstep = 0, neqref = 0
    real    :: zc = 0.
    logical :: l_z = .false.
end type shuffle_null

type :: state_mergeup
    private
    class(parameters), pointer :: p => null()
    class(builder),    pointer :: b => null()
    type(fplane_type),   allocatable :: fpls(:), pouts(:), bpls(:)
    type(reconstructor), allocatable :: recs(:)
    type(exp_samples),   allocatable :: smps(:)
    type(image),         allocatable :: wimgs(:)
    type(ori),           allocatable :: os(:)
    type(sym),           allocatable :: syms(:)
    integer, allocatable :: rows(:), eo(:), vbin(:), emhk(:,:), emsh(:), lslot(:)
    real,    allocatable :: eul(:,:), dfm(:), qcorr(:), emw(:), f2s(:,:), cw(:), msk3(:,:,:), fw(:,:,:), msk0(:,:,:), lmsk(:,:,:,:)
    real(dp), allocatable :: fqh(:), fqk(:)
    integer, allocatable :: fsh(:,:,:)
    integer :: n = 0, nslots = 1, box = 0, nh = 0, nyq = 0, klo = 0, khi = 0, kb = 0, kt = 0, nem = 0, nshem = 0
    logical :: l_exists = .false.
    ! merge-up rule: a pair merges only when shown equivalent at margin
    ! mdelta, never across cnl
    logical, allocatable :: cnl(:,:)
    logical :: l_equiv = .false.
    real    :: mdelta = 0.
    integer :: mtests = 1
    integer :: audit_calls = 0
    integer, allocatable :: jhk(:,:), jsh(:)
    real,    allocatable :: jw(:), f2j(:,:)
    real,    allocatable :: nv(:), con(:), cam(:), cse(:), cpj(:)
    real(dp), allocatable :: cb(:,:), cs(:,:)
    real(dp) :: cfit(3) = 0._dp
    integer :: nj = 0, kj = 0, nfit = 0
    integer, allocatable :: fbhk(:,:), fbsh(:)
    real,    allocatable :: fbw(:), f2b(:,:), cwb(:)
    integer :: nfb = 0, nshfb = 0
    logical :: l_con = .false., l_amp = .false.
  contains
    procedure :: new
    procedure :: focus
    procedure :: site_focus
    procedure :: site_split_em
    procedure, private :: site_refs
    procedure :: leaf_focus
    procedure :: run
    procedure :: test_once
    procedure :: quality_flag
    procedure :: pose_lines
    procedure :: write_weights
    procedure :: site_em
    procedure :: fit_contrast
    procedure :: junk_screen
    procedure :: get_rows
    procedure :: holdout_band
    procedure :: group_halves
    procedure :: kill
    procedure :: audit_em_wrapper
    procedure, private :: consensus
    procedure, private :: screen
    procedure, private :: decide_once
    procedure, private :: null_upto_df
    procedure, private :: pair_stat_df
    procedure, private :: noise_power
    procedure, private :: null_resolve
    procedure, private :: cross_em
    procedure, private :: em_maps
    procedure, private :: write_em_audit
    procedure, private :: estep
    procedure, private :: maps
    procedure, private :: maps_weighted_shared
    procedure, private :: accumulate
    procedure, private :: accumulate_jobs
    procedure, private :: finish
    procedure, private :: dock
    procedure, private :: strata
    procedure, private :: shell_fsc
    procedure, private :: shell_corr
    procedure, private :: shell_dot
    procedure, private :: sym_rots
    procedure, private :: band_density
    procedure, private :: band_wiener
    procedure, private :: band_terms
    procedure          :: half_rel
    procedure          :: half_corr
    procedure          :: dev_gram
end type state_mergeup

contains

    subroutine new( self, params, build, em_lim, contrast, amplitude )
        use simple_matcher_3Drec,   only: init_rec, prep_imgs4rec, cleanup_rec_buffers
        use simple_matcher_ptcl_io, only: prepimgbatch, ptcl_read_session
        use simple_imgarr_utils,    only: alloc_imgarr, dealloc_imgarr
        use simple_rnd,             only: seed_rnd_fixed
        class(state_mergeup),      intent(inout) :: self
        class(parameters), target, intent(inout) :: params
        class(builder),    target, intent(inout) :: build
        real,                      intent(in)    :: em_lim
        logical, optional,         intent(in)    :: contrast, amplitude
        type(fplane_type), allocatable :: fdummy(:)
        type(image),       allocatable :: crop_imgs(:)
        real,              allocatable :: ones(:,:,:), sdn(:)
        type(ptcl_read_session) :: reads
        type(image) :: mskimg
        real    :: dirs(3,NVIEW), v(3), t, z, rz, g, e1, e2
        integer :: i, j, ip, h, k, l, r2, sh, b1, b2, nrows, it, ie
        call self%kill
        self%p => params
        self%b => build
        if( present(contrast)  ) self%l_con = contrast
        if( present(amplitude) ) self%l_amp = amplitude
        if( self%l_con .or. self%l_amp ) write(logfhandle,'(A,L1,A,L1)') '>>> STATE_MERGEUP contrast-aware: maps ', self%l_con, &
            &', E-step amplitude ', self%l_amp
        !$ call omp_set_max_active_levels(1)
        self%nslots = max(2, nthr_glob)
        !$ self%nslots = max(self%nslots, omp_get_max_threads())
        call seed_rnd_fixed(RND_SEED)
        nrows = build%spproj_field%get_noris()
        self%rows = pack([(i, i=1,nrows)], [(build%spproj_field%get_state(i) > 0 .and. &
            &(build%spproj_field%get_eo(i) == 0 .or. build%spproj_field%get_eo(i) == 1), i=1,nrows)])
        self%n = size(self%rows)
        if( self%n < 2 ) THROW_HARD('no particles with a state and an even/odd half')
        self%box = params%box_crop
        self%nh  = self%box / 2 + 1
        self%nyq = self%box / 2
        self%klo = max(nint(real(self%box) * params%smpd_crop / params%hp), 3)
        self%khi = self%nyq - 1
        self%kb  = min(self%nyq - 1, nint(real(self%box) * params%smpd_crop / em_lim))
        allocate(self%eo(self%n), self%vbin(self%n), self%eul(3,self%n), self%dfm(self%n), self%qcorr(self%n))
        do i = 1, NVIEW
            t  = real(i-1) + 0.5
            z  = 1. - 2. * t / real(NVIEW)
            rz = sqrt(max(1. - z*z, 0.))
            g  = PI * (3. - sqrt(5.)) * t
            dirs(:,i) = [rz*cos(g), rz*sin(g), z]
        end do
        do i = 1, self%n
            ip = self%rows(i)
            self%eo(i)    = build%spproj_field%get_eo(ip)
            self%eul(:,i) = build%spproj_field%get_euler(ip)
            self%dfm(i)   = 0.5 * (build%spproj_field%get_dfx(ip) + build%spproj_field%get_dfy(ip))
            self%qcorr(i) = build%spproj_field%get(ip, 'corr')
            e1 = self%eul(1,i) * PI / 180.
            e2 = self%eul(2,i) * PI / 180.
            v  = [sin(e2)*cos(e1), sin(e2)*sin(e1), cos(e2)]
            self%vbin(i) = maxloc(matmul(v, dirs), dim=1)
        end do
        write(logfhandle,'(A)') '>>> STATE_MERGEUP loading '//int2str(self%n)//' particles at box_crop '//&
            &int2str(self%box)//' (Nyquist '//fmtr(real(2.*params%smpd_crop,dp),1)//' A)'
        allocate(self%fpls(self%n))
        call init_rec(params, build, MAXIMGBATCHSZ, fdummy, cropped=params%box_crop < params%box)
        if( params%box_crop < params%box ) call alloc_imgarr(nthr_glob, [params%box_crop,params%box_crop,1], &
            &params%smpd_crop, crop_imgs)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        call reads%new(params)
        allocate(sdn(self%n), source=0.)
        do b1 = 1, self%n, MAXIMGBATCHSZ
            b2 = min(self%n, b1 + MAXIMGBATCHSZ - 1)
            call reads%read(params, build, self%n, self%rows, [b1,b2])
            !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
            do i = b1, b2
                call build%imgbatch(i-b1+1)%norm_noise(build%lmsk, sdn(i))
            end do
            !$omp end parallel do
            call prep_imgs4rec(params, build, b2-b1+1, build%imgbatch(:b2-b1+1), self%rows(b1:b2), self%fpls(b1:b2), &
                &crop_imgs=crop_imgs, kfromto=[1,self%nyq], observation_model=.true., whiten=.false., contrast=self%l_con)
        end do
        call reads%kill
        call cleanup_rec_buffers(build, fdummy)
        if( allocated(crop_imgs) ) call dealloc_imgarr(crop_imgs)
        do ie = 1, 2
            self%nem = 0
            do k = lbound(self%fpls(1)%cmplx_plane,2), 0
                do h = lbound(self%fpls(1)%cmplx_plane,1), ubound(self%fpls(1)%cmplx_plane,1)
                    r2 = h*h + k*k
                    if( r2 < 3 .or. r2 > self%kb*(self%kb+1) ) cycle
                    self%nem = self%nem + 1
                    if( ie == 1 ) cycle
                    self%emhk(:,self%nem) = [h, k]
                    self%emsh(self%nem)   = nint(sqrt(real(r2))) - 1
                    self%emw(self%nem)    = merge(1., 2., k == 0)
                end do
            end do
            if( ie == 1 ) allocate(self%emhk(2,self%nem), self%emsh(self%nem), self%emw(self%nem))
        end do
        self%nshem = maxval(self%emsh)
        allocate(self%f2s(self%nshem,self%n), self%cw(self%nshem), source=0.)
        do j = 1, self%nem
            self%cw(self%emsh(j)) = self%cw(self%emsh(j)) + self%emw(j)
        end do
        ! junk screen band, shells 2..kj (JUNK_LIM): plane samples and per-shell |F|^2 taken before the CTF is applied
        self%kj = min(self%nyq - 1, nint(real(self%box) * params%smpd_crop / JUNK_LIM))
        do ie = 1, 2
            self%nj = 0
            do k = lbound(self%fpls(1)%cmplx_plane,2), 0
                do h = lbound(self%fpls(1)%cmplx_plane,1), ubound(self%fpls(1)%cmplx_plane,1)
                    r2 = h*h + k*k
                    if( r2 < 3 .or. r2 > self%kj*(self%kj+1) ) cycle
                    self%nj = self%nj + 1
                    if( ie == 1 ) cycle
                    self%jhk(:,self%nj) = [h, k]
                    self%jsh(self%nj)   = nint(sqrt(real(r2)))
                    self%jw(self%nj)    = merge(1., 2., k == 0)
                end do
            end do
            if( ie == 1 ) allocate(self%jhk(2,self%nj), self%jsh(self%nj), self%jw(self%nj))
        end do
        allocate(self%f2j(self%kj,self%n), source=0.)
        ! the full-band sample table of the choose rule's membership pass
        do ie = 1, 2
            self%nfb = 0
            do k = lbound(self%fpls(1)%cmplx_plane,2), 0
                do h = lbound(self%fpls(1)%cmplx_plane,1), ubound(self%fpls(1)%cmplx_plane,1)
                    r2 = h*h + k*k
                    if( r2 < 3 .or. r2 > self%khi*(self%khi+1) ) cycle
                    self%nfb = self%nfb + 1
                    if( ie == 1 ) cycle
                    self%fbhk(:,self%nfb) = [h, k]
                    self%fbsh(self%nfb)   = nint(sqrt(real(r2))) - 1
                    self%fbw(self%nfb)    = merge(1., 2., k == 0)
                end do
            end do
            if( ie == 1 ) allocate(self%fbhk(2,self%nfb), self%fbsh(self%nfb), self%fbw(self%nfb))
        end do
        self%nshfb = maxval(self%fbsh)
        allocate(self%f2b(self%nshfb,self%n), self%cwb(self%nshfb), source=0.)
        do j = 1, self%nfb
            self%cwb(self%fbsh(j)) = self%cwb(self%fbsh(j)) + self%fbw(j)
        end do
        !$omp parallel do default(shared) private(i,j,h,k,sh) schedule(static) proc_bind(close)
        do i = 1, self%n
            if( sdn(i) > 0. ) self%fpls(i)%cmplx_plane = sdn(i) * self%fpls(i)%cmplx_plane
            do j = 1, self%nem
                h  = self%emhk(1,j)
                k  = self%emhk(2,j)
                sh = self%emsh(j)
                self%f2s(sh,i) = self%f2s(sh,i) + self%emw(j) * abs(self%fpls(i)%cmplx_plane(h,k))**2
            end do
            do j = 1, self%nj
                self%f2j(self%jsh(j),i) = self%f2j(self%jsh(j),i) + self%jw(j) * &
                    &abs(self%fpls(i)%cmplx_plane(self%jhk(1,j),self%jhk(2,j)))**2
            end do
            do j = 1, self%nfb
                self%f2b(self%fbsh(j),i) = self%f2b(self%fbsh(j),i) + self%fbw(j) * &
                    &abs(self%fpls(i)%cmplx_plane(self%fbhk(1,j),self%fbhk(2,j)))**2
            end do
            self%fpls(i)%cmplx_plane = real(self%fpls(i)%transfer_plane) * self%fpls(i)%cmplx_plane
            deallocate(self%fpls(i)%transfer_plane)
        end do
        !$omp end parallel do
        self%nv = merge(sdn * sdn, 1., sdn > 0.)
        allocate(self%fsh(self%nh,self%box,self%box), self%fw(self%nh,self%box,self%box))
        do l = 1, self%box
            do k = 1, self%box
                do h = 1, self%nh
                    self%fsh(h,k,l) = min(nint(sqrt(real((h-1)**2 + (k-1-merge(self%box,0,k>self%nyq))**2 + &
                        &(l-1-merge(self%box,0,l>self%nyq))**2))), self%nyq)
                    self%fw(h,k,l)  = merge(1., 2., h == 1 .or. h == self%nh)
                end do
            end do
        end do
        allocate(self%fqh(self%nh), self%fqk(self%box))
        do h = 1, self%nh
            self%fqh(h) = real(h-1,dp) / real(params%box,dp)
        end do
        do k = 1, self%box
            self%fqk(k) = real(k-1-merge(self%box,0,k>self%nyq),dp) / real(params%box,dp)
        end do
        self%kt = min(self%kb + RB_MARGIN, self%khi - 2)
        allocate(self%recs(self%nslots), self%os(self%nslots), self%syms(self%nslots), self%smps(self%nslots), &
            &self%pouts(self%nslots), self%wimgs(self%nslots), self%bpls(self%nslots))
        do it = 1, self%nslots
            call self%recs(it)%new_accumulator(params, build%spproj, expand=.true., wthreads=.false.)
            call self%os(it)%new(is_ptcl=.true.)
            call self%syms(it)%new(trim(params%pgrp))
            call self%wimgs(it)%new([self%box,self%box,self%box], params%smpd_crop, wthreads=.false.)
            self%bpls(it) = self%fpls(1)
            self%bpls(it)%cmplx_plane = cmplx(0.,0.)
            self%bpls(it)%ctfsq_plane = 0.
        end do
        allocate(ones(self%box,self%box,self%box), source=1.)
        call mskimg%new([self%box,self%box,self%box], params%smpd_crop, wthreads=.false.)
        call mskimg%set_rmat(ones, .false.)
        call mskimg%mask3D_soft(params%msk_crop, backgr=0.)
        self%msk3 = mskimg%get_rmat()
        call mskimg%kill
        write(logfhandle,'(A)') '>>> STATE_MERGEUP test band shells '//int2str(self%klo)//'..'//int2str(self%khi)//&
            &' ('//fmtr(real(real(self%box)*params%smpd_crop/real(self%klo),dp),1)//'..'//&
            &fmtr(real(real(self%box)*params%smpd_crop/real(self%khi),dp),1)//' A); cross-half EM shells 2..'//&
            &int2str(self%kb)//' ('//int2str(self%nem)//' plane samples)'
        self%l_exists = .true.
    end subroutine new

    subroutine focus( self, fname )
        class(state_mergeup), intent(inout) :: self
        class(string),        intent(in)    :: fname
        type(image) :: vol
        real, allocatable :: fm(:,:,:)
        call vol%read_and_crop(fname, self%p%smpd, self%box, self%p%smpd_crop)
        fm = min(max(vol%get_rmat(), 0.), 1.)
        if( any(shape(fm) /= shape(self%msk3)) ) THROW_HARD('focus mask box does not match the particles; focus')
        self%msk3 = self%msk3 * fm
        self%kt   = min(self%kb + RB_MARGIN_FOCUS, self%khi - 2)
        write(logfhandle,'(A)') '>>> STATE_MERGEUP focus mask '//fname%to_char()//': '//&
            &fmtr(100._dp*real(count(fm > 0.5),dp)/real(size(fm),dp),1)//' % of the crop grid above 0.5'
        call vol%kill
    end subroutine focus

    !> the site split's focus: the node mask times one site mask, reset for every site
    subroutine site_focus( self, fname )
        class(state_mergeup), intent(inout) :: self
        class(string),        intent(in)    :: fname
        type(image) :: vol
        real, allocatable :: fm(:,:,:)
        if( .not. allocated(self%msk0) ) self%msk0 = self%msk3
        call vol%read_and_crop(fname, self%p%smpd, self%box, self%p%smpd_crop)
        fm = min(max(vol%get_rmat(), 0.), 1.)
        if( any(shape(fm) /= shape(self%msk0)) ) THROW_HARD('site mask box does not match the particles; site_focus')
        self%msk3 = self%msk0 * fm
        self%kt   = min(self%kb + RB_MARGIN_FOCUS, self%khi - 2)
        call vol%kill
    end subroutine site_focus

    !> the site split's cross-validated classification: per-half random starts, each half scored against the other
    !! half's site-masked references for niter iterations; labels over the node's rows (0 = in no half) and the
    !! groups' unmasked half maps from the final labels (k + 1: the groups, then every particle)
    subroutine site_split_em( self, k, niter, seed, labels, ma, mb )
        use simple_rnd, only: seed_rnd_fixed, ran3
        class(state_mergeup), intent(inout) :: self
        integer,              intent(in)    :: k, niter, seed
        integer, allocatable, intent(out)   :: labels(:)
        complex, allocatable, intent(out)   :: ma(:,:,:,:), mb(:,:,:,:)
        complex, allocatable :: cons(:,:,:,:)
        integer, allocatable :: ia(:), ib(:), la(:), lb(:), starta(:), startb(:)
        real,    allocatable :: eps(:), wa(:,:), wb(:,:)
        integer :: neq, i, ii
        ia = pack([(ii, ii=1,self%n)], self%eo == 0)
        ib = pack([(ii, ii=1,self%n)], self%eo == 1)
        allocate(la(size(ia)), lb(size(ib)), source=1)
        call self%consensus(ia, ib, la, lb, 1, neq, eps, cons)
        call seed_rnd_fixed(seed)
        allocate(starta(size(ia)), startb(size(ib)))
        do i = 1, size(ia)
            starta(i) = min(k, 1 + int(ran3() * real(k)))
        end do
        do i = 1, size(ib)
            startb(i) = min(k, 1 + int(ran3() * real(k)))
        end do
        call self%cross_em(ia, ib, starta, startb, k, eps, neq, la, lb, l_full=.true., niter=niter, l_site=.true.)
        allocate(wa(size(ia),k), wb(size(ib),k), source=0.)
        do i = 1, size(ia)
            if( la(i) > 0 ) wa(i,la(i)) = 1.
        end do
        do i = 1, size(ib)
            if( lb(i) > 0 ) wb(i,lb(i)) = 1.
        end do
        allocate(ma(self%nh,self%box,self%box,k+1), mb(self%nh,self%box,self%box,k+1))
        call self%em_maps(ia, wa, eps, neq, .true., ma)
        call self%em_maps(ib, wb, eps, neq, .true., mb)
        allocate(labels(self%n), source=0)
        labels(ia) = la
        labels(ib) = lb
    end subroutine site_split_em

    !> a half's site references: its consensus plus the site-masked deviation of each group from the groups' weighted mean
    subroutine site_refs( self, m, w )
        class(state_mergeup), intent(inout) :: self
        complex,              intent(inout) :: m(:,:,:,:)
        real,                 intent(in)    :: w(:,:)
        real, pointer     :: rmat(:,:,:)
        real, allocatable :: r(:,:,:,:), mu(:,:,:)
        real    :: f(size(w,2))
        integer :: g, k
        k = size(w,2)
        f = [(sum(w(:,g)), g=1,k)]
        f = f / max(sum(f), 1.)
        allocate(r(self%box,self%box,self%box,k+1))
        do g = 1, k + 1
            call self%wimgs(1)%set_ft(.true.)
            call self%wimgs(1)%set_cmat(m(:,:,:,g))
            call self%wimgs(1)%ifft
            call self%wimgs(1)%get_rmat_ptr(rmat)
            r(:,:,:,g) = rmat(:self%box,:self%box,:self%box)
        end do
        mu = f(1) * r(:,:,:,1)
        do g = 2, k
            mu = mu + f(g) * r(:,:,:,g)
        end do
        do g = 1, k
            call self%wimgs(1)%set_ft(.false.)
            call self%wimgs(1)%get_rmat_ptr(rmat)
            rmat(:self%box,:self%box,:self%box) = r(:,:,:,k+1) + self%msk3 * (r(:,:,:,g) - mu)
            call self%wimgs(1)%fft
            m(:,:,:,g) = self%wimgs(1)%get_cmat()
        end do
    end subroutine site_refs

    subroutine leaf_focus( self, fname )
        class(state_mergeup), intent(inout) :: self
        class(string),        intent(in)    :: fname
        character(len=4096), allocatable :: files(:)
        integer,             allocatable :: labs(:), sl(:)
        real,                allocatable :: fm(:,:,:)
        character(len=4096) :: line
        type(image) :: vol
        integer :: funit, ios, nl, il, jl, j, ns
        nl = nlines(fname)
        if( nl == 0 ) return
        allocate(files(nl), labs(nl), sl(nl))
        call fopen(funit, fname, status='old', action='read', iostat=ios)
        call fileiochk('state_mergeup: cannot read '//fname%to_char(), ios)
        do il = 1, nl
            read(funit,'(A)') line
            line = adjustl(line)
            j    = index(line, ' ')
            read(line(1:j-1),*) labs(il)
            files(il) = adjustl(line(j+1:))
        end do
        call fclose(funit)
        ns = 0
        do il = 1, nl
            sl(il) = 0
            do jl = 1, il - 1
                if( files(jl) == files(il) ) sl(il) = sl(jl)
            end do
            if( sl(il) > 0 ) cycle
            ns     = ns + 1
            sl(il) = ns
        end do
        allocate(self%lmsk(self%box,self%box,self%box,ns), self%lslot(0:maxval(labs)))
        self%lslot = 0
        do il = 1, nl
            self%lslot(labs(il)) = sl(il)
            if( any(sl(1:il-1) == sl(il)) ) cycle
            call vol%read_and_crop(string(trim(files(il))), self%p%smpd, self%box, self%p%smpd_crop)
            fm = min(max(vol%get_rmat(), 0.), 1.)
            if( any(shape(fm) /= shape(self%msk3)) ) THROW_HARD('focus mask box does not match the particles; leaf_focus')
            self%lmsk(:,:,:,sl(il)) = fm
            call vol%kill
        end do
        self%msk0 = self%msk3
        self%kt   = min(self%kb + RB_MARGIN_FOCUS, self%khi - 2)
        write(logfhandle,'(A)') '>>> STATE_MERGEUP per-state focus masks '//fname%to_char()//': '//int2str(nl)//&
            &' states, '//int2str(ns)//' masks; each group is tested inside the union of its states'' masks'
    end subroutine leaf_focus

    function get_rows( self ) result( rows )
        class(state_mergeup), intent(in) :: self
        integer, allocatable :: rows(:)
        rows = self%rows
    end function get_rows

    function holdout_band( self ) result( b )
        class(state_mergeup), intent(in) :: self
        integer :: b(2)
        b = [self%kb + RB_MARGIN, self%khi]
    end function holdout_band

    real function mask_neff( mfile, smpd, box_crop, smpd_crop, msk_crop, shlo, shhi )
        class(string), intent(in) :: mfile
        real,          intent(in) :: smpd, smpd_crop, msk_crop
        integer,       intent(in) :: box_crop, shlo, shhi
        type(image) :: vol, sph
        real, allocatable :: m(:,:,:), ones(:,:,:)
        real(dp) :: s2, s4
        integer  :: h, k, l, hh, kk, ll, sh, nyq, nband
        call vol%read_and_crop(mfile, smpd, box_crop, smpd_crop)
        m = min(max(vol%get_rmat(), 0.), 1.)
        call vol%kill
        allocate(ones(box_crop,box_crop,box_crop), source=1.)
        call sph%new([box_crop,box_crop,box_crop], smpd_crop, wthreads=.false.)
        call sph%set_rmat(ones, .false.)
        call sph%mask3D_soft(msk_crop, backgr=0.)
        m = m * sph%get_rmat()
        call sph%kill
        s2    = sum(real(m,dp)**2)
        s4    = sum(real(m,dp)**4)
        nyq   = box_crop / 2
        nband = 0
        do l = 1, box_crop
            ll = l - 1 - merge(box_crop, 0, l > nyq)
            do k = 1, box_crop
                kk = k - 1 - merge(box_crop, 0, k > nyq)
                do h = 1, box_crop
                    hh = h - 1 - merge(box_crop, 0, h > nyq)
                    sh = min(nint(sqrt(real(hh*hh + kk*kk + ll*ll))), nyq)
                    if( sh >= shlo .and. sh <= shhi ) nband = nband + 1
                end do
            end do
        end do
        mask_neff = real(0.5_dp * real(nband,dp) * s2 * s2 / max(s4, 1.e-30_dp) / real(box_crop,dp)**3)
    end function mask_neff

    subroutine run( self, tree, labels, units )
        class(state_mergeup), intent(inout) :: self
        class(state_tree),    intent(in)    :: tree
        integer,              intent(in)    :: labels(:)
        integer, allocatable, intent(out)   :: units(:)
        type(ivec), allocatable :: ul(:), nl(:), nodes(:)
        integer,    allocatable :: grp(:), names(:), l2u(:), uptcl(:), cand(:), child(:), newu(:)
        integer,    allocatable :: to(:), merged(:), junk(:), tst(:)
        logical,    allocatable :: isj(:), newj(:)
        integer,    allocatable :: dpairs(:,:), ugrp(:)
        integer :: g, lv, gi, il, u, nid, j, kstar, c, nmerge, ntest, nbefore
        units = max(labels, 0)
        call merge_rule_setup
        allocate(ul(maxval(units)))
        do u = 1, size(ul)
            if( count(units == u) > 0 )then
                ul(u)%v = [u]
            else
                allocate(ul(u)%v(0))
            endif
        end do
        allocate(isj(size(ul)), source=.false.)
        call tree%get_names(names)
        g  = count([(size(ul(u)%v) > 0, u=1,size(ul))])
        lv = 0
        write(logfhandle,'(A)') '>>> STATE_MERGEUP '//int2str(g)//' leaf states, '//int2str(tree%get_nleaves())//&
            &' in the tree'
        do while( g > 1 )
            g  = (g + 1) / 2
            lv = lv + 1
            nbefore = count([(size(ul(u)%v) > 0, u=1,size(ul))])
            call tree%cut(g, grp)
            allocate(l2u(0:max(maxval(names), 1)), source=0)
            do u = 1, size(ul)
                do j = 1, size(ul(u)%v)
                    if( ul(u)%v(j) <= ubound(l2u,1) ) l2u(ul(u)%v(j)) = u
                end do
            end do
            uptcl = units(self%rows)
            allocate(nodes(g))
            ntest = 0
            do gi = 1, g
                allocate(nodes(gi)%v(0))
                do il = 1, size(names)
                    if( grp(il) /= gi ) cycle
                    u = l2u(names(il))
                    if( u == 0 ) cycle
                    if( any(nodes(gi)%v == u) ) cycle
                    if( count(uptcl == u .and. self%eo == 0) == 0 .or. count(uptcl == u .and. self%eo == 1) == 0 ) cycle
                    nodes(gi)%v = [nodes(gi)%v, u]
                end do
                call isort(nodes(gi)%v)
                if( count(.not. isj(nodes(gi)%v)) >= 2 ) ntest = ntest + 1
            end do
            write(logfhandle,'(A)') '>>> level '//int2str(lv)//': cut into '//int2str(g)//' groups: '//int2str(ntest)//&
                &' nodes to test ('//int2str(nbefore)//' units, '//int2str(count(isj))//' junk passed through)'
            allocate(newu(size(units)), source=0)
            allocate(nl(size(ul)), newj(size(ul)))
            newj   = .false.
            nid    = 0
            nmerge = 0
            do gi = 1, g
                if( size(nodes(gi)%v) == 0 ) cycle
                do j = 1, size(nodes(gi)%v)
                    if( isj(nodes(gi)%v(j)) ) call pass(nodes(gi)%v(j), .true.)
                end do
                tst = pack(nodes(gi)%v, .not. isj(nodes(gi)%v))
                if( size(tst) == 0 ) cycle
                if( size(tst) < 2 )then
                    call pass(tst(1), .false.)
                    cycle
                endif
                allocate(cand(self%n), source=0)
                do j = 1, size(tst)
                    where( uptcl == tst(j) ) cand = j
                end do
                write(logfhandle,'(A)') '>>> level '//int2str(lv)//' node '//int2str(gi)//': units '//joined(tst)
                if( allocated(self%lslot) ) call group_mask(tst)
                call link_matrix(tst)
                call self%test_once(cand, size(tst), child, kstar, junk, self%mtests, ugrp=ugrp)
                deallocate(self%cnl)
                if( allocated(self%lslot) ) self%msk3 = self%msk0
                do j = 1, size(junk)
                    call pass(tst(junk(j)), .true.)
                end do
                ! units follow decide_once's groups, so a unit the EM drained (group 0) passes alone
                allocate(to(size(tst)), source=0)
                do j = 1, size(tst)
                    if( any(junk == j) )then
                        to(j) = -1
                        cycle
                    endif
                    to(j) = ugrp(j)
                    if( to(j) == 0 ) call pass(tst(j), .false.)
                end do
                do c = 1, max(1,kstar)
                    if( .not. any(to == c) ) cycle
                    nid = nid + 1
                    allocate(merged(0))
                    do j = 1, size(tst)
                        if( to(j) /= c ) cycle
                        u = tst(j)
                        merged = [merged, ul(u)%v]
                        where( units == u ) newu = nid
                    end do
                    call isort(merged)
                    nl(nid)%v = merged
                    nmerge = nmerge + count(to == c) - 1
                    deallocate(merged)
                end do
                write(logfhandle,'(A)') '>>>   level '//int2str(lv)//' node '//int2str(gi)//': '//&
                    &int2str(size(tst))//' units -> K* '//int2str(kstar)//', '//int2str(size(junk))//' junk'
                deallocate(cand, to, child, junk)
            end do
            units = newu
            ul    = nl(1:nid)
            isj   = newj(1:nid)
            write(logfhandle,'(A)') '>>> level '//int2str(lv)//': '//int2str(nbefore)//' units -> '//int2str(nid)//&
                &' units ('//int2str(nmerge)//' merges, '//int2str(count(isj))//' junk)'
            call write_units('level'//int2str(lv))
            deallocate(newu, nl, newj, nodes, l2u, grp)
        end do
        call write_units('merged')
        write(logfhandle,'(A)') '>>> STATE_MERGEUP DONE: '//int2str(size(ul))//' states, '//int2str(count(isj))//' junk'

    contains

        subroutine merge_rule_setup
            character(len=4096) :: ev
            character(len=:), allocatable :: tdir
            integer, allocatable :: npairs(:,:)
            real    :: ndelta
            integer :: ios, nrs
            self%l_equiv = .true.
            call get_environment_variable('SIMPLE_MERGE_TREEDIR', ev, status=ios)
            if( ios == 0 .and. len_trim(ev) > 0 )then
                tdir = trim(ev)
            else
                tdir = parent_dir(parent_dir(self%p%infile%to_char()))
            endif
            call equiv_tree_inputs(tdir, ndelta, npairs, nrs, self%mtests)
            write(logfhandle,'(A)') '>>> STATE_MERGEUP merge rule inputs from '//tdir//': '//int2str(nrs)//&
                &' node splits kept K* >= 2, smallest effect '//fmtr(real(ndelta,dp),6)//'; '//int2str(size(npairs,2))//&
                &' sibling leaf pairs shown distinct:'//pair_list(npairs)//'; '//int2str(self%mtests)//' split tests in the tree'
            self%mdelta = ndelta
            call move_alloc(npairs, dpairs)
            write(logfhandle,'(A)') '>>> STATE_MERGEUP merge rule: a pair merges only when shown equivalent at margin '//&
                &fmtr(real(self%mdelta,dp),6)//'; '//int2str(size(dpairs,2))//' leaf pairs shown distinct are never merged'
        end subroutine merge_rule_setup

        subroutine link_matrix( uu )
            integer, intent(in) :: uu(:)
            integer :: a, b, k
            allocate(self%cnl(size(uu),size(uu)), source=.false.)
            do k = 1, size(dpairs,2)
                do a = 1, size(uu)
                    if( .not. any(ul(uu(a))%v == dpairs(1,k)) ) cycle
                    do b = 1, size(uu)
                        if( a == b .or. .not. any(ul(uu(b))%v == dpairs(2,k)) ) cycle
                        self%cnl(a,b) = .true.
                        self%cnl(b,a) = .true.
                    end do
                end do
            end do
            write(logfhandle,'(A)') '>>>   merge rule: '//int2str(count(self%cnl) / 2)//' cannot-link unit pairs in this node'
        end subroutine link_matrix

        subroutine pass( uu, lj )
            integer, intent(in) :: uu
            logical, intent(in) :: lj
            nid = nid + 1
            nl(nid)%v = ul(uu)%v
            newj(nid) = lj
            where( units == uu ) newu = nid
        end subroutine pass

        subroutine group_mask( uu )
            integer, intent(in) :: uu(:)
            real, allocatable :: fm(:,:,:)
            integer :: a, b, s
            allocate(fm(self%box,self%box,self%box), source=0.)
            do a = 1, size(uu)
                do b = 1, size(ul(uu(a))%v)
                    s = 0
                    if( ul(uu(a))%v(b) <= ubound(self%lslot,1) ) s = self%lslot(ul(uu(a))%v(b))
                    if( s == 0 )then
                        fm = 1.
                    else
                        fm = max(fm, self%lmsk(:,:,:,s))
                    endif
                end do
            end do
            self%msk3 = self%msk0 * fm
            write(logfhandle,'(A)') '>>>   focus: union of the units'' masks, '//&
                &fmtr(100._dp*real(count(fm > 0.5),dp)/real(size(fm),dp),1)//' % of the crop grid above 0.5'
        end subroutine group_mask

        subroutine write_units( fbody )
            character(len=*), intent(in) :: fbody
            type(string) :: fname
            integer :: funit, ios, ii
            fname = fbody//'_units.txt'
            call fopen(funit, fname, status='replace', action='write', iostat=ios)
            call fileiochk('state_mergeup: cannot write '//fname%to_char(), ios)
            do ii = 1, size(units)
                write(funit,'(I0)') units(ii)
            end do
            call fclose(funit)
            fname = fbody//'_leaves.txt'
            call fopen(funit, fname, status='replace', action='write', iostat=ios)
            call fileiochk('state_mergeup: cannot write '//fname%to_char(), ios)
            write(funit,'(A)') '# unit leaves'
            do ii = 1, size(ul)
                write(funit,'(A)') int2str(ii)//' '//joined(ul(ii)%v, ',')
            end do
            call fclose(funit)
            fname = fbody//'_junk.txt'
            call fopen(funit, fname, status='replace', action='write', iostat=ios)
            call fileiochk('state_mergeup: cannot write '//fname%to_char(), ios)
            do ii = 1, size(isj)
                if( isj(ii) ) write(funit,'(I0)') ii
            end do
            call fclose(funit)
        end subroutine write_units

    end subroutine run

    ! the single-pass split test (audit 6.7): the SSNR junk screen, one cross-half EM from the stage-B labels, one equal N (the
    ! smallest EM group in either half), maps of every group built once at it, every pair statistic read off them in the hold-out
    ! band, one stratified shuffle null at the same N and one z_crit; the weakest kept pair, the null's mean and sd and N as numbers
    subroutine test_once( self, cand, k0, child, kstar, junk, ntests, ugrp, stats )
        class(state_mergeup), intent(inout) :: self
        integer,              intent(in)    :: cand(:), k0, ntests
        integer, allocatable, intent(out)   :: child(:), junk(:)
        integer,              intent(out)   :: kstar
        integer, allocatable, optional, intent(out) :: ugrp(:)
        real,    optional,    intent(out)   :: stats(6)
        complex, allocatable :: cons(:,:,:,:)
        integer, allocatable :: ia(:), ib(:), la(:), lb(:), sta(:), stb(:), pa(:), pb(:), gcom(:), starta(:), startb(:)
        integer, allocatable :: plist(:), gof(:), exa(:), exb(:), sizes(:), order(:), rank(:), rkc(:)
        real,    allocatable :: eps(:)
        type(shuffle_null) :: nl
        real    :: mu2, sd2, tmin, tmin0, thr
        integer :: neq, neq2, s, g, m, ii, jj, kem, ng, nsj, nthr, npair
        ia = pack([(ii, ii=1,self%n)], cand > 0 .and. self%eo == 0)
        ib = pack([(ii, ii=1,self%n)], cand > 0 .and. self%eo == 1)
        la = cand(ia)
        lb = cand(ib)
        write(logfhandle,'(A)') 'A: '//int2str(size(ia))//' node particles (labels '//joined(counts(la, k0))//&
            &'); B: '//int2str(size(ib))//' (labels '//joined(counts(lb, k0))//')'
        call self%consensus(ia, ib, la, lb, k0, neq, eps, cons)
        call self%screen(ia, ib, la, lb, k0, neq, eps, cons, junk)
        gcom = pack([(g, g=1,k0)], [(count(la == g) > 0 .and. count(lb == g) > 0, g=1,k0)])
        ng   = size(gcom)
        sta  = self%strata(ia)
        stb  = self%strata(ib)
        nl%ng  = ng
        nl%neq = neq
        nl%nsh = NSHUF_BONF
        nl%l_z = .true.
        nl%zc  = split_zcrit(ntests)
        nthr   = 1
        !$ nthr = omp_get_max_threads()
        nl%nstep = max(NSHUF_STEP, ceiling(real(nthr) / real(max(1, 2*ng))))
        allocate(nl%nul2(0))
        write(logfhandle,'(A)') 'single pass: equal N '//int2str(neq)//' per group; test band '//int2str(self%kt)//'..'//&
            &int2str(self%khi)//'; up to '//int2str(nl%nsh)//' shuffles; stratified null: '//int2str(nuniq(sta))//' / '//&
            &int2str(nuniq(stb))//' nuisance strata (view x defocus); z_crit '//fmtr(real(nl%zc,dp),2)//' (alpha 1/21 over '//&
            &int2str(ntests)//' split tests)'
        kstar = min(ng, 1)
        mu2   = 0.
        sd2   = 0.
        tmin  = 0.
        tmin0 = 0.
        kem   = ng
        allocate(plist(0), gof(0), exa(0:k0), exb(0:k0))
        exa = 0
        exb = 0
        if( ng >= 2 )then
            ! the cross-half EM from the stage-B labels of every group present in both halves
            allocate(starta(size(ia)), startb(size(ib)), source=0)
            do m = 1, ng
                where( la == gcom(m) ) starta = m
                where( lb == gcom(m) ) startb = m
            end do
            pa = la
            pb = lb
            call self%cross_em(ia, ib, starta, startb, kem, eps, neq, la, lb)
            write(logfhandle,'(A)') 'cross-half EM half A: MAP changed vs start '//fmtr(real(count(la /= starta),dp)/&
                &real(size(la),dp),3)//'; half B '//fmtr(real(count(lb /= startb),dp)/real(size(lb),dp),3)
            ! one common n for the pairs and the null: the smallest EM group in either half, at most the stage-B equal N; the
            ! maps themselves use every particle (draw-free), the shuffles are stored whole
            neq2 = neq
            do g = 1, kem
                if( count(la == g) > 0 .and. count(lb == g) > 0 ) neq2 = min(neq2, count(la == g), count(lb == g))
            end do
            nl%neq    = neq2
            nl%neqref = neq
            allocate(nl%sjidx(nl%nsh * (size(ia) + size(ib))), nl%soff(0:2*ng*nl%nsh))
            nl%soff(0) = 0
            nsj = 0
            do s = 1, nl%nsh
                call add_full_jobs(ia, sshuf(pa, sta), (s-1)*2*ng)
                call add_full_jobs(ib, sshuf(pb, stb), (s-1)*2*ng + ng)
            end do
            call self%null_upto_df(nl, NSHUF, eps, cons)
            call null_thr(nl, mu2, sd2, thr)
            write(logfhandle,'(A)') 'draw-free statistic at the common n '//int2str(neq2)//' (stage-B equal N '//int2str(neq)//&
                &'); shuffle null in the test band at that n ('//int2str(size(nl%nul2))//' pair statistics): mean / sd '//&
                &fmtr(real(mu2,dp),4)//' / '//fmtr(real(sd2,dp),4)//', max '//fmtr(real(maxval(nl%nul2),dp),3)//&
                &' -> threshold '//fmtr(real(thr,dp),3)//' from '//int2str(nl%ndone)//' shuffles'
            if( allocated(self%cnl) ) call link_em_groups
            call self%decide_once(ia, ib, la, lb, kem, eps, neq, neq2, cons, nl, plist, gof, exa, exb, kstar, tmin0, tmin)
            call null_thr(nl, mu2, sd2, thr)
            write(logfhandle,'(A)') 'test band '//int2str(self%kt)//'..'//int2str(self%khi)//': null mean / sd '//&
                &fmtr(real(mu2,dp),4)//' / '//fmtr(real(sd2,dp),4)//' from '//int2str(nl%ndone)//' of '//int2str(nl%nsh)//&
                &' shuffles -> threshold '//fmtr(real(thr,dp),3)//'; weakest pair at the start '//&
                &fmtr(real(tmin0,dp),3)//' (z '//fmtr(real((tmin0 - mu2) / max(sd2, 1.e-12),dp),2)//'), weakest kept pair '//&
                &fmtr(real(tmin,dp),3)//' -> '//merge('PASS', 'FAIL', kstar >= 2)//', K* '//int2str(kstar)
        else
            if( ng == 1 )then
                plist = [1]
                gof   = [1]
                where( la == gcom(1) ) la = 1
                where( lb == gcom(1) ) lb = 1
            endif
        endif
        npair = size(plist)
        allocate(sizes(max(1,kstar)), order(max(1,kstar)), rank(max(1,kstar)), source=0)
        do m = 1, npair
            sizes(gof(m)) = sizes(gof(m)) + count(la == plist(m)) + count(lb == plist(m))
        end do
        order = [(ii, ii=1,max(1,kstar))]
        do ii = 2, kstar
            jj = ii
            do while( jj > 1 )
                if( sizes(order(jj-1)) >= sizes(order(jj)) ) exit
                order(jj-1:jj) = [order(jj), order(jj-1)]
                jj = jj - 1
            end do
        end do
        do ii = 1, max(1,kstar)
            rank(order(ii)) = ii
        end do
        allocate(rkc(0:max(1,kstar)))
        rkc = [(ii, ii=0,max(1,kstar))]
        allocate(child(self%n), source=0)
        do ii = 1, size(ia)
            child(ia(ii)) = group_of(la(ii), exa)
        end do
        do ii = 1, size(ib)
            child(ib(ii)) = group_of(lb(ii), exb)
        end do
        if( kstar >= 2 ) call choose_members
        if( present(ugrp) )then
            ! a unit's group is the final one of the EM group it started (decide_once's merges, the choose renumbering), 0 if drained
            allocate(ugrp(k0), source=0)
            do m = 1, ng
                g = findloc(plist, m, dim=1)
                if( g > 0 ) ugrp(gcom(m)) = rkc(rank(gof(g)))
            end do
        endif
        if( present(stats) ) stats = [tmin, tmin0, mu2, sd2, nl%zc, real(nl%neq)]
        write(logfhandle,'(A)') 'K* = '//int2str(kstar)//'; sizes '//joined([(count(child == ii), ii=1,max(1,kstar))])//&
            &'; junk units '//joined(junk)
        deallocate(cons)

    contains

        integer function group_of( lab, ex )
            integer, intent(in) :: lab, ex(0:)
            integer :: mm, lab_here
            group_of = 0
            lab_here = lab
            if( lab_here <= ubound(ex,1) )then
                if( ex(lab_here) > 0 ) lab_here = ex(lab_here)
            endif
            do mm = 1, npair
                if( plist(mm) == lab_here ) group_of = rank(gof(mm))
            end do
        end function group_of

        ! one shuffle's groups of one half, whole, as consecutive jobs j0+1..j0+ng of the null's storage
        subroutine add_full_jobs( idx, lab, j0 )
            integer, intent(in) :: idx(:), lab(:), j0
            integer, allocatable :: sel(:)
            integer :: gg
            do gg = 1, ng
                sel = pack(idx, lab == gcom(gg))
                nl%sjidx(nsj+1:nsj+size(sel)) = sel
                nsj = nsj + size(sel)
                nl%soff(j0+gg) = nsj
            end do
        end subroutine add_full_jobs

        ! membership by the choose rule (the production default): one full-band cross-half pass from the K* groups and one from
        ! the stage-B states (each state to the group holding most of its particles); the states win when their per-particle
        ! cross-half log-likelihood in the test band exceeds the groups' by 3 SE
        subroutine choose_members
            integer,  allocatable :: ja(:), jb(:), a1(:), b1(:), a2(:), b2(:), smap(:), cnt(:,:), sz(:), rk(:)
            real(dp), allocatable :: pa1(:), pb1(:), pa2(:), pb2(:), dl(:)
            real(dp) :: lb1(2), lb2(2), dmg, dse
            integer  :: ss, gg, mm, nchg, k1
            logical  :: l_two, l_win2
            k1 = kstar
            ja = pack(ia, child(ia) > 0)
            jb = pack(ib, child(ib) > 0)
            allocate(cnt(k0,k1), source=0)
            do mm = 1, size(ja)
                cnt(cand(ja(mm)),child(ja(mm))) = cnt(cand(ja(mm)),child(ja(mm))) + 1
            end do
            do mm = 1, size(jb)
                cnt(cand(jb(mm)),child(jb(mm))) = cnt(cand(jb(mm)),child(jb(mm))) + 1
            end do
            allocate(smap(k0), source=0)
            do ss = 1, k0
                if( sum(cnt(ss,:)) > 0 ) smap(ss) = maxloc(cnt(ss,:), dim=1)
            end do
            l_two = all([(any(smap == gg), gg=1,k1)])
            if( k0 == k1 .and. .not. l_two )then
                smap  = [(ss, ss=1,k0)]
                l_two = .true.
            endif
            nchg   = 0
            l_win2 = .false.
            dmg    = 0._dp
            dse    = 0._dp
            if( l_two )then
                a1 = child(ja)
                b1 = child(jb)
                a2 = smap(cand(ja))
                b2 = smap(cand(jb))
                call self%cross_em(ja, jb, child(ja), child(jb), k1, eps, neq, a1, b1, l_full=.true., kcut=self%kt, &
                    &llband=lb1, niter=1, lpa=pa1, lpb=pb1)
                call self%cross_em(ja, jb, smap(cand(ja)), smap(cand(jb)), k1, eps, neq, a2, b2, l_full=.true., &
                    &kcut=self%kt, llband=lb2, niter=1, lpa=pa2, lpb=pb2)
                a2  = smap(cand(ja))
                b2  = smap(cand(jb))
                dl  = [pa2 - pa1, pb2 - pb1]
                dmg = sum(dl)
                dse = sqrt(real(size(dl),dp)) * sqrt(sum((dl - dmg / real(size(dl),dp))**2) / real(max(size(dl) - 1, 1),dp))
                l_win2 = dmg > 3._dp * dse
                if( l_win2 )then
                    nchg = count(child(ja) /= a2) + count(child(jb) /= b2)
                    child(ja) = a2
                    child(jb) = b2
                endif
            endif
            sz = [(count(child == gg), gg=1,k1)]
            allocate(rk(0:k1), source=0)
            do gg = 1, k1
                if( sz(gg) > 0 ) rk(gg) = count(sz > sz(gg)) + count(sz(1:gg-1) == sz(gg)) + 1
            end do
            child = rk(child)
            rkc   = rk
            kstar = count(sz > 0)
            if( l_two )then
                write(logfhandle,'(A)') 'membership by the choose rule, test band '//int2str(self%kt)//'..'//int2str(self%khi)//&
                    &': cross-half log-likelihood of the raw stage-B states minus the raw EM groups '//fmtr(dmg,1)//', SE '//&
                    &fmtr(dse,1)//' (sqrt(N) x the sd of '//int2str(size(dl))//' per-particle deltas), 3 SE '//fmtr(3._dp*dse,1)//&
                    &' -> '//trim(merge('stage-B states', 'EM groups     ', l_win2))//' (raw); '//int2str(nchg)//' particles change group'
            else
                write(logfhandle,'(A)') 'membership by the choose rule: the stage-B states do not cover the '//int2str(k1)//&
                    &' groups; the EM groups (raw) stand'
            endif
        end subroutine choose_members

        ! the cannot-link matrix over stage-B units becomes one over EM groups: linked when any of their start units were
        subroutine link_em_groups
            logical, allocatable :: cg(:,:)
            integer :: m1, m2
            allocate(cg(kem,kem), source=.false.)
            do m1 = 1, ng
                do m2 = 1, ng
                    if( m1 /= m2 ) cg(m1,m2) = self%cnl(gcom(m1),gcom(m2))
                end do
            end do
            call move_alloc(cg, self%cnl)
        end subroutine link_em_groups

    end subroutine test_once

    ! the pair decision on equal-N maps built once: every present group's maps at min(neq, its size) per half, docked onto the
    ! other half's consensus; T(p,q) in the test band from one FSC curve per pair; the smallest T merges while it sits below the
    ! null threshold (the null extended while it is within reach), the merged group redrawn at neq; tmin0 / tmin = the weakest
    ! pair at the start / among the kept groups
    subroutine decide_once( self, ia, ib, la, lb, kmax, eps, neq, ndraw, cons, nl, plist, gof, exa, exb, kout, tmin0, tmin )
        class(state_mergeup), intent(inout) :: self
        integer,              intent(in)    :: ia(:), ib(:), la(:), lb(:), kmax, neq, ndraw
        real,                 intent(in)    :: eps(0:)
        complex,              intent(in)    :: cons(:,:,:,:)
        type(shuffle_null),   intent(inout) :: nl
        integer, allocatable, intent(out)   :: plist(:), gof(:), exa(:), exb(:)
        integer,              intent(out)   :: kout
        real,                 intent(out)   :: tmin0, tmin
        complex, allocatable :: allm(:,:,:,:)
        integer, allocatable :: presa(:), presb(:), jidx(:), joff(:), mapa(:), mapb(:), groupof(:), nsa(:), nsb(:)
        real,    allocatable :: jscale(:), rmat(:,:), merged(:)
        real(dp),allocatable :: s2(:,:)
        logical, allocatable :: l_alive(:)
        real    :: best, mu, sd, thr
        logical :: l_merge
        integer :: g, npair, p, q, bp, bq, mm
        presa = pack([(g, g=1,kmax)], [(count(la == g) > 0, g=1,kmax)])
        presb = pack([(g, g=1,kmax)], [(count(lb == g) > 0, g=1,kmax)])
        allocate(jidx(0), joff(1), jscale(0), nsa(0), nsb(0))
        joff(1) = 1
        call add_jobs(ia, la, presa, nsa)
        call add_jobs(ib, lb, presb, nsb)
        allocate(allm(self%nh,self%box,self%box,size(jscale)))
        call self%maps(jidx, joff, jscale, eps, .true., .false., allm)
        !$omp parallel do default(shared) private(mm) schedule(dynamic) proc_bind(close)
        do mm = 1, size(jscale)
            call self%dock(allm(:,:,:,mm), cons(:,:,:,merge(2,1,mm <= size(presa))), self%kt, self%khi)
        end do
        !$omp end parallel do
        plist = pack(presa, [(any(presb == presa(g)), g=1,size(presa))])
        npair = size(plist)
        allocate(mapa(npair), mapb(npair), s2(0:self%nyq,npair))
        do mm = 1, npair
            mapa(mm) = findloc(presa, plist(mm), dim=1)
            mapb(mm) = size(presa) + findloc(presb, plist(mm), dim=1)
        end do
        ! every paired group's per-particle noise power per shell from its half difference (nsa / nsb follow presa / presb)
        do mm = 1, npair
            call self%noise_power(allm(:,:,:,mapa(mm)), allm(:,:,:,mapb(mm)), nsa(mapa(mm)), nsb(mapb(mm)-size(presa)), s2(:,mm))
        end do
        allocate(exa(0:kmax), exb(0:kmax), source=0)
        call extras(presa, 0, 1, exa)
        call extras(presb, size(presa), 2, exb)
        groupof = [(p, p=1,npair)]
        allocate(rmat(npair,npair), source=huge(1.))
        !$omp parallel do default(shared) private(p,q) schedule(dynamic) collapse(2) proc_bind(close)
        do p = 1, npair
            do q = 1, npair
                if( q > p ) rmat(p,q) = self%pair_stat_df(allm(:,:,:,mapa(p)), allm(:,:,:,mapa(q)), allm(:,:,:,mapb(p)), &
                    &allm(:,:,:,mapb(q)), s2(:,p), s2(:,q), ndraw, self%kt)
            end do
        end do
        !$omp end parallel do
        tmin0 = huge(tmin0)
        do p = 1, npair
            do q = p + 1, npair
                tmin0 = min(tmin0, rmat(p,q))
            end do
        end do
        if( npair < 2 ) tmin0 = 0.
        write(logfhandle,'(A)') '[single pass, draw-free] full-N maps, statistic at the common n '//int2str(ndraw)//'; test band '//&
            &int2str(self%kt)//'..'//int2str(self%khi)//'; paired '//int2str(npair)//' groups; pair statistics at the start '//&
            &'(sorted): '//sorted_upper(rmat, npair)
        allocate(l_alive(npair), source=.true.)
        allocate(merged(0))
        do
            best = huge(best)
            bp   = 0
            bq   = 0
            do p = 1, npair
                if( .not. l_alive(p) ) cycle
                do q = p+1, npair
                    if( .not. l_alive(q) ) cycle
                    if( rmat(p,q) < best )then
                        best = rmat(p,q)
                        bp   = p
                        bq   = q
                    endif
                end do
            end do
            if( bp == 0 ) exit
            call self%null_resolve(nl, best, eps, cons, l_merge)
            if( .not. l_merge ) exit
            if( self%l_equiv )then
                if( .not. shown_equal(bp, bq, best) )then
                    rmat(bp,bq) = huge(best)
                    cycle
                endif
            endif
            merged = [merged, best]
            l_alive(bq) = .false.
            where( groupof == bq ) groupof = bp
            ! the merged group's full-N maps are the particle-weighted means of its members' (no redraw)
            allm(:,:,:,mapa(bp)) = (real(nsa(mapa(bp))) * allm(:,:,:,mapa(bp)) + real(nsa(mapa(bq))) * allm(:,:,:,mapa(bq))) / &
                &real(nsa(mapa(bp)) + nsa(mapa(bq)))
            allm(:,:,:,mapb(bp)) = (real(nsb(mapb(bp)-size(presa))) * allm(:,:,:,mapb(bp)) + real(nsb(mapb(bq)-size(presa))) * &
                &allm(:,:,:,mapb(bq))) / real(nsb(mapb(bp)-size(presa)) + nsb(mapb(bq)-size(presa)))
            nsa(mapa(bp)) = nsa(mapa(bp)) + nsa(mapa(bq))
            nsb(mapb(bp)-size(presa)) = nsb(mapb(bp)-size(presa)) + nsb(mapb(bq)-size(presa))
            call self%noise_power(allm(:,:,:,mapa(bp)), allm(:,:,:,mapb(bp)), nsa(mapa(bp)), nsb(mapb(bp)-size(presa)), s2(:,bp))
            do q = 1, npair
                if( q == bp .or. .not. l_alive(q) ) cycle
                rmat(min(bp,q),max(bp,q)) = self%pair_stat_df(allm(:,:,:,mapa(bp)), allm(:,:,:,mapa(q)), allm(:,:,:,mapb(bp)), &
                    &allm(:,:,:,mapb(q)), s2(:,bp), s2(:,q), ndraw, self%kt)
            end do
        end do
        allocate(gof(npair))
        do mm = 1, npair
            gof(mm) = count(l_alive(1:groupof(mm)))
        end do
        kout = count(l_alive)
        tmin = huge(tmin)
        do p = 1, npair
            if( .not. l_alive(p) ) cycle
            do q = p + 1, npair
                if( l_alive(q) ) tmin = min(tmin, rmat(p,q))
            end do
        end do
        if( kout < 2 ) tmin = tmin0
        call null_thr(nl, mu, sd, thr)
        write(logfhandle,'(A)') int2str(size(merged))//' merges at '//joinedr(merged)//' -> '//int2str(kout)//' groups'
        deallocate(allm, rmat, l_alive, groupof, merged, mapa, mapb, jidx, joff, jscale, s2, nsa, nsb)

    contains

        ! one-sided equivalence at the margin: r_eq = r_hat(stat + 1.645 sd) at the statistic's n <= mdelta, sd the null's own
        logical function shown_equal( p, q, t )
            integer, intent(in) :: p, q
            real,    intent(in) :: t
            real    :: sg, y, req
            integer :: i, j
            logical :: lnk
            lnk = .false.
            if( allocated(self%cnl) )then
                do i = 1, npair
                    if( groupof(i) /= p ) cycle
                    do j = 1, npair
                        if( groupof(j) == q ) lnk = lnk .or. self%cnl(plist(i), plist(j))
                    end do
                end do
            endif
            call null_thr(nl, mu, sd, thr)
            sg  = max(sd, 1.e-3)
            y   = min(max((t + 1.645 * sg - mu) / (1. - mu), 1.e-4), 0.995)
            req = 2. * y / ((1. - y) * real(max(ndraw, 1)))
            shown_equal = (.not. lnk) .and. req <= self%mdelta
            write(logfhandle,'(A)') 'equivalence: groups '//int2str(p)//' '//int2str(q)//' stat '//fmtr(real(t,dp),3)//' n '//&
                &int2str(ndraw)//' null sd '//fmtr(real(sg,dp),3)//' r_eq '//fmtr(real(req,dp),6)//' margin '//&
                &fmtr(real(self%mdelta,dp),6)//' -> '//trim(merge('merge             ', merge('keep (cannot-link)', &
                &'keep (unresolved) ', lnk), shown_equal))
        end function shown_equal

        ! every particle of every present group, one job per group, eps scaled to the group's size over the consensus N
        subroutine add_jobs( idx, lab, pres, ns )
            integer,              intent(in)    :: idx(:), lab(:), pres(:)
            integer, allocatable, intent(inout) :: ns(:)
            integer, allocatable :: sel(:)
            integer :: gg
            do gg = 1, size(pres)
                sel    = pack(idx, lab == pres(gg))
                jidx   = [jidx, sel]
                joff   = [joff, size(jidx) + 1]
                jscale = [jscale, real(size(sel)) / real(neq)]
                ns     = [ns, size(sel)]
            end do
        end subroutine add_jobs

        subroutine extras( pres, off, hh, ex )
            integer, intent(in)    :: pres(:), off, hh
            integer, intent(inout) :: ex(0:)
            real    :: cc, bestc
            integer :: gg, mmm, ipair
            do gg = 1, size(pres)
                if( any(plist == pres(gg)) ) cycle
                bestc = -huge(bestc)
                do mmm = 1, npair
                    ipair = findloc(pres, plist(mmm), dim=1)
                    cc = self%shell_corr(allm(:,:,:,off+gg) - cons(:,:,:,hh), allm(:,:,:,off+ipair) - cons(:,:,:,hh), &
                        &2, self%klo - 2)
                    if( cc > bestc )then
                        bestc = cc
                        ex(pres(gg)) = plist(mmm)
                    endif
                end do
            end do
        end subroutine extras

    end subroutine decide_once

    ! the quality-split flag (reporting only): each delivered group's mean of three per-particle quality terms, the contrast
    ! and zero-contrast posterior of fit_contrast, the log noise variance the planes were whitened with and the stage-B
    ! correlation, each as a standardised offset d from the node mean (unit: the pooled within-group sd); |d| >= 1 flags
    subroutine quality_flag( self, child, kstar, lines, l_flag )
        class(state_mergeup), intent(inout) :: self
        integer,              intent(in)    :: child(:), kstar
        type(string), allocatable, intent(out) :: lines(:)
        logical,              intent(out)   :: l_flag
        real, allocatable :: x(:,:), mp(:)
        real    :: mc(3), sw(3), mall(3), d(3)
        integer :: g, n, i, k
        character(len=:), allocatable :: ln
        l_flag = .false.
        allocate(lines(0))
        if( .not. allocated(self%cam) ) call self%fit_contrast
        if( kstar < 2 ) return
        allocate(x(self%n,3))
        x(:,1) = self%cam
        x(:,2) = log(max(self%nv, 1.e-12))
        x(:,3) = self%qcorr
        sw = 0.
        n  = 0
        do g = 1, kstar
            if( count(child == g) == 0 ) cycle
            do k = 1, 3
                mc(k) = sum(x(:,k), child == g) / real(count(child == g))
                sw(k) = sw(k) + sum((x(:,k) - mc(k))**2, child == g)
            end do
            n = n + count(child == g)
        end do
        sw   = sqrt(sw / real(max(n - kstar, 1)))
        mall = [(sum(x(:,k), child > 0) / real(max(count(child > 0), 1)), k=1,3)]
        allocate(mp(kstar), source=0.)
        do g = 1, kstar
            if( count(child == g) == 0 ) cycle
            do k = 1, 3
                mc(k) = sum(x(:,k), child == g) / real(count(child == g))
                d(k)  = (mc(k) - mall(k)) / max(sw(k), 1.e-12)
            end do
            mp(g)  = sum(self%cpj, child == g) / real(count(child == g))
            l_flag = l_flag .or. any(abs(d) >= 1.)
            ln = 'group '//int2str(g)//' n '//int2str(count(child == g))//' contrast '//fmtr(real(mc(1),dp),4)//' d '//&
                &fmtr(real(d(1),dp),3)//' | log noise var '//fmtr(real(mc(2),dp),4)//' d '//fmtr(real(d(2),dp),3)//&
                &' | stage-B corr '//fmtr(real(mc(3),dp),4)//' d '//fmtr(real(d(3),dp),3)//' | p_empty '//fmtr(real(mp(g),dp),4)
            lines = [lines, string(ln)]
        end do
        lines = [lines, string('quality split flag '//trim(merge('yes', 'no ', l_flag))//' (within-group sd: contrast '//&
            &fmtr(real(sw(1),dp),4)//', log noise var '//fmtr(real(sw(2),dp),4)//', corr '//fmtr(real(sw(3),dp),4)//&
            &'; |d| >= 1 on any term flags)')]
        do i = 1, size(lines)
            write(logfhandle,'(A)') 'quality: '//lines(i)%to_char()
        end do
    end subroutine quality_flag

    ! report only: the rotation from the node's input pose to the stage-B pose (over the point group), the share that turned more
    ! than POSE_FLIP degrees, and per pair of children (stage-B states where none) the flip-share difference against strata shuffles
    subroutine pose_lines( self, lab, child, kstar, eul_in, pgrp, lines )
        class(state_mergeup),      intent(inout) :: self
        integer,                   intent(in)    :: lab(:), child(:), kstar
        real,         allocatable, intent(in)    :: eul_in(:,:)
        character(len=*),          intent(in)    :: pgrp
        type(string), allocatable, intent(out)   :: lines(:)
        real,    parameter :: POSE_FLIP = 30., DFLIP_MIN = 0.04, DFLIP_Z = 5.
        integer, parameter :: NSHUF_POSE = 200
        type(sym) :: pg
        type(ori) :: o1, o2, o3
        real,    allocatable :: ang(:), dn(:)
        integer, allocatable :: idx(:), st(:), grp(:), pl(:)
        logical, allocatable :: flip(:)
        character(len=512) :: ln
        real    :: amin, fa, fb, fa0, fb0, d, mu, sd, z
        integer :: i, k, g, h, nsym, na, nb, is, ng
        character(len=6) :: kind
        allocate(lines(0))
        if( .not. allocated(eul_in) )then
            lines = [string('# pose: the node''s input poses (../n.simple) were not found, no pose statistics')]
            return
        endif
        call pg%new(pgrp)
        nsym = pg%get_nsym()
        call o1%new(.false.)
        call o2%new(.false.)
        call o3%new(.false.)
        allocate(ang(self%n))
        do i = 1, self%n
            call o1%set_euler(eul_in(:,self%rows(i)))
            call o2%set_euler(self%eul(:,i))
            amin = huge(amin)
            do k = 1, nsym
                call pg%apply(o2, k, o3)
                amin = min(amin, o1%geodesic_dist_trace(o3))
            end do
            ang(i) = amin * 180. / PI
        end do
        flip = ang > POSE_FLIP
        write(ln,'(A,F8.5,A,I0,A,I0,A,F8.3)') '# pose node flip30 ', real(count(flip)) / real(self%n), ' nact ', self%n, &
            &' nsym ', nsym, ' mean_angle ', sum(ang) / real(self%n)
        lines = [lines, string(trim(ln))]
        call o1%kill
        call o2%kill
        call o3%kill
        call pg%kill
        if( kstar >= 2 )then
            grp  = child
            ng   = kstar
            kind = 'split '
        else
            grp  = lab
            ng   = maxval(lab)
            kind = 'states'
        endif
        allocate(dn(NSHUF_POSE))
        do g = 1, ng - 1
            do h = g + 1, ng
                na = count(grp == g)
                nb = count(grp == h)
                if( min(na, nb) < 50 ) cycle
                idx = pack([(i, i=1,self%n)], grp == g .or. grp == h)
                pl  = merge(1, 2, grp(idx) == g)
                st  = self%strata(idx)
                d   = dflip(pl)
                fa0 = fa
                fb0 = fb
                do is = 1, NSHUF_POSE
                    dn(is) = dflip(sshuf(pl, st))
                end do
                mu = sum(dn) / real(NSHUF_POSE)
                sd = sqrt(sum((dn - mu)**2) / real(NSHUF_POSE - 1))
                z  = (d - mu) / max(sd, 1.e-9)
                write(ln,'(A,I0,A,I0,A,I0,1X,I0,A,F7.4,1X,F7.4,A,F7.4,A,F8.2,A,F7.4,1X,F7.4,A)') '# pose '//kind//' pair ', g, '|', h, ' n ', na, nb, &
                    &' flips ', fa0, fb0, ' dflip ', d, ' z ', z, ' null ', mu, sd, ' flag '//merge('yes', 'no ', &
                    &d >= DFLIP_MIN .and. z >= DFLIP_Z)
                lines = [lines, string(trim(ln))]
            end do
        end do

    contains

        real function dflip( l )
            integer, intent(in) :: l(:)
            fa = real(count(flip(idx) .and. l == 1)) / real(max(1, count(l == 1)))
            fb = real(count(flip(idx) .and. l == 2)) / real(max(1, count(l == 2)))
            dflip = abs(fa - fb)
        end function dflip

    end subroutine pose_lines

    subroutine consensus( self, ia, ib, la, lb, k0, neq, eps, cons )
        class(state_mergeup), intent(inout) :: self
        integer,              intent(in)    :: ia(:), ib(:), la(:), lb(:), k0
        integer,              intent(out)   :: neq
        real,    allocatable, intent(out)   :: eps(:)
        complex, allocatable, intent(inout) :: cons(:,:,:,:)
        real,    allocatable :: epsa(:), epsb(:)
        integer, allocatable :: js(:)
        neq = min(minnz(la, k0), minnz(lb, k0))
        call self%accumulate_jobs([ia, ib], [1, size(ia)+1, size(ia)+size(ib)+1], .false., js)
        call self%recs(js(1))%get_rho_shell_medians(epsa)
        call self%recs(js(2))%get_rho_shell_medians(epsb)
        eps = EPS_REG * 0.5 * (epsa * real(neq) / real(size(ia)) + epsb * real(neq) / real(size(ib)))
        if( .not. allocated(cons) ) allocate(cons(self%nh,self%box,self%box,2))
        call self%finish(js(1), eps * real(size(ia)) / real(neq), .true., cons(:,:,:,1))
        call self%finish(js(2), eps * real(size(ib)) / real(neq), .true., cons(:,:,:,2))
    end subroutine consensus

    subroutine screen( self, ia, ib, la, lb, k0, neq, eps, cons, junk )
        class(state_mergeup), intent(inout) :: self
        integer, allocatable, intent(inout) :: ia(:), ib(:), la(:), lb(:)
        integer,              intent(in)    :: k0
        integer,              intent(inout) :: neq
        real,    allocatable, intent(inout) :: eps(:)
        complex, allocatable, intent(inout) :: cons(:,:,:,:)
        integer, allocatable, intent(out)   :: junk(:)
        complex, allocatable :: um(:,:,:,:), sc(:,:,:,:)
        integer, allocatable :: units(:), jidx(:), joff(:), sa(:), sb(:), part(:), js(:)
        real,    allocatable :: jscale(:), met(:,:,:)
        logical, allocatable :: flag(:), keepa(:), keepb(:)
        character(len=:), allocatable :: tail
        integer :: nu, u, nfl, i
        logical :: lkeep
        units = pack([(u, u=1,k0)], [(count(la == u) > 0 .and. count(lb == u) > 0, u=1,k0)])
        nu    = size(units)
        allocate(junk(0))
        if( nu == 0 ) return
        allocate(jidx(0), joff(1), jscale(0))
        joff(1) = 1
        do u = 1, nu
            call add_all(ia, la, units(u))
        end do
        do u = 1, nu
            call add_all(ib, lb, units(u))
        end do
        allocate(um(self%nh,self%box,self%box,2*nu), met(3,3,nu), flag(nu))
        call self%maps(jidx, joff, jscale, eps, .true., .false., um)
        call flags(cons, 0.5 * real(size(ia) + size(ib)))
        nfl = count(flag)
        if( nfl > 0 .and. nfl < nu )then
            sa = pack(ia, [(any(units == la(i) .and. .not. flag), i=1,size(ia))])
            sb = pack(ib, [(any(units == lb(i) .and. .not. flag), i=1,size(ib))])
            allocate(sc(self%nh,self%box,self%box,2))
            call self%accumulate_jobs([sa, sb], [1, size(sa)+1, size(sa)+size(sb)+1], .false., js)
            call self%finish(js(1), eps * real(size(sa)) / real(neq), .true., sc(:,:,:,1))
            call self%finish(js(2), eps * real(size(sb)) / real(neq), .true., sc(:,:,:,2))
            call flags(sc, 0.5 * real(size(sa) + size(sb)))
        endif
        junk = pack(units, flag)
        part = pack(units, .not. flag .and. min(met(2,2,:), met(2,3,:)) < JUNK_PROJ)
        ! junk_screen=no reports what the screen flags and sets nothing aside
        lkeep = self%p%junk_screen .eq. 'no'
        do u = 1, nu
            tail = ''
            if( flag(u) .and. lkeep ) tail = ' (kept: junk_screen=no)'
            write(logfhandle,'(A)') '  junk screen unit '//int2str(units(u))//': projection on the screened consensus '//&
                &fmtr(real(met(2,2,u),dp),2)//' / '//fmtr(real(met(2,3,u),dp),2)//', rel. SSNR '//&
                &fmtr(real(met(3,2,u),dp),2)//' / '//fmtr(real(met(3,3,u),dp),2)//' -> '//merge('JUNK','ok  ',flag(u))//tail
        end do
        if( size(part) > 0 ) write(logfhandle,'(A)') '  partial units (projection < '//fmtr(real(JUNK_PROJ,dp),2)//&
            &' with the node signal; tested as composition, not junk): '//joined(part)
        if( lkeep )then
            deallocate(junk)
            allocate(junk(0))
        endif
        if( size(junk) > 0 .and. size(junk) < nu )then
            keepa = [(.not. any(junk == la(i)), i=1,size(la))]
            keepb = [(.not. any(junk == lb(i)), i=1,size(lb))]
            ia = pack(ia, keepa)
            la = pack(la, keepa)
            ib = pack(ib, keepb)
            lb = pack(lb, keepb)
            call self%consensus(ia, ib, la, lb, k0, neq, eps, cons)
            write(logfhandle,'(A)') 'junk screen: units '//joined(junk)//' set aside; the split test runs on '//&
                &int2str(size(ia))//' / '//int2str(size(ib))//' particles, equal N '//int2str(neq)
        endif

    contains

        subroutine add_all( idx, lab, uu )
            integer, intent(in) :: idx(:), lab(:), uu
            integer, allocatable :: sel(:)
            sel    = pack(idx, lab == uu)
            jidx   = [jidx, sel]
            joff   = [joff, size(jidx) + 1]
            jscale = [jscale, real(max(size(sel),1)) / real(neq)]
        end subroutine add_all

        subroutine flags( c, ncv )
            complex, intent(in) :: c(:,:,:,:)
            real,    intent(in) :: ncv
            integer :: uu
            do uu = 1, nu
                call unit_metrics(um(:,:,:,uu), um(:,:,:,nu+uu), c(:,:,:,1), c(:,:,:,2), &
                    &0.5 * real(count(la == units(uu)) + count(lb == units(uu))), ncv, met(:,:,uu))
                flag(uu) = met(3,3,uu) < JUNK_SSNR
            end do
        end subroutine flags

        subroutine unit_metrics( xa, xb, ca, cb, ng, nc, mt )
            complex, intent(in)  :: xa(:,:,:), xb(:,:,:), ca(:,:,:), cb(:,:,:)
            real,    intent(in)  :: ng, nc
            real,    intent(out) :: mt(3,3)
            real(dp) :: pg(3), pc(3), xc(3), aa(3), bb(3), ra(3), rb(3), w, fg, fc
            integer  :: h, k, l, sh, ib_, lo(3), hi(3)
            lo = [2, self%klo, self%kb + 1]
            hi = [self%klo - 1, self%kb, self%khi]
            pg = 0._dp
            pc = 0._dp
            xc = 0._dp
            aa = 0._dp
            bb = 0._dp
            ra = 0._dp
            rb = 0._dp
            do l = 1, self%box
                do k = 1, self%box
                    do h = 1, self%nh
                        sh = self%fsh(h,k,l)
                        do ib_ = 1, 3
                            if( sh < lo(ib_) .or. sh > hi(ib_) ) cycle
                            w = real(self%fw(h,k,l),dp)
                            pg(ib_) = pg(ib_) + w * real(xb(h,k,l)*conjg(xa(h,k,l)),dp)
                            pc(ib_) = pc(ib_) + w * real(cb(h,k,l)*conjg(ca(h,k,l)),dp)
                            xc(ib_) = xc(ib_) + w * 0.5_dp * (real(cb(h,k,l)*conjg(xa(h,k,l)),dp) + &
                                &real(ca(h,k,l)*conjg(xb(h,k,l)),dp))
                            aa(ib_) = aa(ib_) + w * real(abs(xa(h,k,l)),dp)**2
                            bb(ib_) = bb(ib_) + w * real(abs(xb(h,k,l)),dp)**2
                            ra(ib_) = ra(ib_) + w * real(abs(ca(h,k,l)),dp)**2
                            rb(ib_) = rb(ib_) + w * real(abs(cb(h,k,l)),dp)**2
                        end do
                    end do
                end do
            end do
            do ib_ = 1, 3
                fg = pg(ib_) / max(sqrt(aa(ib_) * bb(ib_)), 1.e-30_dp)
                fc = pc(ib_) / max(sqrt(ra(ib_) * rb(ib_)), 1.e-30_dp)
                mt(1,ib_) = real(pg(ib_) / sign(max(abs(pc(ib_)), 1.e-30_dp), pc(ib_)))
                mt(2,ib_) = real(xc(ib_) / sign(max(abs(pc(ib_)), 1.e-30_dp), pc(ib_)))
                mt(3,ib_) = real((max(fg,0._dp) / max(1._dp - fg, 1.e-3_dp) / real(ng,dp)) / &
                    &max(fc / max(1._dp - fc, 1.e-3_dp) / real(nc,dp), 1.e-12_dp))
            end do
        end subroutine unit_metrics

    end subroutine screen

    ! the shuffle null of the draw-free statistic: every shuffled group's maps from all its particles, the statistic at the
    ! common n nl%neq; nul2 holds the test-band values
    subroutine null_upto_df( self, nl, n, eps, cons )
        class(state_mergeup), intent(inout) :: self
        type(shuffle_null),   intent(inout) :: nl
        integer,              intent(in)    :: n
        real,                 intent(in)    :: eps(0:)
        complex,              intent(in)    :: cons(:,:,:,:)
        complex,  allocatable :: nm(:,:,:,:)
        real(dp), allocatable :: s2(:,:)
        integer,  allocatable :: joff(:), nsz(:)
        real,     allocatable :: jscale(:), vals(:)
        integer :: nbatch, s1, s2b, nb, nj, npr, ipr, s, ii, g, q, j, ng, j0
        ng     = nl%ng
        npr    = ng * (ng - 1) / 2
        nbatch = max(1, min(nl%nsh, MAXJOBS / max(1, 2*ng)))
        do s1 = nl%ndone + 1, n, nbatch
            s2b = min(n, s1 + nbatch - 1)
            nb  = s2b - s1 + 1
            nj  = 2 * ng * nb
            j0  = (s1 - 1) * 2 * ng
            joff   = [(nl%soff(j0+j) - nl%soff(j0) + 1, j=0,nj)]
            nsz    = [(nl%soff(j0+j) - nl%soff(j0+j-1), j=1,nj)]
            jscale = [(real(nsz(j)) / real(nl%neqref), j=1,nj)]
            allocate(nm(self%nh,self%box,self%box,nj), s2(0:self%nyq,nj))
            call self%maps(nl%sjidx(nl%soff(j0)+1:nl%soff(j0+nj)), joff, jscale, eps, .true., .false., nm)
            !$omp parallel do default(shared) private(j) schedule(dynamic) proc_bind(close)
            do j = 1, nj
                call self%dock(nm(:,:,:,j), cons(:,:,:,2-mod((j-1)/ng,2)), self%kt, self%khi)
            end do
            !$omp end parallel do
            !$omp parallel do default(shared) private(s,g) schedule(dynamic) collapse(2) proc_bind(close)
            do s = 0, nb - 1
                do g = 1, ng
                    call self%noise_power(nm(:,:,:,s*2*ng+g), nm(:,:,:,s*2*ng+ng+g), nsz(s*2*ng+g), nsz(s*2*ng+ng+g), s2(:,s*2*ng+g))
                end do
            end do
            !$omp end parallel do
            allocate(vals(npr*nb))
            !$omp parallel do default(shared) private(ipr,s,g,q,ii) schedule(dynamic) proc_bind(close)
            do ipr = 1, npr*nb
                s  = (ipr - 1) / npr
                ii = ipr - s * npr
                call pair_index(ii, ng, g, q)
                vals(ipr) = self%pair_stat_df(nm(:,:,:,s*2*ng+g), nm(:,:,:,s*2*ng+q), nm(:,:,:,s*2*ng+ng+g), &
                    &nm(:,:,:,s*2*ng+ng+q), s2(:,s*2*ng+g), s2(:,s*2*ng+q), nl%neq, self%kt)
            end do
            !$omp end parallel do
            nl%nul2 = [nl%nul2, vals]
            deallocate(nm, s2, vals)
        end do
        nl%ndone = max(nl%ndone, n)
    end subroutine null_upto_df

    ! a group's per-particle noise power per shell from its half difference: |A - B|^2 summed over the shell, over 1/nA + 1/nB
    subroutine noise_power( self, a, b, na, nb, s2 )
        class(state_mergeup), intent(in)  :: self
        complex,              intent(in)  :: a(:,:,:), b(:,:,:)
        integer,              intent(in)  :: na, nb
        real(dp),             intent(out) :: s2(0:)
        integer :: h, k, l, sh
        s2 = 0._dp
        do l = 1, self%box
            do k = 1, self%box
                do h = 1, self%nh
                    sh = self%fsh(h,k,l)
                    s2(sh) = s2(sh) + self%fw(h,k,l) * real(abs(a(h,k,l) - b(h,k,l)),dp)**2
                end do
            end do
        end do
        s2 = s2 / (1._dp / real(max(na,1),dp) + 1._dp / real(max(nb,1),dp))
    end subroutine noise_power

    ! the draw-free pair statistic: per shell the between-group signal from the cross-half product of the full-N difference maps
    ! (the noise cancels), the noise of a difference of two n-particle maps from the groups' per-particle noise powers, the
    ! expected cross-half FSC at n, and the peak of its 3-shell means over the band
    real function pair_stat_df( self, ap, aq, bp, bq, s2p, s2q, n, kl )
        class(state_mergeup), intent(in) :: self
        complex,              intent(in) :: ap(:,:,:), aq(:,:,:), bp(:,:,:), bq(:,:,:)
        real(dp),             intent(in) :: s2p(0:), s2q(0:)
        integer,              intent(in) :: n, kl
        real(dp) :: c(0:self%nyq), nn(0:self%nyq)
        real     :: f(0:self%nyq)
        complex  :: da, db
        integer  :: h, k, l, sh
        c = 0._dp
        do l = 1, self%box
            do k = 1, self%box
                do h = 1, self%nh
                    sh = self%fsh(h,k,l)
                    da = ap(h,k,l) - aq(h,k,l)
                    db = bp(h,k,l) - bq(h,k,l)
                    c(sh) = c(sh) + self%fw(h,k,l) * (real(real(da),dp)*real(real(db),dp) + real(aimag(da),dp)*real(aimag(db),dp))
                end do
            end do
        end do
        nn = (s2p + s2q) / real(max(n,1),dp)
        f  = real(c / max(max(c, 0._dp) + nn, 1.e-30_dp))
        pair_stat_df = -huge(1.)
        do k = kl, max(self%khi - 2, kl)
            pair_stat_df = max(pair_stat_df, sum(f(k:k+2)) / 3.)
        end do
    end function pair_stat_df

    ! more shuffles only while r lies within NSHUF_KSIG spreads of the threshold the full nsh shuffles could still reach
    subroutine null_resolve( self, nl, r, eps, cons, l_merge )
        class(state_mergeup), intent(inout) :: self
        type(shuffle_null),   intent(inout) :: nl
        real,                 intent(in)    :: r, eps(0:)
        complex,              intent(in)    :: cons(:,:,:,:)
        logical,              intent(out)   :: l_merge
        real :: mu, sd, thr
        do
            call null_thr(nl, mu, sd, thr)
            if( nl%ndone >= nl%nsh ) exit
            if( abs(r - thr) > NSHUF_KSIG * sd * null_spread(nl%ndone, nl%nsh, nl%zc) ) exit
            call self%null_upto_df(nl, min(nl%nsh, nl%ndone + nl%nstep), eps, cons)
        end do
        l_merge = r < thr
    end subroutine null_resolve

    !> l_site: the site split's design, masked references (the focus mask) and half B's groups renumbered to half A's
    !! by the cross-half correlation of their deviations before every iteration
    subroutine cross_em( self, ia, ib, starta, startb, k, eps, neq, la, lb, l_full, ll, kcut, llband, niter, lpa, lpb, l_site )
        class(state_mergeup), intent(inout) :: self
        integer,              intent(in)    :: ia(:), ib(:), starta(:), startb(:), k, neq
        real,                 intent(in)    :: eps(0:)
        integer,              intent(inout) :: la(:), lb(:)
        logical,  optional,   intent(in)    :: l_full
        real(dp), optional,   intent(out)   :: ll
        integer,  optional,   intent(in)    :: kcut
        real(dp), optional,   intent(out)   :: llband(2)
        integer,  optional,   intent(in)    :: niter
        real(dp), allocatable, optional, intent(out) :: lpa(:), lpb(:)
        logical,  optional,   intent(in)    :: l_site
        complex,             allocatable :: ma(:,:,:,:), mb(:,:,:,:)
        type(reconstructor), allocatable :: refa(:), refb(:)
        real,    allocatable :: wa(:,:), wb(:,:), wan(:,:), wbn(:,:), fl(:,:)
        integer, allocatable :: mpa(:), mpb(:), ca(:), cb(:), perm(:)
        logical, allocatable :: l_alive(:)
        character(len=:), allocatable :: line
        real(dp) :: lla, llb, lba(2), lbb(2), cc(k,k), cbest
        real     :: chga, chgb
        integer  :: it, g, i, kband, nit, h, gb
        logical  :: l_fb, l_st
        l_fb = .false.
        if( present(l_full) ) l_fb = l_full
        l_st = .false.
        if( present(l_site) ) l_st = l_site
        kband = merge(self%khi, self%kb, l_fb)
        nit   = NEM_ITERS
        if( present(niter) ) nit = niter
        allocate(wa(size(ia),k), wb(size(ib),k), source=0.)
        do i = 1, size(ia)
            if( starta(i) > 0 ) wa(i,starta(i)) = 1.
        end do
        do i = 1, size(ib)
            if( startb(i) > 0 ) wb(i,startb(i)) = 1.
        end do
        allocate(l_alive(k), source=.true.)
        allocate(ma(self%nh,self%box,self%box,k+1), mb(self%nh,self%box,self%box,k+1), fl(0:self%nyq,k))
        allocate(refa(k), refb(k))
        do g = 1, k
            call refa(g)%new_accumulator(self%p, self%b%spproj, expand=.true., wthreads=.false.)
            call refb(g)%new_accumulator(self%p, self%b%spproj, expand=.true., wthreads=.false.)
        end do
        write(logfhandle,'(A)') 'cross-half EM: '//int2str(k)//' groups, '//int2str(nit)//' iterations, shells 2..'//&
            &int2str(kband)//', '//int2str(merge(self%nfb, self%nem, l_fb))//' plane samples per particle'
        do it = 1, nit
            if( self%p%merge_insert_audit == 'yes' )then
                self%audit_calls = self%audit_calls + 1
                call self%em_maps(ia, wa, eps, neq, l_fb, ma, 'cross_'//int2str(self%audit_calls)//'_a')
                call self%em_maps(ib, wb, eps, neq, l_fb, mb, 'cross_'//int2str(self%audit_calls)//'_b')
            else
                call self%em_maps(ia, wa, eps, neq, l_fb, ma)
                call self%em_maps(ib, wb, eps, neq, l_fb, mb)
            endif
            if( l_st )then
                call self%site_refs(ma, wa)
                call self%site_refs(mb, wb)
            endif
            if( l_st .and. k > 1 )then
                ! half B's groups renumbered to half A's: greedy on the cross-half correlation of the group deviations
                do g = 1, k
                    do h = 1, k
                        cc(g,h) = self%shell_dot(ma(:,:,:,g) - ma(:,:,:,k+1), mb(:,:,:,h) - mb(:,:,:,k+1), 2, kband) / &
                            &sqrt(max(self%shell_dot(ma(:,:,:,g) - ma(:,:,:,k+1), ma(:,:,:,g) - ma(:,:,:,k+1), 2, kband) * &
                            &self%shell_dot(mb(:,:,:,h) - mb(:,:,:,k+1), mb(:,:,:,h) - mb(:,:,:,k+1), 2, kband), 1.e-300_dp))
                    end do
                end do
                allocate(perm(k), source=0)
                do g = 1, k
                    cbest = -huge(1._dp)
                    gb    = 0
                    do h = 1, k
                        if( any(perm == h) ) cycle
                        if( cc(g,h) > cbest )then
                            cbest = cc(g,h)
                            gb    = h
                        endif
                    end do
                    perm(g) = gb
                end do
                if( any(perm /= [(g, g=1,k)]) )then
                    write(logfhandle,'(A)') '  half B groups renumbered to half A: '//joined(perm)
                    mb(:,:,:,1:k) = mb(:,:,:,perm)
                    wb(:,1:k)     = wb(:,perm)
                endif
                deallocate(perm)
            endif
            ! no docking on the site path: a shift per group map turns the two references into a consensus-gradient dipole
            if( .not. l_st )then
                !$omp parallel do default(shared) private(g) schedule(dynamic) proc_bind(close)
                do g = 1, 2*k
                    if( g <= k )then
                        call self%dock(ma(:,:,:,g), mb(:,:,:,k+1), 2, kband)
                    else
                        call self%dock(mb(:,:,:,g-k), ma(:,:,:,k+1), 2, kband)
                    endif
                end do
                !$omp end parallel do
            endif
            !$omp parallel do default(shared) private(g) schedule(dynamic) proc_bind(close)
            do g = 1, k
                call filtered_refs(g)
            end do
            !$omp end parallel do
            if( present(kcut) .and. present(llband) )then
                call self%estep(ia, wa, refb, l_alive, wan, l_full=l_fb, ll=lla, kcut=kcut, llband=lba, lpp=lpa, l_flat=l_st)
                call self%estep(ib, wb, refa, l_alive, wbn, l_full=l_fb, ll=llb, kcut=kcut, llband=lbb, lpp=lpb, l_flat=l_st)
            else
                call self%estep(ia, wa, refb, l_alive, wan, l_full=l_fb, ll=lla, l_flat=l_st)
                call self%estep(ib, wb, refa, l_alive, wbn, l_full=l_fb, ll=llb, l_flat=l_st)
            endif
            if( l_st )then
                ! the site path classifies hard: soft weights at a noise-level site contrast pull every group map to the consensus
                call harden(wan)
                call harden(wbn)
            endif
            mpa  = maxloc(wa, dim=2)
            mpb  = maxloc(wb, dim=2)
            chga = real(count(maxloc(wan,dim=2) /= mpa)) / real(size(ia))
            chgb = real(count(maxloc(wbn,dim=2) /= mpb)) / real(size(ib))
            line = '  iteration '//int2str(it)//': MAP changed A '//fmtr(real(chga,dp),3)//' B '//fmtr(real(chgb,dp),3)
            wa  = wan
            wb  = wbn
            mpa = maxloc(wa, dim=2)
            mpb = maxloc(wb, dim=2)
            ca  = [(count(mpa == g), g=1,k)]
            cb  = [(count(mpb == g), g=1,k)]
            do g = 1, k
                if( .not. l_alive(g) ) cycle
                ! no size floor: a group leaves only when a half holds none of it; the re-test beyond the EM band decides
                if( min(ca(g), cb(g)) == 0 .and. count(l_alive) > 1 )then
                    l_alive(g) = .false.
                    write(logfhandle,'(A)') '  group '//int2str(g)//' drained ('//int2str(ca(g))//' / '//int2str(cb(g))//&
                        &' MAP particles) and dropped'
                endif
            end do
            line = line//'; weights A'
            do g = 1, k
                line = line//' '//fmtr(real(sum(wa(:,g)),dp)/real(size(ia),dp),3)
            end do
            line = line//'; deviation FSC per group (shell 4/8/12)'
            do g = 1, k
                if( .not. l_alive(g) ) cycle
                line = line//' '//fmtr(real(fl(min(4,self%nyq),g),dp),2)//'/'//fmtr(real(fl(min(8,self%nyq),g),dp),2)//&
                    &'/'//fmtr(real(fl(min(12,self%nyq),g),dp),2)
            end do
            write(logfhandle,'(A)') line
            ! the site path runs to convergence (SITE_CONV of the assignments changing), SITE_ITS_MIN to niter iterations
            if( l_st .and. it >= SITE_ITS_MIN .and. max(chga, chgb) < SITE_CONV ) exit
        end do
        la = maxloc(wa, dim=2, mask=spread(l_alive, 1, size(ia)))
        lb = maxloc(wb, dim=2, mask=spread(l_alive, 1, size(ib)))
        if( present(ll) ) ll = lla + llb
        if( present(kcut) .and. present(llband) ) llband = lba + lbb
        do g = 1, k
            call refa(g)%kill
            call refb(g)%kill
        end do
        deallocate(ma, mb, refa, refb, wa, wb, fl, l_alive)

    contains

        subroutine filtered_refs( gg )
            integer, intent(in) :: gg
            real :: f(0:self%nyq), fs(0:self%nyq)
            integer :: kk
            call self%shell_fsc(ma(:,:,:,gg) - ma(:,:,:,k+1), mb(:,:,:,gg) - mb(:,:,:,k+1), f)
            do kk = 0, self%nyq
                fs(kk) = (f(max(kk-1,0)) + f(kk) + f(min(kk+1,self%nyq))) / 3.
            end do
            fs = min(max(fs, 0.), 1.)
            fl(:,gg) = fs
            if( l_st ) fs = 1.
            call set_ref(refa(gg), ma(:,:,:,k+1), ma(:,:,:,gg), fs)
            call set_ref(refb(gg), mb(:,:,:,k+1), mb(:,:,:,gg), fs)
        end subroutine filtered_refs

        subroutine set_ref( ref, c0, cg, fs )
            type(reconstructor), intent(inout) :: ref
            complex,             intent(in)    :: c0(:,:,:), cg(:,:,:)
            real,                intent(in)    :: fs(0:)
            complex, allocatable :: r(:,:,:)
            integer :: h, kk, l
            allocate(r(self%nh,self%box,self%box))
            do l = 1, self%box
                do kk = 1, self%box
                    do h = 1, self%nh
                        r(h,kk,l) = c0(h,kk,l) + fs(self%fsh(h,kk,l)) * (cg(h,kk,l) - c0(h,kk,l))
                    end do
                end do
            end do
            call ref%set_ft(.true.)
            call ref%set_cmat(r)
            call ref%expand_exp
        end subroutine set_ref

        subroutine harden( w )
            real, intent(inout) :: w(:,:)
            integer :: i, g
            do i = 1, size(w,1)
                g = maxloc(w(i,:), dim=1)
                w(i,:) = 0.
                w(i,g) = 1.
            end do
        end subroutine harden

    end subroutine cross_em

    ! Cross-half EM and its acceptance driver use this same weight-to-map wrapper.
    !> l_mask: the maps masked by the focus mask (the site split's references), only on the direct insertion path
    subroutine em_maps( self, idx, w, eps, neq, l_full, m, audit_tag )
        class(state_mergeup), intent(inout) :: self
        integer,              intent(in) :: idx(:), neq
        real,                 intent(in) :: w(:,:), eps(0:)
        logical,              intent(in) :: l_full
        complex,              intent(inout) :: m(:,:,:,:)
        character(len=*), optional, intent(in) :: audit_tag
        integer, allocatable :: jidx(:), joff(:)
        real, allocatable :: jscale(:), weights(:)
        integer :: gg, nn, nk, coverage(3)
        nn = size(idx)
        nk = size(w,2)
        jscale = [[(max(sum(w(:,gg)), 1.)/real(neq), gg=1,nk)], real(nn)/real(neq)]
        coverage = 0
        if( self%p%merge_insert_reuse == 'yes' )then
            call self%maps_weighted_shared(idx, w, jscale, eps, .not. l_full, m, coverage)
        else
            jidx = [[(idx, gg=1,nk)], idx]
            joff = [(1 + gg*nn, gg=0,nk+1)]
            weights = [reshape(w, [nn*nk]), [(1., gg=1,nn)]]
            call self%maps(jidx, joff, jscale, eps, .false., .not. l_full, m, weights, coverage(1))
        endif
        if( present(audit_tag) ) call self%write_em_audit(audit_tag, idx, w, eps, jscale, l_full, m, coverage)
    end subroutine em_maps

    ! Diagnostic replay of the production weight-to-map wrapper; source project is read-only.
    subroutine audit_em_wrapper( self )
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use simple_syslib, only: dir_exists, simple_mkdir
        class(state_mergeup), intent(inout) :: self
        class(parameters), pointer :: params
        class(builder), pointer :: build
        integer, allocatable :: ia(:), ib(:), ga(:), gb(:)
        real, allocatable :: wa(:,:), wb(:,:), eps(:)
        complex, allocatable :: ma(:,:,:,:), mb(:,:,:,:)
        character(len=4096) :: line, paths(5)
        character(len=17), parameter :: lives(3) = ['fresh            ', 'reset_same_object', 'kill_recreate    ']
        character(len=:), allocatable :: tag
        integer :: fu, io, na, nb, nk, nthr, slots, nit, band, neq, j, life, it, old_threads, observed
        logical :: old_dynamic, old_contrast, old_amplitude
        if( self%p%merge_insert_audit /= 'yes' ) THROW_HARD('EM wrapper audit requires explicit opt-in')
        if( .not. self%l_exists ) THROW_HARD('EM wrapper audit needs initialized particles')
        params => self%p
        build => self%b
        old_contrast = self%l_con
        old_amplitude = self%l_amp
        call fopen(fu, params%merge_audit_request, status='old', action='read', iostat=io)
        call fileiochk('cannot read EM wrapper request', io)
        read(fu,'(A)',iostat=io) line
        if( io /= 0 ) THROW_HARD('invalid EM wrapper request schema')
        if( trim(line) /= 'simple-em-wrapper-v1' ) THROW_HARD('invalid EM wrapper request schema')
        read(fu,*,iostat=io) na, nb, nk, nthr, slots, nit, band, neq
        if( io /= 0 ) THROW_HARD('invalid EM wrapper controls')
        if( min(na,nb) < 0 .or. max(na,nb) > self%n .or. nk < 1 .or. nk > 8 ) THROW_HARD('invalid EM wrapper dimensions')
        if( nthr < 1 .or. nthr > 8 .or. slots /= 8 .or. self%nslots /= slots ) THROW_HARD('EM wrapper audit needs eight prepared slots and one to eight threads')
        if( nit /= 2 .or. band < 0 .or. band > 1 .or. neq < 1 ) THROW_HARD('invalid EM wrapper controls')
        do j = 1, size(paths)
            read(fu,'(A)',iostat=io) paths(j)
            if( io /= 0 .or. len_trim(paths(j)) == 0 ) THROW_HARD('missing EM wrapper input path')
        end do
        call fclose(fu)
        call read_rows(trim(paths(1)), na, 0, ia, ga)
        call read_rows(trim(paths(2)), nb, 1, ib, gb)
        call read_weights(trim(paths(3)), na, wa)
        call read_weights(trim(paths(4)), nb, wb)
        if( nlines(string(trim(paths(5)))) /= self%nyq+1 ) THROW_HARD('epsilon_shape_mismatch')
        allocate(eps(0:self%nyq))
        call fopen(fu, string(trim(paths(5))), status='old', action='read', iostat=io)
        call fileiochk('cannot read EM epsilon', io)
        do j = 0, self%nyq
            read(fu,*,iostat=io) eps(j)
            if( io /= 0 ) THROW_HARD('invalid EM epsilon')
        end do
        call fclose(fu)
        if( any(.not. ieee_is_finite(eps)) .or. any(eps < 0.) ) THROW_HARD('invalid EM epsilon')
        if( dir_exists(string('merge_insert_audit')) ) THROW_HARD('EM audit output must be fresh')
        call simple_mkdir(string('merge_insert_audit'))
        allocate(ma(self%nh,self%box,self%box,nk+1), mb(self%nh,self%box,self%box,nk+1))
        old_threads = 1
        old_dynamic = .false.
        !$ old_threads = omp_get_max_threads()
        !$ old_dynamic = omp_get_dynamic()
        !$ call omp_set_dynamic(.false.)
        !$ call omp_set_num_threads(nthr)
        observed = 1
        !$omp parallel default(shared)
        !$omp single
        !$ observed = omp_get_num_threads()
        !$omp end single
        !$omp end parallel
        if( observed /= nthr ) THROW_HARD('EM audit observed thread count differs from request')
        do life = 1, size(lives)
            if( life == 3 )then
                call self%kill()
                call self%new(params, build, EM_LIM, contrast=old_contrast, amplitude=old_amplitude)
                if( params%pcg_mskfile%strlen() > 0 ) call self%focus(params%pcg_mskfile)
                if( self%nslots /= slots ) THROW_HARD('EM audit slot count changed after recreation')
                if( any(self%rows(ia) /= ga) .or. any(self%rows(ib) /= gb) ) THROW_HARD('EM audit rows changed after recreation')
            endif
            do it = 1, nit
                tag = trim(lives(life))//'_iter_'//int2str(it)
                call self%em_maps(ia, wa, eps, neq, band == 0, ma, tag//'_a')
                call self%em_maps(ib, wb, eps, neq, band == 0, mb, tag//'_b')
            end do
        end do
        !$ call omp_set_num_threads(old_threads)
        !$ call omp_set_dynamic(old_dynamic)
        call fopen(fu, string('merge_insert_audit/complete.txt'), status='new', action='write', iostat=io)
        call fileiochk('cannot publish EM wrapper completion', io)
        write(fu,'(A)') 'simple-em-wrapper-v1'
        write(fu,'(A)') 'production_path state_mergeup/em_maps'
        write(fu,'(A,I0)') 'observed_threads ', observed
        write(fu,'(A)') 'iterations repeat frozen wrapper inputs; no E-step or docking in stress replay'
        call fclose(fu)

    contains

        subroutine read_rows( path, n, half, idx, global_rows )
            character(len=*), intent(in) :: path
            integer, intent(in) :: n, half
            integer, allocatable, intent(out) :: idx(:), global_rows(:)
            integer :: f, status, r
            if( nlines(string(path)) /= n ) THROW_HARD('row_shape_mismatch')
            allocate(idx(n), global_rows(n))
            call fopen(f, string(path), status='old', action='read', iostat=status)
            call fileiochk('cannot read EM row input', status)
            do r = 1, n
                read(f,*,iostat=status) global_rows(r)
                if( status /= 0 ) THROW_HARD('out_of_range_row')
                idx(r) = findloc(self%rows, global_rows(r), dim=1)
                if( idx(r) == 0 ) THROW_HARD('out_of_range_row')
                if( self%eo(idx(r)) /= half ) THROW_HARD('row_half_mismatch')
                if( any(idx(:r-1) == idx(r)) ) THROW_HARD('duplicate_row')
            end do
            call fclose(f)
        end subroutine read_rows

        subroutine read_weights( path, n, w )
            character(len=*), intent(in) :: path
            integer, intent(in) :: n
            real, allocatable, intent(out) :: w(:,:)
            character(len=4096) :: record
            real :: extra
            integer :: f, status, r
            if( nlines(string(path)) /= n ) THROW_HARD('weight_shape_mismatch')
            allocate(w(n,nk))
            call fopen(f, string(path), status='old', action='read', iostat=status)
            call fileiochk('cannot read EM weight input', status)
            do r = 1, n
                read(f,'(A)',iostat=status) record
                if( status /= 0 ) THROW_HARD('weight_shape_mismatch')
                read(record,*,iostat=status) w(r,:)
                if( status /= 0 ) THROW_HARD('weight_shape_mismatch')
                read(record,*,iostat=status) w(r,:), extra
                if( status >= 0 ) THROW_HARD('weight_shape_mismatch')
            end do
            call fclose(f)
            if( any(.not. ieee_is_finite(w)) .or. any(w < 0.) ) THROW_HARD('invalid_weights')
        end subroutine read_weights

    end subroutine audit_em_wrapper

    ! Native-order arrays after compression/reduction/finish and before docking; no project writes.
    subroutine write_em_audit( self, tag, idx, w, eps, jscale, l_full, m, coverage )
        use simple_syslib, only: dir_exists, simple_mkdir
        class(state_mergeup), intent(in) :: self
        character(len=*), intent(in) :: tag
        integer, intent(in) :: idx(:), coverage(3)
        real, intent(in) :: w(:,:), eps(0:), jscale(:)
        logical, intent(in) :: l_full
        complex, intent(in) :: m(:,:,:,:)
        character(len=:), allocatable :: prefix
        character(len=4) :: order
        integer :: fu, io, j, p, nthr
        order = transfer(1, order)
        if( storage_size(1) /= 32 .or. storage_size(1.) /= 32 .or. storage_size(cmplx(0.,0.)) /= 64 ) &
            THROW_HARD('EM audit requires int32, float32 and complex64 storage')
        if( iachar(order(1:1)) /= 1 ) THROW_HARD('EM audit requires little-endian storage')
        if( tag == 'cross_1_a' .and. dir_exists(string('merge_insert_audit')) ) THROW_HARD('EM audit output must be fresh')
        if( .not. dir_exists(string('merge_insert_audit')) ) call simple_mkdir(string('merge_insert_audit'))
        prefix = 'merge_insert_audit/'//tag
        nthr = 1
        !$ nthr = omp_get_max_threads()
        call fopen(fu, string(prefix//'.meta'), status='new', action='write', iostat=io)
        call fileiochk('cannot write EM audit metadata', io)
        write(fu,'(A)') 'simple-em-wrapper-v1'
        write(fu,'(A)') 'phase post_finish_pre_dock'
        write(fu,'(A)') 'storage little_endian Fortran_order int32 float32 complex64'
        write(fu,'(A,4(1X,I0))') 'shape', shape(m)
        write(fu,'(A,2(1X,I0))') 'weights_shape', shape(w)
        write(fu,'(A,2(1X,I0))') 'threads_slots', nthr, self%nslots
        write(fu,'(A,3(1X,I0))') 'chunks_tiles_batches', coverage
        write(fu,'(A,I0)') 'full_band ', merge(1,0,l_full)
        write(fu,'(A)') 'insert_reuse '//self%p%merge_insert_reuse
        write(fu,'(A)') 'pgrp '//trim(self%p%pgrp)
        call fclose(fu)
        call open_raw('.maps.c64')
        write(fu,iostat=io) m
        call close_raw()
        call open_raw('.weights.f32')
        write(fu,iostat=io) w
        call close_raw()
        call open_raw('.eps.f32')
        write(fu,iostat=io) eps
        call close_raw()
        call open_raw('.jscale.f32')
        write(fu,iostat=io) jscale
        call close_raw()
        call open_raw('.rows.i32')
        write(fu,iostat=io) self%rows(idx)
        call close_raw()
        call open_raw('.halves.i32')
        write(fu,iostat=io) self%eo(idx)
        call close_raw()
        call open_raw('.operators.bin')
        write(fu,iostat=io) self%p%box, self%box, self%p%smpd, self%p%smpd_crop, OSMPL_PAD_FAC, &
            &self%kb, self%emhk, self%eul(:,idx)
        call close_raw()
        call open_raw('.planes.bin')
        do j = 1, size(idx)
            p = idx(j)
            write(fu,iostat=io) lbound(self%fpls(p)%cmplx_plane), ubound(self%fpls(p)%cmplx_plane), &
                &self%fpls(p)%nyq, self%fpls(p)%frlims, self%fpls(p)%shconst, &
                &self%fpls(p)%cmplx_plane, self%fpls(p)%ctfsq_plane
            call fileiochk('cannot write EM input plane', io)
        end do
        io = 0
        call close_raw()

    contains

        subroutine open_raw( suffix )
            character(len=*), intent(in) :: suffix
            open(newunit=fu, file=prefix//suffix, status='new', action='write', access='stream', form='unformatted', iostat=io)
            call fileiochk('cannot create EM audit array', io)
        end subroutine open_raw

        subroutine close_raw()
            call fileiochk('cannot write EM audit array', io)
            close(fu,iostat=io)
            call fileiochk('cannot close EM audit array', io)
        end subroutine close_raw

    end subroutine write_em_audit

    subroutine site_em( self, seed, sel, mfile, post, nstar, curve, devs, score, l_all )
        use simple_srch_sort_loc, only: hpsort
        class(state_mergeup), intent(inout) :: self
        logical,              intent(in)    :: seed(:), sel(:)
        class(string),        intent(in)    :: mfile
        real,                 intent(in)    :: score(:)
        real,    allocatable, intent(out)   :: post(:), curve(:,:)
        integer,              intent(out)   :: nstar
        complex, allocatable, intent(out)   :: devs(:,:,:,:)
        logical,              intent(in)    :: l_all
        type(image)          :: vol
        complex, allocatable :: mp(:,:,:,:,:), rf(:,:,:,:,:), tm(:,:,:,:), da(:,:,:), db(:,:,:)
        real,    allocatable :: sm(:,:,:), wa(:,:), wb(:,:), rot(:,:,:,:), fw2(:), key(:), cn(:)
        real(dp),allocatable :: epsb(:), epsf(:)
        complex(dp), allocatable :: nma(:,:,:), nmb(:,:,:), nm2(:,:,:), nmp(:,:,:)
        real(dp),    allocatable :: dna(:,:,:), dnb(:,:,:), dn2(:,:,:), dnp(:,:,:)
        integer, allocatable :: ia(:), ib(:), oa(:), ob(:), fhk(:,:)
        real(dp) :: bta(self%nshem), btb(self%nshem), s2a(self%nshem), s2b(self%nshem), pra(2), prb(2)
        real    :: f(0:self%nyq), cbest, chga, chgb
        integer :: it, g, ig, nh2, i, ngrid, nfull, h, k, r2, nsym, npa, npb
        call vol%read_and_crop(mfile, self%p%smpd, self%box, self%p%smpd_crop)
        sm = min(max(vol%get_rmat(), 0.), 1.) * self%msk3
        call vol%kill
        rot  = self%sym_rots()
        nsym = size(rot,3)
        if( allocated(self%con) )then
            cn = self%con
        else
            allocate(cn(self%n), source=1.)
        endif
        nfull = 0
        do ig = 1, 2
            do k = lbound(self%fpls(1)%cmplx_plane,2), 0
                do h = lbound(self%fpls(1)%cmplx_plane,1), ubound(self%fpls(1)%cmplx_plane,1)
                    r2 = h*h + k*k
                    if( r2 < 3 .or. r2 > self%khi*(self%khi+1) ) cycle
                    nfull = nfull + 1
                    if( ig == 1 ) cycle
                    fhk(:,nfull) = [h, k]
                    fw2(nfull)   = merge(1., 2., k == 0)
                end do
            end do
            if( ig == 1 ) allocate(fhk(2,nfull), fw2(nfull))
            if( ig == 1 ) nfull = 0
        end do
        ia = pack([(i, i=1,self%n)], self%eo == 0 .and. sel)
        ib = pack([(i, i=1,self%n)], self%eo == 1 .and. sel)
        allocate(wa(size(ia),2), wb(size(ib),2))
        if( size(score) > 0 )then
            wa(:,2) = score(ia)
            wb(:,2) = score(ib)
        else
            wa(:,2) = merge(1., 0., seed(ia))
            wb(:,2) = merge(1., 0., seed(ib))
        endif
        wa(:,1) = 1. - wa(:,2)
        wb(:,1) = 1. - wb(:,2)
        epsb = 0.5_dp * (density(ia, self%emhk, self%emw) + density(ib, self%emhk, self%emw))
        epsf = 0.5_dp * (density(ia, fhk, fw2) + density(ib, fhk, fw2))
        allocate(mp(self%nh,self%box,self%box,3,2), rf(self%nh,self%box,self%box,2,2))
        if( size(score) > 0 )then
            write(logfhandle,'(A)') 'site ranking by a given score (no EM), the deviation inside the site mask'
        else
            write(logfhandle,'(A)') 'site EM: 2 classes from a '//int2str(count(seed))//'-particle start, until under '//&
                &fmtr(100._dp*real(SITE_CONV,dp),1)//' % of the assignments change ('//int2str(SITE_ITS_MIN)//' to '//&
                &int2str(SITE_ITS)//' iterations), shells 2..'//int2str(self%kb)//', nearest-neighbour maps, the deviation '//&
                &'inside the site mask ('//fmtr(100._dp*real(count(sm > 0.5),dp)/real(max(1,count(self%msk3 > 0.5)),dp),1)//&
                &' % of the sphere)'
        endif
        ! the all-particle sums are fixed and the maps linear in the weights: class 1 = all - class 2
        call nn_acc(ia, [(1., i=1,size(ia))], self%emhk, self%emw, nma, dna)
        call nn_acc(ib, [(1., i=1,size(ib))], self%emhk, self%emw, nmb, dnb)
        call nn_norm(nma, dna, epsb, real(size(ia),dp), mp(:,:,:,3,1))
        call nn_norm(nmb, dnb, epsb, real(size(ib),dp), mp(:,:,:,3,2))
        do it = 1, merge(0, SITE_ITS, size(score) > 0)
            call nn_acc(ia, wa(:,2), self%emhk, self%emw, nm2, dn2)
            call nn_norm(nma - nm2, dna - dn2, epsb, sum(real(wa(:,1),dp)), mp(:,:,:,1,1))
            call nn_norm(nm2, dn2, epsb, sum(real(wa(:,2),dp)), mp(:,:,:,2,1))
            call nn_acc(ib, wb(:,2), self%emhk, self%emw, nm2, dn2)
            call nn_norm(nmb - nm2, dnb - dn2, epsb, sum(real(wb(:,1),dp)), mp(:,:,:,1,2))
            call nn_norm(nm2, dn2, epsb, sum(real(wb(:,2),dp)), mp(:,:,:,2,2))
            do g = 1, 2
                call self%shell_fsc(mp(:,:,:,g,1) - mp(:,:,:,3,1), mp(:,:,:,g,2) - mp(:,:,:,3,2), f)
                f = min(max(f, 0.), 1.)
                rf(:,:,:,g,1) = mp(:,:,:,3,1) + masked(shell_weighted(mp(:,:,:,g,1) - mp(:,:,:,3,1), f))
                rf(:,:,:,g,2) = mp(:,:,:,3,2) + masked(shell_weighted(mp(:,:,:,g,2) - mp(:,:,:,3,2), f))
            end do
            call nn_estep(ia, wa, rf(:,:,:,:,2), it, chga, bta, s2a, pra)
            call nn_estep(ib, wb, rf(:,:,:,:,1), it, chgb, btb, s2b, prb)
            ! the assignments stall for a few iterations before a growing class takes off
            if( it >= SITE_ITS_MIN .and. max(chga, chgb) < SITE_CONV ) exit
        end do
        if( size(score) == 0 ) write(logfhandle,'(A)') '  site EM '//&
            &trim(merge('converged after   ', 'stopped at the cap', it <= SITE_ITS))//' '//int2str(min(it, SITE_ITS))//' iterations'
        allocate(post(self%n), source=0.)
        post(ia) = wa(:,2)
        post(ib) = wb(:,2)
        if( l_all .and. size(score) == 0 )then
            call nn_score(pack([(i, i=1,self%n)], self%eo == 0 .and. .not. sel), rf(:,:,:,:,2), bta, s2a, pra)
            call nn_score(pack([(i, i=1,self%n)], self%eo == 1 .and. .not. sel), rf(:,:,:,:,1), btb, s2b, prb)
            write(logfhandle,'(A)') '  scored '//int2str(count(.not. sel))//' particles outside the '//int2str(count(sel))//&
                &'-particle fit with the converged maps'
        endif
        deallocate(mp, rf)
        oa  = [(i, i=1,size(ia))]
        key = -post(ia)
        call hpsort(key, oa)
        ob  = [(i, i=1,size(ib))]
        key = -post(ib)
        call hpsort(key, ob)
        allocate(tm(self%nh,self%box,self%box,4))
        allocate(devs(self%nh,self%box,self%box,2), source=cmplx(0.,0.))
        call nn_map(ia, [(1., i=1,size(ia))], fhk, fw2, epsf, tm(:,:,:,1))
        call nn_map(ib, [(1., i=1,size(ib))], fhk, fw2, epsf, tm(:,:,:,2))
        ngrid = max(1, count(SITE_GRID <= int(0.7 * real(count(sel)))))
        allocate(curve(4,ngrid), source=0.)
        cbest = -huge(cbest)
        nstar = SITE_GRID(1)
        npa   = 0
        npb   = 0
        do ig = 1, ngrid
            nh2 = SITE_GRID(ig) / 2
            ! the grid's top-N sets are nested: each one adds its new members to the previous sums
            call nn_grow(ia, oa, min(nh2,size(ia)), npa, nma, dna)
            call nn_norm(nma, dna, epsf, real(npa,dp), tm(:,:,:,3))
            call nn_grow(ib, ob, min(nh2,size(ib)), npb, nmb, dnb)
            call nn_norm(nmb, dnb, epsf, real(npb,dp), tm(:,:,:,4))
            da = masked(tm(:,:,:,3) - tm(:,:,:,1))
            db = masked(tm(:,:,:,4) - tm(:,:,:,2))
            curve(1,ig) = real(SITE_GRID(ig))
            curve(2,ig) = self%shell_corr(da, db, self%kb + RB_MARGIN, self%khi)
            curve(3,ig) = real(sqrt(max(self%shell_dot(da, db, self%kb + RB_MARGIN, self%khi), 0._dp)))
            write(logfhandle,'(A)') '  N '//int2str(SITE_GRID(ig))//': hold-out cross-half correlation of the site deviation '//&
                &fmtr(real(curve(2,ig),dp),3)
            if( curve(2,ig) > cbest )then
                cbest = curve(2,ig)
                nstar = SITE_GRID(ig)
                devs(:,:,:,1) = da
                devs(:,:,:,2) = db
            endif
        end do
        if( curve(3,1) > 0. )then
            curve(4,:) = curve(3,:) * curve(1,:) / curve(3,1)
        else
            curve(4,:) = real(nstar)
        endif
        write(logfhandle,'(A)') '  N* '//int2str(nstar)//' (hold-out shells '//int2str(self%kb + RB_MARGIN)//'..'//&
            &int2str(self%khi)//', correlation '//fmtr(real(cbest,dp),3)//'), equivalent pure members '//&
            &int2str(nint(maxval(curve(4,:))))
        deallocate(tm, rot)

    contains

        function density( pidx, hk, sw ) result( epsv )
            integer, intent(in)   :: pidx(:), hk(:,:)
            real,    intent(in)   :: sw(:)
            real(dp), allocatable :: epsv(:)
            real(dp), allocatable :: den(:,:,:)
            real,     allocatable :: vals(:)
            real(dp) :: c2v
            integer  :: ii, j, is, c1(3), c2(3), s, cnt
            logical  :: cnj
            allocate(den(self%nh,self%box,self%box), source=0._dp)
            !$omp parallel do default(shared) private(ii,j,is,c1,c2,cnj,c2v) schedule(static) reduction(+:den) proc_bind(close)
            do ii = 1, size(pidx)
                do j = 1, size(sw)
                    c2v = real(sw(j) * cn(pidx(ii))**2, dp) * real(self%fpls(pidx(ii))%ctfsq_plane(hk(1,j),hk(2,j)), dp)
                    do is = 1, nsym
                        call plane_cell(rot(:,:,is,pidx(ii)), hk(1,j), hk(2,j), self%nyq, self%box, c1, c2, cnj)
                        if( c1(1) == 0 ) cycle
                        if( c2(1) > 0 )then
                            den(c1(1),c1(2),c1(3)) = den(c1(1),c1(2),c1(3)) + 0.5_dp * c2v
                            den(c2(1),c2(2),c2(3)) = den(c2(1),c2(2),c2(3)) + 0.5_dp * c2v
                        else
                            den(c1(1),c1(2),c1(3)) = den(c1(1),c1(2),c1(3)) + c2v
                        endif
                    end do
                end do
            end do
            !$omp end parallel do
            allocate(epsv(0:self%nyq), source=1._dp)
            do s = 0, self%nyq
                vals = real(pack(den, self%fsh == s .and. den > 0._dp))
                cnt  = size(vals)
                if( cnt == 0 ) cycle
                call hpsort(vals)
                epsv(s) = real(vals((cnt + 1) / 2), dp) / real(size(pidx), dp)
            end do
        end function density

        subroutine nn_map( pidx, pw, hk, sw, epsv, out )
            integer,  intent(in)  :: pidx(:), hk(:,:)
            real,     intent(in)  :: pw(:), sw(:)
            real(dp), intent(in)  :: epsv(0:)
            complex,  intent(out) :: out(:,:,:)
            complex(dp), allocatable :: num(:,:,:)
            real(dp),    allocatable :: den(:,:,:)
            call nn_acc(pidx, pw, hk, sw, num, den)
            call nn_norm(num, den, epsv, sum(real(pw,dp)), out)
        end subroutine nn_map

        subroutine nn_grow( pidx, ord, k, np, num, den )
            integer,                  intent(in)    :: pidx(:), ord(:), k
            integer,                  intent(inout) :: np
            complex(dp), allocatable, intent(inout) :: num(:,:,:)
            real(dp),    allocatable, intent(inout) :: den(:,:,:)
            integer :: jj
            if( np == 0 )then
                call nn_acc(pidx(ord(1:k)), [(1., jj=1,k)], fhk, fw2, num, den)
            else if( k > np )then
                call nn_acc(pidx(ord(np+1:k)), [(1., jj=np+1,k)], fhk, fw2, nm2, dn2)
                num = num + nm2
                den = den + dn2
            endif
            np = max(np, k)
        end subroutine nn_grow

        subroutine nn_norm( num, den, epsv, wpw, out )
            complex(dp), intent(in)  :: num(:,:,:)
            real(dp),    intent(in)  :: den(:,:,:), epsv(0:), wpw
            complex,     intent(out) :: out(:,:,:)
            real(dp) :: wsum
            integer  :: hh, kk, ll
            wsum = max(wpw, 1._dp)
            do ll = 1, self%box
                do kk = 1, self%box
                    do hh = 1, self%nh
                        out(hh,kk,ll) = cmplx(num(hh,kk,ll) / (den(hh,kk,ll) + epsv(self%fsh(hh,kk,ll)) * wsum), kind=kind(out))
                    end do
                end do
            end do
        end subroutine nn_norm

        subroutine nn_acc( pidx, pw, hk, sw, num, den )
            integer,                  intent(in)  :: pidx(:), hk(:,:)
            real,                     intent(in)  :: pw(:), sw(:)
            complex(dp), allocatable, intent(out) :: num(:,:,:)
            real(dp),    allocatable, intent(out) :: den(:,:,:)
            complex(dp) :: v, vj
            real(dp)    :: c2v
            integer     :: ii, j, is, c1(3), c2(3)
            logical     :: cnj
            complex(dp), allocatable :: numl(:,:,:)
            real(dp),    allocatable :: denl(:,:,:)
            allocate(numl(self%nh,self%box,self%box), source=cmplx(0._dp,0._dp,dp))
            allocate(denl(self%nh,self%box,self%box), source=0._dp)
            !$omp parallel do default(shared) private(ii,j,is,c1,c2,cnj,v,vj,c2v) schedule(static) reduction(+:numl,denl) &
            !$omp& proc_bind(close)
            do ii = 1, size(pidx)
                if( pw(ii) <= 0. ) cycle
                do j = 1, size(sw)
                    vj  = real(pw(ii) * sw(j) * cn(pidx(ii)), dp) * cmplx(self%fpls(pidx(ii))%cmplx_plane(hk(1,j),hk(2,j)), kind=dp)
                    c2v = real(pw(ii) * sw(j) * cn(pidx(ii))**2, dp) * real(self%fpls(pidx(ii))%ctfsq_plane(hk(1,j),hk(2,j)), dp)
                    do is = 1, nsym
                        call plane_cell(rot(:,:,is,pidx(ii)), hk(1,j), hk(2,j), self%nyq, self%box, c1, c2, cnj)
                        if( c1(1) == 0 ) cycle
                        v = vj
                        if( cnj ) v = conjg(v)
                        if( c2(1) > 0 )then
                            numl(c1(1),c1(2),c1(3)) = numl(c1(1),c1(2),c1(3)) + 0.5_dp * v
                            denl(c1(1),c1(2),c1(3)) = denl(c1(1),c1(2),c1(3)) + 0.5_dp * c2v
                            numl(c2(1),c2(2),c2(3)) = numl(c2(1),c2(2),c2(3)) + 0.5_dp * conjg(v)
                            denl(c2(1),c2(2),c2(3)) = denl(c2(1),c2(2),c2(3)) + 0.5_dp * c2v
                        else
                            numl(c1(1),c1(2),c1(3)) = numl(c1(1),c1(2),c1(3)) + v
                            denl(c1(1),c1(2),c1(3)) = denl(c1(1),c1(2),c1(3)) + c2v
                        endif
                    end do
                end do
            end do
            !$omp end parallel do
            call move_alloc(numl, num)
            call move_alloc(denl, den)
        end subroutine nn_acc

        subroutine nn_estep( pidx, w, refs, iter, chg, beta, sig2, pr )
            integer,  intent(in)    :: pidx(:), iter
            real,     intent(inout) :: w(:,:)
            complex,  intent(in)    :: refs(:,:,:,:)
            real,     intent(out)   :: chg
            real(dp), intent(out)   :: beta(:), sig2(:), pr(2)
            real, allocatable :: aa(:,:,:), bb(:,:,:), wn(:,:)
            integer, allocatable :: mpc(:)
            real(dp) :: sa(self%nshem), sb(self%nshem), sf(self%nshem), sc(self%nshem)
            integer  :: ii, gg
            call nn_terms(pidx, refs, aa, bb)
            mpc = maxloc(w, dim=2)
            do gg = 1, 2
                pr(gg) = log(max(real(sum(w(:,gg)),dp) / real(size(pidx),dp), real(SITE_PFLOOR,dp)))
            end do
            sa = 0._dp
            sb = 0._dp
            sf = 0._dp
            sc = 0._dp
            !$omp parallel do default(shared) private(ii) schedule(static) reduction(+:sa,sb,sf,sc) proc_bind(close)
            do ii = 1, size(pidx)
                sa = sa + real(aa(:,mpc(ii),ii),dp)
                sb = sb + real(bb(:,mpc(ii),ii),dp)
                sf = sf + real(self%f2s(:,pidx(ii)),dp)
                sc = sc + real(self%cw,dp)
            end do
            !$omp end parallel do
            beta = min(max(sa / max(sb, 1.e-30_dp), 0._dp), 2._dp)
            sig2 = max((sf - 2._dp*beta*sa + beta*beta*sb) / max(sc, 1._dp), 1.e-12_dp)
            wn   = nn_post(aa, bb, beta, sig2, pr)
            chg  = real(count((wn(:,2) > 0.5) .neqv. (w(:,2) > 0.5))) / real(max(1, size(pidx)))
            write(logfhandle,'(A)') '  iteration '//int2str(iter)//': argmax changed '//fmtr(real(chg,dp),3)//&
                &', site class weight '//fmtr(real(sum(wn(:,2)),dp)/real(size(pidx),dp),3)//&
                &', model scale at shells 2 / '//int2str(self%kb)//' '//fmtr(beta(1),2)//' / '//fmtr(beta(self%nshem),2)
            w = wn
        end subroutine nn_estep

        subroutine nn_score( pidx, refs, beta, sig2, pr )
            integer,  intent(in) :: pidx(:)
            complex,  intent(in) :: refs(:,:,:,:)
            real(dp), intent(in) :: beta(:), sig2(:), pr(2)
            real, allocatable :: aa(:,:,:), bb(:,:,:), wn(:,:)
            if( size(pidx) == 0 ) return
            call nn_terms(pidx, refs, aa, bb)
            wn         = nn_post(aa, bb, beta, sig2, pr)
            post(pidx) = wn(:,2)
        end subroutine nn_score

        subroutine nn_terms( pidx, refs, aa, bb )
            integer,           intent(in)  :: pidx(:)
            complex,           intent(in)  :: refs(:,:,:,:)
            real, allocatable, intent(out) :: aa(:,:,:), bb(:,:,:)
            complex  :: pv, dv
            real     :: wc, wc2, c2v
            integer  :: ii, gg, j, is, c1(3), c2(3), sh
            logical  :: cnj
            allocate(aa(self%nshem,2,size(pidx)), bb(self%nshem,2,size(pidx)), source=0.)
            !$omp parallel do default(shared) private(ii,gg,j,is,c1,c2,cnj,sh,pv,dv,wc,wc2,c2v) schedule(static) proc_bind(close)
            do ii = 1, size(pidx)
                wc  = cn(pidx(ii)) / real(nsym)
                wc2 = cn(pidx(ii))**2 / real(nsym)
                do j = 1, self%nem
                    sh  = self%emsh(j)
                    dv  = self%fpls(pidx(ii))%cmplx_plane(self%emhk(1,j),self%emhk(2,j))
                    c2v = self%fpls(pidx(ii))%ctfsq_plane(self%emhk(1,j),self%emhk(2,j))
                    do is = 1, nsym
                        call plane_cell(rot(:,:,is,pidx(ii)), self%emhk(1,j), self%emhk(2,j), self%nyq, self%box, c1, c2, cnj)
                        if( c1(1) == 0 ) cycle
                        do gg = 1, 2
                            pv = refs(c1(1),c1(2),c1(3),gg)
                            if( cnj ) pv = conjg(pv)
                            aa(sh,gg,ii) = aa(sh,gg,ii) + self%emw(j) * wc * (real(pv)*real(dv) + aimag(pv)*aimag(dv))
                            bb(sh,gg,ii) = bb(sh,gg,ii) + self%emw(j) * wc2 * c2v * abs(pv)**2
                        end do
                    end do
                end do
            end do
            !$omp end parallel do
        end subroutine nn_terms

        function nn_post( aa, bb, beta, sig2, pr ) result( wn )
            real,     intent(in) :: aa(:,:,:), bb(:,:,:)
            real(dp), intent(in) :: beta(:), sig2(:), pr(2)
            real, allocatable :: wn(:,:)
            real(dp) :: s(2), smax
            integer  :: ii, gg
            allocate(wn(size(aa,3),2))
            !$omp parallel do default(shared) private(ii,gg,s,smax) schedule(static) proc_bind(close)
            do ii = 1, size(aa,3)
                do gg = 1, 2
                    s(gg) = 0.5_dp * sum((2._dp*beta*real(aa(:,gg,ii),dp) - beta*beta*real(bb(:,gg,ii),dp)) / sig2) + pr(gg)
                end do
                smax = maxval(s)
                wn(ii,:) = real(exp(s - smax))
                wn(ii,:) = wn(ii,:) / sum(wn(ii,:))
            end do
            !$omp end parallel do
        end function nn_post

        function shell_weighted( cm, fs ) result( cmo )
            complex, intent(in)  :: cm(:,:,:)
            real,    intent(in)  :: fs(0:)
            complex, allocatable :: cmo(:,:,:)
            integer :: hh, kk, ll
            allocate(cmo(self%nh,self%box,self%box))
            do ll = 1, self%box
                do kk = 1, self%box
                    do hh = 1, self%nh
                        cmo(hh,kk,ll) = fs(self%fsh(hh,kk,ll)) * cm(hh,kk,ll)
                    end do
                end do
            end do
        end function shell_weighted

        function masked( cm ) result( cmo )
            complex, intent(in)  :: cm(:,:,:)
            complex, allocatable :: cmo(:,:,:)
            real, pointer :: rmat(:,:,:)
            call self%wimgs(1)%set_ft(.true.)
            call self%wimgs(1)%set_cmat(cm)
            call self%wimgs(1)%ifft
            call self%wimgs(1)%get_rmat_ptr(rmat)
            rmat = rmat * sm
            call self%wimgs(1)%fft
            cmo = self%wimgs(1)%get_cmat()
        end function masked

    end subroutine site_em

    ! per-particle contrast against the other half's mean map at its own noise level: zero-contrast atom vs log-normal genuine
    subroutine fit_contrast( self )
        class(state_mergeup), intent(inout) :: self
        type(ivec)            :: hx(2)
        complex,  allocatable :: mp(:,:,:,:)
        real,     allocatable :: rot(:,:,:,:), a0(:,:), b0(:,:), x(:), v(:), snr(:)
        real(dp), allocatable :: epsv(:)
        logical,  allocatable :: ok(:)
        real(dp) :: sa(self%kj), sb(self%kj), sr(self%kj), sc(self%kj), beta(self%kj), s2(self%kj), a(self%kj,1), b(self%kj,1)
        real(dp) :: tq(NQ_CON), lgq(NQ_CON), den, r, et, et2, ea
        real     :: cw(self%kj)
        integer  :: e, i, ii, j, its, funit, ios
        rot = self%sym_rots()
        do e = 1, 2
            hx(e)%v = pack([(i, i=1,self%n)], self%eo == e - 1)
        end do
        epsv = 0.5_dp * (self%band_density(rot, hx(1)%v) + self%band_density(rot, hx(2)%v))
        allocate(mp(self%nh,self%box,self%box,2))
        do e = 1, 2
            call self%band_wiener(rot, hx(e)%v, epsv, mp(:,:,:,e))
        end do
        allocate(a0(self%kj,self%n), b0(self%kj,self%n))
        !$omp parallel do default(shared) private(i,a,b) schedule(static) proc_bind(close)
        do i = 1, self%n
            call self%band_terms(rot, i, mp(:,:,:,2-self%eo(i):2-self%eo(i)), a, b)
            a0(:,i) = real(a(:,1))
            b0(:,i) = real(b(:,1))
        end do
        !$omp end parallel do
        deallocate(mp, rot)
        cw = 0.
        do j = 1, self%nj
            cw(self%jsh(j)) = cw(self%jsh(j)) + self%jw(j)
        end do
        if( allocated(self%con) ) deallocate(self%con, self%cam, self%cse, self%cpj)
        if( allocated(self%cb)  ) deallocate(self%cb, self%cs)
        allocate(self%cb(self%kj,2), self%cs(self%kj,2), source=0._dp)
        allocate(self%con(self%n), source=1.)
        allocate(self%cam(self%n), self%cse(self%n), self%cpj(self%n), source=0.)
        allocate(ok(self%n), source=.false.)
        do e = 1, 2
            sa = 0._dp
            sb = 0._dp
            sr = 0._dp
            sc = 0._dp
            do ii = 1, size(hx(e)%v)
                i  = hx(e)%v(ii)
                sa = sa + real(a0(:,i),dp) / real(self%nv(i),dp)
                sb = sb + real(b0(:,i),dp) / real(self%nv(i),dp)
            end do
            beta = min(max(sa / max(sb, 1.e-30_dp), 0._dp), 2._dp)
            do ii = 1, size(hx(e)%v)
                i  = hx(e)%v(ii)
                sr = sr + (real(self%f2j(:,i),dp) - 2._dp*beta*a0(:,i) + beta*beta*b0(:,i)) / real(self%nv(i),dp)
                sc = sc + real(cw,dp)
            end do
            s2 = max(sr / max(sc, 1._dp), 1.e-12_dp)
            self%cb(:,e) = beta
            self%cs(:,e) = s2
            do ii = 1, size(hx(e)%v)
                i   = hx(e)%v(ii)
                den = sum(beta * beta * b0(:,i) / s2)
                if( den <= 0._dp ) cycle
                self%cam(i) = real(sum(beta * a0(:,i) / s2) / den)
                self%cse(i) = real(sqrt(real(self%nv(i),dp) / den))
                ok(i)       = .true.
            end do
        end do
        x = pack(self%cam, ok)
        v = pack(self%cse**2, ok)
        self%cfit = 0._dp
        self%nfit = 0
        if( size(x) > 1 )then
            do j = 1, NQ_CON
                tq(j) = -6._dp + 12._dp * real(j-1,dp) / real(NQ_CON-1,dp)
            end do
            lgq = -0.5_dp * tq * tq
            lgq = lgq - log(sum(exp(lgq)))
            call zl_fit
            do i = 1, self%n
                if( .not. ok(i) ) cycle
                call zl_post(real(self%cam(i),dp), real(self%cse(i),dp)**2, self%cfit, tq, lgq, r, et, et2, ea)
                self%cpj(i) = real(r)
                self%con(i) = real((1._dp - r) * ea)
            end do
        endif
        snr = pack(self%cam / max(self%cse, 1.e-12), ok)
        write(logfhandle,'(A)') '>>> STATE_MERGEUP contrast of '//int2str(count(ok))//' particles against the other half''s '//&
            &'mean map (shells 2..'//int2str(self%kj)//'), each at its own noise level: zero-contrast weight '//&
            &fmtr(self%cfit(1),4)//', genuine contrast median '//fmtr(exp(self%cfit(2)),3)//' and log-sd '//&
            &fmtr(self%cfit(3),3)//' ('//int2str(self%nfit)//' EM iterations); contrast / se 5 / 50 / 95 % '//&
            &fmtr(real(quantr(snr,0.05),dp),1)//' / '//fmtr(real(quantr(snr,0.5),dp),1)//' / '//&
            &fmtr(real(quantr(snr,0.95),dp),1)//'; '//int2str(count(self%cpj >= 0.5))//' at zero contrast (posterior >= 0.5)'
        call fopen(funit, string('site_contrast.txt'), status='replace', action='write', iostat=ios)
        call fileiochk('fit_contrast: cannot write site_contrast.txt', ios)
        write(funit,'(A)') '# project row, posterior-mean contrast'
        do i = 1, self%n
            write(funit,'(I0,1X,ES12.5)') self%rows(i), max(self%con(i), 1.e-6)
        end do
        call fclose(funit)

    contains

        subroutine zl_fit
            real(dp) :: p(3), pn(3), sr1, s1, st, st2, med
            integer  :: ii
            med  = real(medianr(pack(x, x > 0.)), dp)
            p(1) = 0.5_dp
            p(2) = log(max(med, 1.e-6_dp))
            p(3) = max(real(MAD_SD,dp) * real(medianr(abs(x - real(med))),dp) / max(med, 1.e-6_dp), 0.05_dp)
            do its = 1, JUNK_ITS
                sr1 = 0._dp
                s1  = 0._dp
                st  = 0._dp
                st2 = 0._dp
                !$omp parallel do default(shared) private(ii,r,et,et2,ea) schedule(static) reduction(+:sr1,s1,st,st2) &
                !$omp& proc_bind(close)
                do ii = 1, size(x)
                    call zl_post(real(x(ii),dp), real(v(ii),dp), p, tq, lgq, r, et, et2, ea)
                    sr1 = sr1 + r
                    s1  = s1  + (1._dp - r)
                    st  = st  + (1._dp - r) * et
                    st2 = st2 + (1._dp - r) * et2
                end do
                !$omp end parallel do
                pn(1) = sr1 / real(size(x),dp)
                pn(2) = p(2) + p(3) * st / max(s1, 1.e-12_dp)
                pn(3) = p(3) * sqrt(max(st2 / max(s1, 1.e-12_dp) - (st / max(s1, 1.e-12_dp))**2, 1.e-12_dp))
                self%nfit = its
                if( maxval(abs(pn - p)) < 1.e-6_dp )then
                    p = pn
                    exit
                endif
                p = pn
            end do
            self%cfit = p
        end subroutine zl_fit

    end subroutine fit_contrast

    ! per-node junk screen at the seeding labels: amplitudes against every state's other-half reference, 3-component mixture cut
    subroutine junk_screen( self, seeds )
        use simple_srch_sort_loc, only: hpsort
        class(state_mergeup), intent(inout) :: self
        integer,              intent(inout) :: seeds(:)
        type(ivec)            :: hx(2)
        complex,  allocatable :: mp(:,:,:,:,:), rf(:,:,:,:)
        real,     allocatable :: rot(:,:,:,:), fs(:,:), amp(:), var(:), post(:), x(:), v(:), key(:), res(:), rz(:), wg(:)
        real(dp), allocatable :: epsv(:)
        integer,  allocatable :: lab(:), best(:), cnt(:,:), ord(:), idx(:), asg(:)
        logical,  allocatable :: ljunk(:), lj(:), lres(:), lamp(:)
        character(len=:), allocatable :: line
        real(dp) :: mu(3), s2(3), pw(3), gap(3), xm, tc(3), lgc(3)
        real     :: f(0:self%nyq), mv, med, sdr
        integer  :: nk, i, g, e, c, cg, h, k, l, funit, rep, ios, nf
        nk = maxval([0, seeds])
        if( nk == 0 )then
            write(logfhandle,'(A)') '>>> JUNK_SCREEN no seed state, nothing screened'
            return
        endif
        if( .not. allocated(self%cpj) ) call self%fit_contrast
        lab = seeds
        call fopen(rep, string('site_junk_report.txt'), status='replace', action='write', iostat=ios)
        call fileiochk('junk_screen: cannot write site_junk_report.txt', ios)
        rot = self%sym_rots()
        allocate(cnt(nk,2))
        do e = 1, 2
            hx(e)%v  = pack([(i, i=1,self%n)], self%eo == e - 1)
            cnt(:,e) = [(count(lab(hx(e)%v) == g), g=1,nk)]
        end do
        call say(int2str(nk)//' seed states, '//int2str(self%n)//' particles, band shells 2..'//int2str(self%kj)//' ('//&
            &fmtr(real(real(self%box)*self%p%smpd_crop/real(self%kj),dp),1)//' A, '//int2str(self%nj)//' plane samples), '//&
            &'every particle against the other half''s state references')
        epsv = 0.5_dp * (self%band_density(rot, hx(1)%v) + self%band_density(rot, hx(2)%v))
        allocate(mp(self%nh,self%box,self%box,0:nk,2))
        do e = 1, 2
            call self%band_wiener(rot, hx(e)%v, epsv, mp(:,:,:,0,e))
            do g = 1, nk
                call self%band_wiener(rot, pack(hx(e)%v, lab(hx(e)%v) == g), epsv, mp(:,:,:,g,e))
            end do
        end do
        allocate(fs(0:self%nyq,nk))
        do g = 1, nk
            call self%shell_fsc(mp(:,:,:,g,1) - mp(:,:,:,0,1), mp(:,:,:,g,2) - mp(:,:,:,0,2), f)
            fs(:,g) = min(max(f, 0.), 1.)
        end do
        allocate(rf(self%nh,self%box,self%box,nk))
        allocate(amp(self%n), var(self%n), post(self%n), res(self%n), rz(self%n), source=0.)
        allocate(best(self%n), source=0)
        do e = 1, 2
            do g = 1, nk
                do l = 1, self%box
                    do k = 1, self%box
                        do h = 1, self%nh
                            rf(h,k,l,g) = mp(h,k,l,0,3-e) + fs(self%fsh(h,k,l),g) * (mp(h,k,l,g,3-e) - mp(h,k,l,0,3-e))
                        end do
                    end do
                end do
            end do
            call score(hx(e)%v, cnt(:,3-e) > 0, e)
        end do
        deallocate(mp, rf)
        x = pack(amp, best > 0)
        v = pack(var, best > 0)
        allocate(ljunk(self%n), lamp(self%n), source=.false.)
        allocate(asg(self%n), source=0)
        allocate(lj(3), source=.false.)
        if( size(x) > 0 )then
            xm = sum(real(x,dp)) / real(size(x),dp)
            mu = real([quantr(x, 0.02), quantr(x, 0.5), quantr(x, 0.5)], dp)
            s2 = [0.1_dp, 0.25_dp, 1._dp] * sum((real(x,dp) - xm)**2) / real(size(x),dp)
            pw = [0.05_dp, 0.75_dp, 0.2_dp]
            call xd_fit
            mv = medianr(v)
            cg = maxloc(pw, dim=1)
            do c = 1, 3
                gap(c) = (mu(cg) - mu(c)) / sqrt(0.5_dp * (s2(cg) + s2(c)) + real(mv,dp))
            end do
            lj = [(c /= cg .and. pw(c) >= JUNK_MINW .and. gap(c) >= JUNK_GAP, c=1,3)]
            do i = 1, self%n
                if( best(i) == 0 ) cycle
                tc     = s2 + real(var(i),dp)
                lgc    = log(max(pw, 1.e-300_dp)) - 0.5_dp * (real(amp(i),dp) - mu)**2 / tc - 0.5_dp * log(tc)
                asg(i) = maxloc(lgc, dim=1)
            end do
            if( any(lj) )then
                do i = 1, self%n
                    if( best(i) > 0 ) post(i) = jpost(amp(i), var(i))
                end do
            endif
            lamp = best > 0 .and. post >= 0.5
            key = real(mu)
            ord = [1, 2, 3]
            call hpsort(key, ord)
            line = ''
            do c = 1, 3
                line = line//'  '//fmtr(mu(ord(c)),2)//' / '//fmtr(sqrt(s2(ord(c))),2)//' / '//fmtr(pw(ord(c)),3)
                if( ord(c) == cg )then
                    line = line//' (genuine)'
                else
                    line = line//' ('//fmtr(gap(ord(c)),1)//' sd below)'
                endif
            end do
            call say('mixture of the amplitudes (mean / intrinsic sd / weight):'//line//'; median noise sd '//&
                &fmtr(real(sqrt(mv),dp),2))
            line = ' none'
            if( any(lj) )then
                line = ''
                do c = 1, 3
                    if( lj(ord(c)) ) line = line//' '//fmtr(mu(ord(c)),2)
                end do
            endif
            call say('amplitude-rule junk components (weight >= '//fmtr(real(JUNK_MINW,dp),3)//', >= '//&
                &int2str(nint(JUNK_GAP))//' sd below the genuine one):'//line//'; '//int2str(count(lamp))//&
                &' particles (not applied)')
        endif
        ljunk = self%cpj >= 0.5
        do e = 1, 2
            idx = pack(hx(e)%v, best(hx(e)%v) > 0)
            if( size(idx) == 0 ) cycle
            med = medianr(res(idx))
            sdr = MAD_SD * medianr(abs(res(idx) - med))
            if( sdr > 0. ) rz(idx) = (res(idx) - med) / sdr
        end do
        lres = best > 0 .and. rz >= HALVES_KSIG
        wg   = (1. - self%cpj) * merge(0., 1., lres)
        call say('diagnostic, not applied: contrast against the other half''s mean map at each particle''s own noise level, '//&
            &'zero-contrast weight '//fmtr(self%cfit(1),4)//', genuine contrast median / log-sd '//fmtr(exp(self%cfit(2)),3)//&
            &' / '//fmtr(self%cfit(3),3)//'; '//int2str(count(self%cpj >= 0.5))//' at zero contrast, '//int2str(count(lres))//&
            &' with a residual '//int2str(nint(HALVES_KSIG))//' robust sd above their half''s median, '//&
            &int2str(count(wg < 0.5))//' with a genuine weight under 0.5')
        nf = count(ljunk)
        call say(int2str(nf)//' of '//int2str(self%n)//' particles flagged junk at zero contrast (posterior >= 0.5; '//&
            &fmtr(real(nf,dp)/real(self%n,dp),4)//&
            &'; half 1 '//fmtr(real(count(ljunk(hx(1)%v)),dp)/real(max(1,size(hx(1)%v)),dp),4)//', half 2 '//&
            &fmtr(real(count(ljunk(hx(2)%v)),dp)/real(max(1,size(hx(2)%v)),dp),4)//'): label 0 at this node and below')
        if( count(best == 0) > 0 ) call say(int2str(count(best == 0))//' particles not scored (no state in the other half)')
        do g = 1, nk
            x  = pack(amp, lab == g .and. best > 0)
            nf = count(lab == g .and. ljunk)
            call say('state '//int2str(g)//': '//int2str(nf)//' of '//int2str(count(lab == g))//' flagged ('//&
                &fmtr(real(nf,dp)/real(max(1,count(lab == g)),dp),4)//'), amplitude 5 % / median '//&
                &fmtr(real(quantr(x, 0.05),dp),2)//' / '//fmtr(real(quantr(x, 0.5),dp),2))
        end do
        call fclose(rep)
        call fopen(funit, string('site_junk.txt'), status='replace', action='write', iostat=ios)
        call fileiochk('junk_screen: cannot write site_junk.txt', ios)
        write(funit,'(A)') '# project row, junk (1) or not (0) at zero contrast, amplitude, noise variance, best state, '//&
            &'amplitude-rule junk posterior (not applied), contrast, contrast se, posterior-mean contrast, zero-contrast '//&
            &'posterior, residual z, genuine weight, amplitude component (1 lowest), amplitude-rule junk (not applied)'
        do i = 1, self%n
            c = 0
            if( asg(i) > 0 ) c = count(mu < mu(asg(i))) + 1
            write(funit,'(I0,1X,I0,1X,F9.4,1X,ES11.4,1X,I0,1X,F7.4,3(1X,ES11.4),1X,F7.4,1X,ES11.4,1X,F7.4,2(1X,I0))') &
                &self%rows(i), merge(1, 0, ljunk(i)), amp(i), var(i), best(i), post(i), self%cam(i), self%cse(i), &
                &self%con(i), self%cpj(i), rz(i), wg(i), c, merge(1, 0, lamp(i))
        end do
        call fclose(funit)
        where( ljunk ) seeds = 0

    contains

        subroutine say( str )
            character(len=*), intent(in) :: str
            write(logfhandle,'(A)') '>>> JUNK_SCREEN '//str
            write(rep,'(A)') str
        end subroutine say

        subroutine score( pidx, loff, e )
            integer, intent(in) :: pidx(:), e
            logical, intent(in) :: loff(:)
            real(dp) :: sa(self%kj), sb(self%kj), sf(self%kj), sc(self%kj), beta(self%kj), wsh(self%kj)
            real(dp) :: a(self%kj,nk), b(self%kj,nk), num, den, s, smax, tot, rmin
            real     :: cw(self%kj)
            integer  :: ii, i, j, g, nd
            cw = 0.
            do j = 1, self%nj
                cw(self%jsh(j)) = cw(self%jsh(j)) + self%jw(j)
            end do
            sa = 0._dp
            sb = 0._dp
            sf = 0._dp
            sc = 0._dp
            !$omp parallel do default(shared) private(ii,i,g,a,b) schedule(static) reduction(+:sa,sb,sf,sc) proc_bind(close)
            do ii = 1, size(pidx)
                i = pidx(ii)
                g = lab(i)
                if( g < 1 ) cycle
                if( .not. loff(g) ) cycle
                call self%band_terms(rot, i, rf(:,:,:,g:g), a(:,1:1), b(:,1:1))
                sa = sa + a(:,1)
                sb = sb + b(:,1)
                sf = sf + real(self%f2j(:,i),dp)
                sc = sc + real(cw,dp)
            end do
            !$omp end parallel do
            beta = min(max(sa / max(sb, 1.e-30_dp), 0._dp), 2._dp)
            wsh  = 1._dp / max((sf - 2._dp*beta*sa + beta*beta*sb) / max(sc, 1._dp), 1.e-12_dp)
            !$omp parallel do default(shared) private(ii,i,g,a,b,num,den,s,smax,tot,rmin) schedule(static) proc_bind(close)
            do ii = 1, size(pidx)
                i = pidx(ii)
                call self%band_terms(rot, i, rf, a, b)
                tot  = sum(wsh * real(self%f2j(:,i),dp))
                rmin = tot
                smax = -huge(smax)
                do g = 1, nk
                    if( .not. loff(g) ) cycle
                    num = sum(wsh * beta * a(:,g))
                    den = sum(wsh * beta * beta * b(:,g))
                    if( num > 0._dp .and. den > 0._dp ) rmin = min(rmin, tot - num * num / den)
                    s   = num - 0.5_dp * den
                    if( s <= smax ) cycle
                    smax    = s
                    best(i) = g
                    amp(i)  = real(num / max(den, 1.e-12_dp))
                    var(i)  = real(1._dp / max(den, 1.e-12_dp))
                end do
                res(i) = real(rmin) / self%nv(i)
            end do
            !$omp end parallel do
            nd = count([(best(pidx(ii)) /= lab(pidx(ii)), ii=1,size(pidx))])
            call say('half '//int2str(e)//': '//int2str(size(pidx))//' particles, best state differs from the label for '//&
                &fmtr(real(nd,dp)/real(max(1,size(pidx)),dp),3)//', amplitude median '//&
                &fmtr(real(medianr(pack(amp(pidx), best(pidx) > 0)),dp),2)//', model scale at shells 2 / '//&
                &int2str(self%kj)//' '//fmtr(beta(2),2)//' / '//fmtr(beta(self%kj),2))
        end subroutine score

        subroutine xd_fit
            real(dp) :: t(3), lg(3), r(3), bm(3), sr(3), srb(3), srb2(3), srv(3), wc(3), xi, vi
            integer  :: it, ii
            do it = 1, JUNK_ITS
                sr   = 0._dp
                srb  = 0._dp
                srb2 = 0._dp
                srv  = 0._dp
                !$omp parallel do default(shared) private(ii,xi,vi,t,lg,r,bm) schedule(static) &
                !$omp& reduction(+:sr,srb,srb2,srv) proc_bind(close)
                do ii = 1, size(x)
                    xi   = real(x(ii),dp)
                    vi   = real(v(ii),dp)
                    t    = s2 + vi
                    lg   = log(max(pw, 1.e-300_dp)) - 0.5_dp * (xi - mu)**2 / t - 0.5_dp * log(t)
                    r    = exp(lg - maxval(lg))
                    r    = r / sum(r)
                    bm   = mu + s2 / t * (xi - mu)
                    sr   = sr   + r
                    srb  = srb  + r * bm
                    srb2 = srb2 + r * bm * bm
                    srv  = srv  + r * s2 * vi / t
                end do
                !$omp end parallel do
                wc = max(sr, 1.e-12_dp)
                pw = sr / real(size(x),dp)
                mu = srb / wc
                s2 = max((srb2 - 2._dp * mu * srb + mu * mu * sr + srv) / wc, 1.e-8_dp)
            end do
        end subroutine xd_fit

        real function jpost( xi, vi )
            real, intent(in) :: xi, vi
            real(dp) :: t(3), lg(3), r(3)
            t     = s2 + real(vi,dp)
            lg    = log(max(pw, 1.e-300_dp)) - 0.5_dp * (real(xi,dp) - mu)**2 / t - 0.5_dp * log(t)
            r     = exp(lg - maxval(lg))
            jpost = real(sum(r, mask=lj) / sum(r))
        end function jpost

    end subroutine junk_screen

    function sym_rots( self ) result( rot )
        class(state_mergeup), intent(inout) :: self
        real, allocatable :: rot(:,:,:,:)
        real    :: es(3)
        integer :: i, is, nsym
        nsym = self%syms(1)%get_nsym()
        allocate(rot(3,3,nsym,self%n))
        do i = 1, self%n
            call self%os(1)%set_euler(self%eul(:,i))
            rot(:,:,1,i) = self%os(1)%get_mat()
            do is = 2, nsym
                call self%syms(1)%apply(self%eul(:,i), is, es)
                call self%os(1)%set_euler(es)
                rot(:,:,is,i) = self%os(1)%get_mat()
            end do
        end do
    end function sym_rots

    function band_density( self, rot, pidx ) result( ev )
        use simple_srch_sort_loc, only: hpsort
        class(state_mergeup), intent(in) :: self
        real,                 intent(in) :: rot(:,:,:,:)
        integer,              intent(in) :: pidx(:)
        real(dp), allocatable :: ev(:)
        real(dp), allocatable :: den(:,:,:)
        real,     allocatable :: vals(:)
        real(dp) :: c2v
        integer  :: ii, j, is, c1(3), c2(3), s, cnt
        logical  :: cnj
        allocate(den(self%nh,self%box,self%box), source=0._dp)
        !$omp parallel do default(shared) private(ii,j,is,c1,c2,cnj,c2v) schedule(static) reduction(+:den) proc_bind(close)
        do ii = 1, size(pidx)
            do j = 1, self%nj
                c2v = real(self%jw(j) * self%fpls(pidx(ii))%ctfsq_plane(self%jhk(1,j),self%jhk(2,j)), dp)
                do is = 1, size(rot,3)
                    call plane_cell(rot(:,:,is,pidx(ii)), self%jhk(1,j), self%jhk(2,j), self%nyq, self%box, c1, c2, cnj)
                    if( c1(1) == 0 ) cycle
                    if( c2(1) > 0 )then
                        den(c1(1),c1(2),c1(3)) = den(c1(1),c1(2),c1(3)) + 0.5_dp * c2v
                        den(c2(1),c2(2),c2(3)) = den(c2(1),c2(2),c2(3)) + 0.5_dp * c2v
                    else
                        den(c1(1),c1(2),c1(3)) = den(c1(1),c1(2),c1(3)) + c2v
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        allocate(ev(0:self%nyq), source=1._dp)
        do s = 0, self%nyq
            vals = real(pack(den, self%fsh == s .and. den > 0._dp))
            cnt  = size(vals)
            if( cnt == 0 ) cycle
            call hpsort(vals)
            ev(s) = real(vals((cnt + 1) / 2), dp) / real(max(1, size(pidx)), dp)
        end do
    end function band_density

    subroutine band_wiener( self, rot, pidx, epsv, out )
        class(state_mergeup), intent(in)  :: self
        real,                 intent(in)  :: rot(:,:,:,:)
        integer,              intent(in)  :: pidx(:)
        real(dp),             intent(in)  :: epsv(0:)
        complex,              intent(out) :: out(:,:,:)
        complex(dp), allocatable :: num(:,:,:)
        real(dp),    allocatable :: den(:,:,:)
        complex(dp) :: cv, cj
        real(dp)    :: c2v
        integer     :: ii, j, is, c1(3), c2(3), hh, kk, ll
        logical     :: cnj
        allocate(num(self%nh,self%box,self%box), source=cmplx(0._dp,0._dp,dp))
        allocate(den(self%nh,self%box,self%box), source=0._dp)
        !$omp parallel do default(shared) private(ii,j,is,c1,c2,cnj,cv,cj,c2v) schedule(static) reduction(+:num,den) &
        !$omp& proc_bind(close)
        do ii = 1, size(pidx)
            do j = 1, self%nj
                cj  = real(self%jw(j),dp) * cmplx(self%fpls(pidx(ii))%cmplx_plane(self%jhk(1,j),self%jhk(2,j)), kind=dp)
                c2v = real(self%jw(j) * self%fpls(pidx(ii))%ctfsq_plane(self%jhk(1,j),self%jhk(2,j)), dp)
                do is = 1, size(rot,3)
                    call plane_cell(rot(:,:,is,pidx(ii)), self%jhk(1,j), self%jhk(2,j), self%nyq, self%box, c1, c2, cnj)
                    if( c1(1) == 0 ) cycle
                    cv = cj
                    if( cnj ) cv = conjg(cv)
                    if( c2(1) > 0 )then
                        num(c1(1),c1(2),c1(3)) = num(c1(1),c1(2),c1(3)) + 0.5_dp * cv
                        den(c1(1),c1(2),c1(3)) = den(c1(1),c1(2),c1(3)) + 0.5_dp * c2v
                        num(c2(1),c2(2),c2(3)) = num(c2(1),c2(2),c2(3)) + 0.5_dp * conjg(cv)
                        den(c2(1),c2(2),c2(3)) = den(c2(1),c2(2),c2(3)) + 0.5_dp * c2v
                    else
                        num(c1(1),c1(2),c1(3)) = num(c1(1),c1(2),c1(3)) + cv
                        den(c1(1),c1(2),c1(3)) = den(c1(1),c1(2),c1(3)) + c2v
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        do ll = 1, self%box
            do kk = 1, self%box
                do hh = 1, self%nh
                    out(hh,kk,ll) = cmplx(num(hh,kk,ll) / (den(hh,kk,ll) + epsv(self%fsh(hh,kk,ll)) * &
                        &real(max(1, size(pidx)),dp)), kind=kind(out))
                end do
            end do
        end do
    end subroutine band_wiener

    subroutine band_terms( self, rot, i, refs, a, b )
        class(state_mergeup), intent(in)  :: self
        real,                 intent(in)  :: rot(:,:,:,:)
        integer,              intent(in)  :: i
        complex,              intent(in)  :: refs(:,:,:,:)
        real(dp),             intent(out) :: a(:,:), b(:,:)
        complex :: pv, dv
        real    :: c2v, ws
        integer :: j, is, g, sh, c1(3), c2(3)
        logical :: cnj
        a = 0._dp
        b = 0._dp
        do j = 1, self%nj
            sh  = self%jsh(j)
            dv  = self%fpls(i)%cmplx_plane(self%jhk(1,j),self%jhk(2,j))
            c2v = self%fpls(i)%ctfsq_plane(self%jhk(1,j),self%jhk(2,j))
            ws  = self%jw(j) / real(size(rot,3))
            do is = 1, size(rot,3)
                call plane_cell(rot(:,:,is,i), self%jhk(1,j), self%jhk(2,j), self%nyq, self%box, c1, c2, cnj)
                if( c1(1) == 0 ) cycle
                do g = 1, size(refs,4)
                    pv = refs(c1(1),c1(2),c1(3),g)
                    if( cnj ) pv = conjg(pv)
                    a(sh,g) = a(sh,g) + ws * (real(pv)*real(dv) + aimag(pv)*aimag(dv))
                    b(sh,g) = b(sh,g) + ws * c2v * abs(pv)**2
                end do
            end do
        end do
    end subroutine band_terms

    subroutine merge_poles( posts, thr, minseed, sp, grp, seeds )
        real,                 intent(in)  :: posts(:,:), thr
        integer,              intent(in)  :: minseed
        real,    allocatable, intent(out) :: sp(:,:)
        integer, allocatable, intent(out) :: grp(:), seeds(:)
        real(dp), allocatable :: rk(:,:), gpost(:,:)
        integer,  allocatable :: lab(:), plab(:)
        integer :: n, ne, a, b, ga, gb, npole
        n  = size(posts,1)
        ne = size(posts,2)
        allocate(rk(n,ne), sp(ne,ne))
        do a = 1, ne
            rk(:,a) = avg_ranks(posts(:,a))
            rk(:,a) = rk(:,a) - sum(rk(:,a)) / real(n,dp)
        end do
        do a = 1, ne
            do b = 1, ne
                sp(a,b) = real(sum(rk(:,a) * rk(:,b)) / max(sqrt(sum(rk(:,a)**2) * sum(rk(:,b)**2)), 1.e-300_dp))
            end do
        end do
        grp = [(a, a=1,ne)]
        do a = 1, ne
            do b = a + 1, ne
                if( sp(a,b) < thr ) cycle
                ga = grp(a)
                gb = grp(b)
                where( grp == gb ) grp = ga
            end do
        end do
        allocate(lab(ne), source=0)
        npole = 0
        do a = 1, ne
            if( lab(a) > 0 ) cycle
            npole = npole + 1
            where( grp == grp(a) ) lab = npole
        end do
        grp = lab
        allocate(gpost(n,npole), source=0._dp)
        do a = 1, ne
            gpost(:,grp(a)) = gpost(:,grp(a)) + real(posts(:,a),dp) / real(count(grp == grp(a)),dp)
        end do
        call seed_poles(gpost, minseed, .true., seeds, plab)
        grp = plab(grp)
    end subroutine merge_poles

    subroutine seed_poles( gpost, minseed, l_rest, seeds, plab )
        real(dp),             intent(in)  :: gpost(:,:)
        integer,              intent(in)  :: minseed
        logical,              intent(in)  :: l_rest
        integer, allocatable, intent(out) :: seeds(:), plab(:)
        integer, allocatable :: cnt(:)
        logical, allocatable :: alive(:)
        logical :: lrest
        integer :: n, np, i, a
        n  = size(gpost,1)
        np = size(gpost,2)
        allocate(seeds(n))
        allocate(alive(np), source=.true.)
        lrest = l_rest
        do
            do i = 1, n
                seeds(i) = np + 1
                if( .not. any(alive) ) cycle
                a = maxloc(gpost(i,:), dim=1, mask=alive)
                if( .not. lrest .or. gpost(i,a) >= 0.5_dp ) seeds(i) = a
            end do
            cnt = [(count(seeds == a), a=1,np+1)]
            if( count(alive) + merge(1, 0, lrest) <= 1 ) exit
            a = minloc(cnt, dim=1, mask=[alive, lrest])
            if( cnt(a) >= minseed ) exit
            if( a > np )then
                lrest = .false.
            else
                alive(a) = .false.
            endif
        end do
        allocate(plab(np+1), source=0)
        do a = 1, np
            if( alive(a) ) plab(a) = maxval(plab) + 1
        end do
        if( any(seeds == np + 1) ) plab(np+1) = maxval(plab) + 1
        seeds = plab(seeds)
        plab  = plab(1:np)
    end subroutine seed_poles

    ! cross-half correlation of deviation map a's even half with map b's odd half, from the hold-out-band Gram matrices
    real function gram_corr( gee, goo, geo, a, b )
        real(dp), intent(in) :: gee(:,:), goo(:,:), geo(:,:)
        integer,  intent(in) :: a, b
        gram_corr = real(geo(a,b) / max(sqrt(gee(a,a) * goo(b,b)), 1.e-300_dp))
    end function gram_corr

    ! Fisher z of a correlation over neff independent coefficients
    real function fisher_z( r, neff )
        real, intent(in) :: r, neff
        fisher_z = atanh(min(max(r, -0.999999), 0.999999)) * sqrt(max(neff - 3., 1.))
    end function fisher_z

    ! the principal-angle test: the cross-half correlation of map b's residual after the span of the groups' directions (each
    ! group's direction = the mean of its members' unit deviations) is projected out, in each half with that half's directions
    real function gram_resid( gee, goo, geo, b, gmem ) result( r )
        real(dp),   intent(in) :: gee(:,:), goo(:,:), geo(:,:)
        integer,    intent(in) :: b
        type(ivec), intent(in) :: gmem(:)
        real(dp), allocatable :: we(:,:), wo(:,:), ge(:,:), go(:,:), gi(:,:), re(:), ro(:)
        integer :: m, k, j, i, a, err
        m = size(gee,1)
        k = size(gmem)
        allocate(we(m,k), wo(m,k), source=0._dp)
        do j = 1, k
            do i = 1, size(gmem(j)%v)
                a = gmem(j)%v(i)
                we(a,j) = we(a,j) + 1._dp / (sqrt(gee(a,a)) * real(size(gmem(j)%v),dp))
                wo(a,j) = wo(a,j) + 1._dp / (sqrt(goo(a,a)) * real(size(gmem(j)%v),dp))
            end do
        end do
        allocate(re(m), ro(m), source=0._dp)
        re(b) = 1._dp
        ro(b) = 1._dp
        ge = matmul(transpose(we), matmul(gee, we))
        go = matmul(transpose(wo), matmul(goo, wo))
        do j = 1, k
            ge(j,j) = ge(j,j) * (1._dp + 1.e-9_dp)
            go(j,j) = go(j,j) * (1._dp + 1.e-9_dp)
        end do
        allocate(gi(k,k))
        call matinv(ge, gi, k, err)
        if( err == 0 ) re = re - matmul(we, matmul(gi, matmul(transpose(we), gee(:,b))))
        call matinv(go, gi, k, err)
        if( err == 0 ) ro = ro - matmul(wo, matmul(gi, matmul(transpose(wo), goo(:,b))))
        r = real(dot_product(re, matmul(geo, ro)) / max(sqrt(max(dot_product(re, matmul(gee, re)), 0._dp) * &
            &max(dot_product(ro, matmul(goo, ro)), 0._dp)), 1.e-300_dp))
    end function gram_resid

    ! the subspace test (audit 6.3) on the endpoints' cross-half deviation maps: an endpoint is accepted when its even-odd
    ! correlation is HALVES_KSIG Fisher SEs above zero; two accepted endpoints are one pole when they correlate positively and
    ! neither's residual after the other reproduces; groups merge closest first when every cross pair is such a pair; a pole in
    ! the span of two or more other positively correlated poles is redundant and dropped (zred)
    subroutine merge_poles_subspace( gee, goo, geo, neff, posts, minseed, grp, seeds, zself, zdist, rho, zres, zred, l_ok )
        real(dp),             intent(in)  :: gee(:,:), goo(:,:), geo(:,:)
        real,                 intent(in)  :: neff, posts(:,:)
        integer,              intent(in)  :: minseed
        integer, allocatable, intent(out) :: grp(:), seeds(:)
        real,    allocatable, intent(out) :: zself(:), zdist(:), rho(:,:), zres(:,:), zred(:)
        logical, optional,    intent(in)  :: l_ok(:)
        type(ivec), allocatable :: mem(:), oth(:)
        real(dp),   allocatable :: gpost(:,:)
        integer,    allocatable :: lab(:), plab(:), gl(:), ord(:)
        logical,    allocatable :: acc(:), same(:,:)
        type(ivec) :: one
        real    :: best, zmin, gz
        integer :: ne, n, a, b, ba, bb, npole, g, h, k, i, j
        logical :: l_all
        ne = size(gee,1)
        n  = size(posts,1)
        allocate(zself(ne), zdist(ne), zred(ne), rho(ne,ne), zres(ne,ne), acc(ne), same(ne,ne))
        zdist = huge(1.)
        zred  = huge(1.)
        rho   = 0.
        zres  = huge(1.)
        same  = .false.
        do a = 1, ne
            zself(a) = fisher_z(gram_corr(gee, goo, geo, a, a), neff)
        end do
        acc = zself >= HALVES_KSIG
        if( present(l_ok) ) acc = acc .and. l_ok
        allocate(one%v(1))
        do a = 1, ne
            do b = a + 1, ne
                if( .not. (acc(a) .and. acc(b)) ) cycle
                rho(a,b) = 0.5 * (gram_corr(gee, goo, geo, a, b) + gram_corr(gee, goo, geo, b, a)) / &
                    &sqrt(max(gram_corr(gee, goo, geo, a, a) * gram_corr(gee, goo, geo, b, b), 1.e-12))
                rho(b,a) = rho(a,b)
                one%v(1) = b
                zres(a,b) = fisher_z(gram_resid(gee, goo, geo, a, [one]), neff - 1.)
                one%v(1) = a
                zres(b,a) = fisher_z(gram_resid(gee, goo, geo, b, [one]), neff - 1.)
                same(a,b) = rho(a,b) > 0. .and. zres(a,b) < HALVES_KSIG .and. zres(b,a) < HALVES_KSIG
                same(b,a) = same(a,b)
            end do
        end do
        lab = [(a, a=1,ne)]
        where( .not. acc ) lab = 0
        do
            best = -huge(best)
            ba   = 0
            bb   = 0
            do a = 1, ne
                do b = a + 1, ne
                    if( .not. same(a,b) .or. lab(a) == 0 .or. lab(b) == 0 .or. lab(a) == lab(b) ) cycle
                    if( rho(a,b) <= best ) cycle
                    l_all = .true.
                    do i = 1, ne
                        do j = 1, ne
                            if( lab(i) == lab(a) .and. lab(j) == lab(b) ) l_all = l_all .and. same(i,j)
                        end do
                    end do
                    if( .not. l_all ) cycle
                    best = rho(a,b)
                    ba   = a
                    bb   = b
                end do
            end do
            if( ba == 0 ) exit
            where( lab == lab(bb) ) lab = lab(ba)
        end do
        do a = 1, ne
            if( .not. acc(a) ) cycle
            do b = 1, ne
                if( b == a .or. .not. acc(b) .or. lab(b) == lab(a) .or. rho(a,b) <= 0. ) cycle
                zdist(a) = min(zdist(a), max(zres(a,b), zres(b,a)))
            end do
            if( zdist(a) == huge(1.) ) zdist(a) = zself(a)
        end do
        ! redundancy, weakest group first
        gl = pack([(a, a=1,ne)], [(lab(a) > 0 .and. .not. any(lab(1:a-1) == lab(a)), a=1,ne)])
        allocate(mem(size(gl)))
        do g = 1, size(gl)
            mem(g)%v = pack([(a, a=1,ne)], lab == lab(gl(g)))
        end do
        ord = [(g, g=1,size(gl))]
        do g = 2, size(gl)
            h = g
            do while( h > 1 )
                if( group_z(ord(h-1)) <= group_z(ord(h)) ) exit
                ord(h-1:h) = [ord(h), ord(h-1)]
                h = h - 1
            end do
        end do
        do i = 1, size(ord)
            g = ord(i)
            if( lab(mem(g)%v(1)) == 0 ) cycle
            allocate(oth(0))
            do h = 1, size(gl)
                if( h == g .or. lab(mem(h)%v(1)) == 0 ) cycle
                gz = sum([(sum(rho(mem(g)%v, mem(h)%v(k))), k=1,size(mem(h)%v))]) / real(size(mem(g)%v) * size(mem(h)%v))
                if( gz > 0. ) oth = [oth, mem(h)]
            end do
            if( size(oth) >= 2 )then
                zmin = huge(zmin)
                do k = 1, size(mem(g)%v)
                    zred(mem(g)%v(k)) = fisher_z(gram_resid(gee, goo, geo, mem(g)%v(k), oth), neff - real(size(oth)))
                    zmin = min(zmin, zred(mem(g)%v(k)))
                end do
                if( zmin < HALVES_KSIG ) lab(mem(g)%v) = 0
            endif
            deallocate(oth)
        end do
        allocate(grp(ne), source=0)
        npole = 0
        do a = 1, ne
            if( lab(a) == 0 .or. grp(a) > 0 ) cycle
            npole = npole + 1
            where( lab == lab(a) ) grp = npole
        end do
        allocate(seeds(n), source=0)
        if( npole == 0 ) return
        allocate(gpost(n,npole), source=0._dp)
        do a = 1, ne
            if( grp(a) == 0 ) cycle
            gpost(:,grp(a)) = gpost(:,grp(a)) + real(posts(:,a),dp) / real(count(grp == grp(a)),dp)
        end do
        call seed_poles(gpost, minseed, .true., seeds, plab)
        do a = 1, ne
            if( grp(a) > 0 ) grp(a) = plab(grp(a))
        end do

    contains

        real function group_z( g )
            integer, intent(in) :: g
            group_z = sum(zself(mem(g)%v)) / real(size(mem(g)%v))
        end function group_z

    end subroutine merge_poles_subspace

    subroutine estep( self, idx, w, refs, l_alive, wnew, l_full, ll, kcut, llband, lpp, l_flat )
        use simple_srch_sort_loc, only: hpsort
        class(state_mergeup), target, intent(inout) :: self
        integer,              intent(in)    :: idx(:)
        real,                 intent(in)    :: w(:,:)
        type(reconstructor),  intent(in)    :: refs(:)
        logical,              intent(in)    :: l_alive(:)
        real, allocatable,    intent(out)   :: wnew(:,:)
        logical,  optional,   intent(in)    :: l_full
        real(dp), optional,   intent(out)   :: ll
        integer,  optional,   intent(in)    :: kcut
        real(dp), optional,   intent(out)   :: llband(2)
        real(dp), allocatable, optional, intent(out) :: lpp(:)
        logical,  optional,   intent(in)    :: l_flat
        real(dp), allocatable :: lpl(:)
        real,    allocatable :: aa(:,:,:), bb(:,:,:), gap(:), agap(:)
        integer, allocatable :: gidx(:)
        complex, allocatable :: projected(:,:)
        integer, allocatable :: mp(:)
        integer, pointer :: bhk(:,:), bsh(:)
        real,    pointer :: bw(:), bf2(:,:), bcw(:)
        real(dp) :: sa(max(self%nshem,self%nshfb)), sb(max(self%nshem,self%nshfb)), sf(max(self%nshem,self%nshfb)), &
            &sc(max(self%nshem,self%nshfb)), scl(max(self%nshem,self%nshfb)), sig2(max(self%nshem,self%nshfb))
        real(dp) :: prior(size(refs)), s(size(refs)), smax, lse, num, sbd(size(refs),2), lseb(2), bmax
        complex  :: pv, wv
        integer  :: n, k, i, ip, g, j, h, kk, sh, it, nyq_save, nsb, nsmp, kband, ic, ib, nthr
        logical  :: l_fb, l_bd, l_reuse, l_fl
        l_fb = .false.
        if( present(l_full) ) l_fb = l_full
        l_fl = .false.
        if( present(l_flat) ) l_fl = l_flat
        l_fl = l_fl .and. size(refs) == 2
        if( l_fl ) allocate(gap(size(idx)), source=0.)
        if( l_fb )then
            bhk => self%fbhk
            bsh => self%fbsh
            bw  => self%fbw
            bf2 => self%f2b
            bcw => self%cwb
            nsb   = self%nshfb
            nsmp  = self%nfb
            kband = self%khi
        else
            bhk => self%emhk
            bsh => self%emsh
            bw  => self%emw
            bf2 => self%f2s
            bcw => self%cw
            nsb   = self%nshem
            nsmp  = self%nem
            kband = self%kb
        endif
        n  = size(idx)
        k  = size(refs)
        mp = maxloc(w, dim=2)
        prior = 0._dp
        if( .not. l_fl )then
            do g = 1, k
                prior(g) = log(max(real(sum(w(:,g)),dp) / real(n,dp), 1.e-6_dp))
            end do
        endif
        allocate(aa(nsb,k,n), bb(nsb,k,n), source=0.)
        l_reuse = self%p%merge_geom_reuse == 'yes'
        if( l_reuse )then
            nthr = 1
            !$ nthr = omp_get_max_threads()
            if( nthr > size(self%smps) ) THROW_HARD('projection geometry scratch is smaller than the OpenMP team')
            ! The half-disc fits inside this rectangle, including unscored low shells.
            allocate(projected((2*kband+1)*(kband+1),nthr))
        endif
        sa = 0._dp
        sb = 0._dp
        sf = 0._dp
        sc = 0._dp
        !$omp parallel do default(shared) private(i,ip,it,g,j,h,kk,sh,pv,wv,nyq_save) schedule(static)&
        !$omp& reduction(+:sa,sb,sf,sc) proc_bind(close)
        do i = 1, n
            it = omp_thread()
            ip = idx(i)
            call self%os(it)%set_euler(self%eul(:,ip))
            nyq_save = self%fpls(ip)%nyq
            self%fpls(ip)%nyq = OSMPL_PAD_FAC * kband
            ! cross_em creates all refs on one lattice with one kernel; only their values differ.
            if( l_reuse ) call self%smps(it)%new_plane(refs(1), self%os(it), self%fpls(ip))
            do g = 1, k
                if( .not. l_alive(g) .and. g /= mp(i) ) cycle
                if( l_reuse )then
                    call self%smps(it)%gather(refs(g), projected(:,it))
                    call self%smps(it)%put_plane(projected(:,it), self%fpls(ip), self%pouts(it), .false.)
                else
                    call refs(g)%project_fplane(self%os(it), self%fpls(ip), self%pouts(it), apply_ctf_amp=.false., &
                        &samples=self%smps(it))
                endif
                do j = 1, nsmp
                    h  = bhk(1,j)
                    kk = bhk(2,j)
                    sh = bsh(j)
                    pv = self%pouts(it)%cmplx_plane(h,kk)
                    wv = self%fpls(ip)%cmplx_plane(h,kk)
                    aa(sh,g,i) = aa(sh,g,i) + bw(j) * (real(pv)*real(wv) + aimag(pv)*aimag(wv))
                    bb(sh,g,i) = bb(sh,g,i) + bw(j) * self%fpls(ip)%ctfsq_plane(h,kk) * abs(pv)**2
                end do
            end do
            self%fpls(ip)%nyq = nyq_save
            sa(:nsb) = sa(:nsb) + real(aa(:,mp(i),i),dp)
            sb(:nsb) = sb(:nsb) + real(bb(:,mp(i),i),dp)
            sf(:nsb) = sf(:nsb) + real(bf2(:,ip),dp)
            sc(:nsb) = sc(:nsb) + real(bcw,dp)
        end do
        !$omp end parallel do
        scl(:nsb)  = sa(:nsb) / max(sb(:nsb), 1.e-30_dp)
        sig2(:nsb) = max((sf(:nsb) - 2._dp*scl(:nsb)*sa(:nsb) + scl(:nsb)*scl(:nsb)*sb(:nsb)) / max(sc(:nsb), 1._dp), 1.e-12_dp)
        allocate(wnew(n,k))
        ! with kcut, the same model's log-likelihood also over the shells below kcut and from kcut up (sh = shell - 1)
        l_bd = present(kcut) .and. present(llband)
        ic   = nsb + 1
        if( l_bd ) ic = min(max(kcut - 1, 1), nsb + 1)
        lse  = 0._dp
        lseb = 0._dp
        allocate(lpl(n), source=0._dp)
        !$omp parallel do default(shared) private(i,g,s,smax,num,sbd,ib,bmax) schedule(static) reduction(+:lse,lseb) proc_bind(close)
        do i = 1, n
            do g = 1, k
                if( l_alive(g) .and. self%l_amp )then
                    ! l_amp: the particle's amplitude in the group is fitted (at least 0), so a group cannot win on brightness
                    num  = sum(scl(:nsb)*real(aa(:,g,i),dp) / sig2(:nsb))
                    s(g) = 0.5_dp * max(num, 0._dp)**2 / max(sum(scl(:nsb)**2*real(bb(:,g,i),dp) / sig2(:nsb)), 1.e-30_dp) + prior(g)
                else if( l_alive(g) )then
                    s(g) = 0.5_dp * sum((2._dp*scl(:nsb)*real(aa(:,g,i),dp) - scl(:nsb)*scl(:nsb)*real(bb(:,g,i),dp)) / &
                        &sig2(:nsb)) + prior(g)
                else
                    s(g) = -huge(1._dp)
                endif
                if( l_bd .and. l_alive(g) )then
                    sbd(g,1) = band_score(g, i, 1, ic - 1)
                    sbd(g,2) = band_score(g, i, ic, nsb)
                endif
            end do
            if( l_fl ) gap(i) = real(s(1) - s(2))
            smax = maxval(s)
            do g = 1, k
                wnew(i,g) = 0.
                if( l_alive(g) ) wnew(i,g) = real(exp(s(g) - smax))
            end do
            if( sum(wnew(i,:)) > 0. )then
                lse = lse + smax + log(real(sum(wnew(i,:)),dp))
                wnew(i,:) = wnew(i,:) / sum(wnew(i,:))
            else
                wnew(i,:) = w(i,:)
            endif
            if( l_bd )then
                do ib = 1, 2
                    bmax = maxval(sbd(:,ib), mask=l_alive)
                    num  = 0._dp
                    do g = 1, k
                        if( l_alive(g) ) num = num + exp(sbd(g,ib) - bmax)
                    end do
                    lseb(ib) = lseb(ib) + bmax + log(num)
                end do
                lpl(i) = bmax + log(num) - 0.5_dp * sum(real(bf2(ic:nsb,idx(i)),dp) / sig2(ic:nsb) + real(bcw(ic:nsb),dp) * &
                    &log(sig2(ic:nsb)))
            endif
        end do
        !$omp end parallel do
        if( present(ll) ) ll = lse - 0.5_dp * sum(sf(:nsb) / sig2(:nsb) + sc(:nsb) * log(sig2(:nsb)))
        if( l_bd )then
            llband(1) = lseb(1) - 0.5_dp * sum(sf(:ic-1) / sig2(:ic-1) + sc(:ic-1) * log(sig2(:ic-1)))
            llband(2) = lseb(2) - 0.5_dp * sum(sf(ic:nsb) / sig2(ic:nsb) + sc(ic:nsb) * log(sig2(ic:nsb)))
        endif
        if( present(lpp) ) call move_alloc(lpl, lpp)
        write(logfhandle,'(A)') '  cross-half E-step: model scale at shells 2 / '//int2str(kband)//' '//&
            &fmtr(scl(1),4)//' / '//fmtr(scl(nsb),4)
        if( l_fl )then
            ! the site path's evidence per particle, as xsbg.py reports it: the two-group score gap's P10 / median / P90 and median |gap|
            gidx = [(i, i=1,n)]
            agap = abs(gap)
            call hpsort(gap, gidx)
            gidx = [(i, i=1,n)]
            call hpsort(agap, gidx)
            write(logfhandle,'(A)') '  cross-half E-step: score gap P10 '//fmtr(real(gap(max(1,n/10)),dp),2)//' median '//&
                &fmtr(real(gap((n+1)/2),dp),2)//' P90 '//fmtr(real(gap(max(1,(9*n)/10)),dp),2)//'; median |gap| '//&
                &fmtr(real(agap((n+1)/2),dp),2)
        endif
        deallocate(aa, bb, mp)

    contains

        real(dp) function band_score( gg, ii, s0, s1 )
            integer, intent(in) :: gg, ii, s0, s1
            real(dp) :: nm
            if( self%l_amp )then
                nm = sum(scl(s0:s1)*real(aa(s0:s1,gg,ii),dp) / sig2(s0:s1))
                band_score = 0.5_dp * max(nm, 0._dp)**2 / max(sum(scl(s0:s1)**2*real(bb(s0:s1,gg,ii),dp) / sig2(s0:s1)), 1.e-30_dp)
            else
                band_score = 0.5_dp * sum((2._dp*scl(s0:s1)*real(aa(s0:s1,gg,ii),dp) - scl(s0:s1)**2*real(bb(s0:s1,gg,ii),dp)) / &
                    &sig2(s0:s1))
            endif
            band_score = band_score + prior(gg)
        end function band_score

    end subroutine estep

    subroutine accumulate( self, idx, it, l_band, w )
        class(state_mergeup), intent(inout) :: self
        integer,              intent(in)    :: idx(:), it
        logical,              intent(in)    :: l_band
        real, optional,       intent(in)    :: w(:)
        integer :: i, j, h, k
        call self%recs(it)%reset
        call self%recs(it)%reset_exp
        do i = 1, size(idx)
            if( present(w) )then
                if( w(i) <= 0. ) cycle
            endif
            call self%os(it)%set_euler(self%eul(:,idx(i)))
            if( l_band )then
                do j = 1, self%nem
                    h = self%emhk(1,j)
                    k = self%emhk(2,j)
                    self%bpls(it)%cmplx_plane(h,k) = self%fpls(idx(i))%cmplx_plane(h,k)
                    self%bpls(it)%ctfsq_plane(h,k) = self%fpls(idx(i))%ctfsq_plane(h,k)
                end do
                if( present(w) )then
                    call self%recs(it)%insert_plane_oversamp(self%syms(it), self%os(it), self%bpls(it), w=w(i))
                else
                    call self%recs(it)%insert_plane_oversamp(self%syms(it), self%os(it), self%bpls(it))
                endif
            else if( present(w) )then
                call self%recs(it)%insert_plane_oversamp(self%syms(it), self%os(it), self%fpls(idx(i)), w=w(i))
            else
                call self%recs(it)%insert_plane_oversamp(self%syms(it), self%os(it), self%fpls(idx(i)))
            endif
        end do
        call self%recs(it)%compress_exp
    end subroutine accumulate

    ! job j's particles jidx(joff(j):joff(j+1)-1) cut into chunks accumulated in parallel; its sum ends in slot js(j)
    subroutine accumulate_jobs( self, jidx, joff, l_band, js, jw, nc_used )
        class(state_mergeup), intent(inout) :: self
        integer,              intent(in)    :: jidx(:), joff(:)
        logical,              intent(in)    :: l_band
        integer, allocatable, intent(out)   :: js(:)
        real, optional,       intent(in)    :: jw(:)
        integer, optional,    intent(out)   :: nc_used
        integer, allocatable :: cfrom(:), cto(:)
        integer :: nj, nc, nthr, j, c, ic
        nj   = size(joff) - 1
        nthr = 1
        !$ nthr = omp_get_max_threads()
        nc = max(1, min(self%nslots, max(nthr, nj)) / nj)
        if( present(nc_used) ) nc_used = nc
        if( nj * nc > self%nslots ) THROW_HARD('more jobs than slots; accumulate_jobs')
        allocate(cfrom(nj*nc), cto(nj*nc), js(nj))
        do j = 1, nj
            js(j) = (j - 1) * nc + 1
            do c = 1, nc
                ic = js(j) + c - 1
                cfrom(ic) = joff(j) + ((c - 1) * (joff(j+1) - joff(j))) / nc
                cto(ic)   = joff(j) + (c * (joff(j+1) - joff(j))) / nc - 1
            end do
        end do
        !$omp parallel do default(shared) private(ic) schedule(static,1) proc_bind(close)
        do ic = 1, nj * nc
            if( present(jw) )then
                call self%accumulate(jidx(cfrom(ic):cto(ic)), ic, l_band, jw(cfrom(ic):cto(ic)))
            else
                call self%accumulate(jidx(cfrom(ic):cto(ic)), ic, l_band)
            endif
        end do
        !$omp end parallel do
        !$omp parallel do default(shared) private(j,c) schedule(static,1) proc_bind(close)
        do j = 1, nj
            do c = 2, nc
                call self%recs(js(j))%sum_reduce(self%recs(js(j)+c-1))
            end do
        end do
        !$omp end parallel do
    end subroutine accumulate_jobs

    subroutine finish( self, it, epsv, l_mask, cm )
        class(state_mergeup), intent(inout) :: self
        integer,              intent(in)    :: it
        real,                 intent(in)    :: epsv(0:)
        logical,              intent(in)    :: l_mask
        complex,              intent(out)   :: cm(:,:,:)
        real, pointer :: rmat(:,:,:)
        complex, allocatable :: cu(:,:,:), ch(:,:,:)
        call self%recs(it)%add_shell_rho(epsv)
        call self%recs(it)%sampl_dens_correct
        call self%recs(it)%ifft
        call self%recs(it)%deapodize
        call self%wimgs(it)%copy_fast(self%recs(it))
        if( l_mask )then
            call self%wimgs(it)%get_rmat_ptr(rmat)
            rmat = rmat * self%msk3
        endif
        call self%wimgs(it)%fft
        cm = self%wimgs(it)%get_cmat()
    end subroutine finish

    subroutine maps( self, jidx, joff, jscale, eps, l_mask, l_band, cms, jw, nc_used )
        class(state_mergeup), intent(inout) :: self
        integer,              intent(in)    :: jidx(:), joff(:)
        real,                 intent(in)    :: jscale(:), eps(0:)
        logical,              intent(in)    :: l_mask, l_band
        complex,              intent(inout) :: cms(:,:,:,:)
        real, optional,       intent(in)    :: jw(:)
        integer, optional,    intent(out)   :: nc_used
        integer, allocatable :: js(:)
        integer :: j, it, nthr
        nthr = 1
        !$ nthr = omp_get_max_threads()
        if( present(nc_used) ) nc_used = 1
        if( size(jscale) < nthr .and. size(jscale) <= self%nslots )then
            call self%accumulate_jobs(jidx, joff, l_band, js, jw, nc_used)
            !$omp parallel do default(shared) private(j) schedule(static,1) proc_bind(close)
            do j = 1, size(jscale)
                call self%finish(js(j), eps * jscale(j), l_mask, cms(:,:,:,j))
            end do
            !$omp end parallel do
            return
        endif
        !$omp parallel do default(shared) private(j,it) schedule(dynamic,1) proc_bind(close)
        do j = 1, size(jscale)
            it = omp_thread()
            if( present(jw) )then
                call self%accumulate(jidx(joff(j):joff(j+1)-1), it, l_band, jw(joff(j):joff(j+1)-1))
            else
                call self%accumulate(jidx(joff(j):joff(j+1)-1), it, l_band)
            endif
            call self%finish(it, eps * jscale(j), l_mask, cms(:,:,:,j))
        end do
        !$omp end parallel do
    end subroutine maps

    ! EM targets share an ordered row list. Keep the original per-map chunks and reduction order.
    subroutine maps_weighted_shared( self, idx, w, jscale, eps, l_band, cms, coverage )
        use simple_reconstructor, only: insert_planes_targets, MAX_INSERT_TARGETS, MAX_INSERT_PLANES
        class(state_mergeup), intent(inout) :: self
        integer,              intent(in)    :: idx(:)
        real,                 intent(in)    :: w(:,:), jscale(:), eps(0:)
        logical,              intent(in)    :: l_band
        complex,              intent(inout) :: cms(:,:,:,:)
        integer, optional,    intent(out) :: coverage(3)
        class(ori), allocatable :: poses(:)
        type(fplane_type), allocatable :: planes(:)
        integer :: slots(MAX_INSERT_TARGETS), nj, nc, nthr, first, nt, q, c, i, j, h, kk, lo, hi, it, ib, nr, begin
        real :: weights(MAX_INSERT_TARGETS,MAX_INSERT_PLANES)
        nj = size(jscale)
        if( size(w,1) /= size(idx) .or. size(w,2)+1 /= nj ) THROW_HARD('EM weight shape; maps_weighted_shared')
        nthr = 1
        !$ nthr = omp_get_max_threads()
        nc = 1
        if( nj < nthr .and. nj <= self%nslots ) nc = max(1, min(self%nslots, max(nthr,nj)) / nj)
        if( present(coverage) ) coverage = [nc, 0, 0]
        allocate(poses(MAX_INSERT_PLANES), planes(MAX_INSERT_PLANES))
        do ib = 1, MAX_INSERT_PLANES
            call poses(ib)%new(is_ptcl=.true.)
            if( l_band )then
                planes(ib) = self%fpls(1)
                planes(ib)%cmplx_plane = cmplx(0.,0.)
                planes(ib)%ctfsq_plane = 0.
            endif
        end do
        first = 1
        do while( first <= nj )
            if( present(coverage) ) coverage(2) = coverage(2) + 1
            nt = min(MAX_INSERT_TARGETS, self%nslots/nc, nj-first+1)
            if( nt < 1 ) THROW_HARD('no target slots; maps_weighted_shared')
            do it = 1, nt*nc
                call self%recs(it)%reset
                call self%recs(it)%reset_exp
            end do
            do c = 1, nc
                slots(:nt) = [((q-1)*nc+c, q=1,nt)]
                it = slots(1)
                lo = ((c-1)*size(idx))/nc + 1
                hi = (c*size(idx))/nc
                do begin = lo, hi, MAX_INSERT_PLANES
                    nr = min(MAX_INSERT_PLANES,hi-begin+1)
                    do ib = 1, nr
                        i = begin+ib-1
                        do q = 1, nt
                            j = first+q-1
                            weights(q,ib) = 1.
                            if( j < nj ) weights(q,ib) = w(i,j)
                        end do
                        call poses(ib)%set_euler(self%eul(:,idx(i)))
                        if( l_band )then
                            do j = 1, self%nem
                                h  = self%emhk(1,j)
                                kk = self%emhk(2,j)
                                planes(ib)%cmplx_plane(h,kk) = self%fpls(idx(i))%cmplx_plane(h,kk)
                                planes(ib)%ctfsq_plane(h,kk) = self%fpls(idx(i))%ctfsq_plane(h,kk)
                            end do
                        else
                            planes(ib) = self%fpls(idx(i))
                        endif
                    end do
                    call insert_planes_targets(self%recs, slots(:nt), self%syms(it), poses, planes, weights(:nt,:nr), nr)
                    if( present(coverage) ) coverage(3) = coverage(3) + 1
                end do
            end do
            !$omp parallel do default(shared) private(q,c,it) schedule(static,1) proc_bind(close)
            do q = 1, nt
                it = (q-1)*nc+1
                do c = 0, nc-1
                    call self%recs(it+c)%compress_exp
                end do
                do c = 1, nc-1
                    call self%recs(it)%sum_reduce(self%recs(it+c))
                end do
                call self%finish(it, eps*jscale(first+q-1), .false., cms(:,:,:,first+q-1))
            end do
            !$omp end parallel do
            first = first + nt
        end do
        do ib = 1, MAX_INSERT_PLANES
            call poses(ib)%kill()
        end do
    end subroutine maps_weighted_shared

    function strata( self, idx ) result( st )
        use simple_srch_sort_loc, only: hpsort
        class(state_mergeup), intent(in) :: self
        integer,              intent(in) :: idx(:)
        integer, allocatable :: st(:)
        real,    allocatable :: dfs(:)
        real    :: q(2), pos
        integer :: i, n, lo
        n   = size(idx)
        dfs = self%dfm(idx)
        call hpsort(dfs)
        do i = 1, 2
            pos  = real(n-1) * real(i) / 3.
            lo   = int(pos)
            q(i) = dfs(lo+1) + (pos - real(lo)) * (dfs(min(lo+2,n)) - dfs(lo+1))
        end do
        allocate(st(n))
        do i = 1, n
            st(i) = (self%vbin(idx(i)) - 1) * 3 + count(q < self%dfm(idx(i)))
        end do
    end function strata

    subroutine dock( self, cm, cref, b1, b2 )
        class(state_mergeup), intent(in)    :: self
        complex,              intent(inout) :: cm(:,:,:)
        complex,              intent(in)    :: cref(:,:,:)
        integer,              intent(in)    :: b1, b2
        real(dp),    parameter :: TWOPI = 2._dp * real(PI,dp)
        complex(dp) :: ex(self%nh), ey(self%box), ez(self%box), z
        real(dp)    :: t(3), g(3), hm(3,3), st(3), q(3), det, zr, s11, s12, s13, s22, s23, s33
        integer     :: it, h, k, l, sh
        t = 0._dp
        do it = 1, NDOCK_ITS
            call phases
            g   = 0._dp
            s11 = 0._dp
            s12 = 0._dp
            s13 = 0._dp
            s22 = 0._dp
            s23 = 0._dp
            s33 = 0._dp
            do l = 1, self%box
                do k = 1, self%box
                    do h = 1, self%nh
                        sh = self%fsh(h,k,l)
                        if( sh < b1 .or. sh > b2 ) cycle
                        z = real(self%fw(h,k,l),dp) * cmplx(cm(h,k,l),kind=dp) * conjg(cmplx(cref(h,k,l),kind=dp)) * &
                            &ex(h) * ey(k) * ez(l)
                        q   = [self%fqh(h), self%fqk(k), self%fqk(l)]
                        g   = g + aimag(z) * q
                        zr  = real(z,dp)
                        s11 = s11 + zr * q(1) * q(1)
                        s12 = s12 + zr * q(1) * q(2)
                        s13 = s13 + zr * q(1) * q(3)
                        s22 = s22 + zr * q(2) * q(2)
                        s23 = s23 + zr * q(2) * q(3)
                        s33 = s33 + zr * q(3) * q(3)
                    end do
                end do
            end do
            g  = TWOPI * g
            hm = -TWOPI**2 * reshape([s11, s12, s13, s12, s22, s23, s13, s23, s33], [3,3])
            det = hm(1,1)*(hm(2,2)*hm(3,3) - hm(2,3)*hm(3,2)) - hm(1,2)*(hm(2,1)*hm(3,3) - hm(2,3)*hm(3,1)) + &
                &hm(1,3)*(hm(2,1)*hm(3,2) - hm(2,2)*hm(3,1))
            if( abs(det) < tiny(det) ) exit
            st(1) = (g(1)*(hm(2,2)*hm(3,3) - hm(2,3)*hm(3,2)) - hm(1,2)*(g(2)*hm(3,3) - hm(2,3)*g(3)) + &
                &hm(1,3)*(g(2)*hm(3,2) - hm(2,2)*g(3))) / det
            st(2) = (hm(1,1)*(g(2)*hm(3,3) - hm(2,3)*g(3)) - g(1)*(hm(2,1)*hm(3,3) - hm(2,3)*hm(3,1)) + &
                &hm(1,3)*(hm(2,1)*g(3) - g(2)*hm(3,1))) / det
            st(3) = (hm(1,1)*(hm(2,2)*g(3) - g(2)*hm(3,2)) - hm(1,2)*(hm(2,1)*g(3) - g(2)*hm(3,1)) + &
                &g(1)*(hm(2,1)*hm(3,2) - hm(2,2)*hm(3,1))) / det
            t = t - st
            if( maxval(abs(st)) < 1.e-4_dp ) exit
        end do
        if( any(t /= t) .or. maxval(abs(t)) > real(DOCK_MAX,dp) ) return
        call phases
        do l = 1, self%box
            do k = 1, self%box
                do h = 1, self%nh
                    cm(h,k,l) = cmplx(cmplx(cm(h,k,l),kind=dp) * ex(h) * ey(k) * ez(l))
                end do
            end do
        end do

    contains

        subroutine phases
            integer :: i
            do i = 1, self%nh
                ex(i) = exp(cmplx(0._dp, -TWOPI * self%fqh(i) * t(1), kind=dp))
            end do
            do i = 1, self%box
                ey(i) = exp(cmplx(0._dp, -TWOPI * self%fqk(i) * t(2), kind=dp))
                ez(i) = exp(cmplx(0._dp, -TWOPI * self%fqk(i) * t(3), kind=dp))
            end do
        end subroutine phases

    end subroutine dock

    subroutine shell_fsc( self, x, y, f )
        class(state_mergeup), intent(in)  :: self
        complex,              intent(in)  :: x(:,:,:), y(:,:,:)
        real,                 intent(out) :: f(0:)
        real(dp) :: sxy(0:self%nyq), sxx(0:self%nyq), syy(0:self%nyq)
        integer  :: h, k, l, sh
        sxy = 0._dp
        sxx = 0._dp
        syy = 0._dp
        do l = 1, self%box
            do k = 1, self%box
                do h = 1, self%nh
                    sh = self%fsh(h,k,l)
                    sxy(sh) = sxy(sh) + self%fw(h,k,l) * (real(real(x(h,k,l)),dp)*real(real(y(h,k,l)),dp) + &
                        &real(aimag(x(h,k,l)),dp)*real(aimag(y(h,k,l)),dp))
                    sxx(sh) = sxx(sh) + self%fw(h,k,l) * real(abs(x(h,k,l)),dp)**2
                    syy(sh) = syy(sh) + self%fw(h,k,l) * real(abs(y(h,k,l)),dp)**2
                end do
            end do
        end do
        f = real(sxy / max(sqrt(sxx*syy), 1.e-30_dp))
    end subroutine shell_fsc

    real function half_rel( self, a )
        class(state_mergeup), intent(in) :: self
        complex,              intent(in) :: a(:,:,:,:)
        real :: r
        r        = self%shell_corr(a(:,:,:,1), a(:,:,:,2), self%kb + RB_MARGIN, self%khi)
        half_rel = 2. * r / max(1. + r, 1.e-6)
    end function half_rel

    real function half_corr( self, a, b )
        class(state_mergeup), intent(in) :: self
        complex,              intent(in) :: a(:,:,:,:), b(:,:,:,:)
        half_corr = self%shell_corr(a(:,:,:,1) + a(:,:,:,2), b(:,:,:,1) + b(:,:,:,2), self%kb + RB_MARGIN, self%khi)
    end function half_corr

    ! Gram matrices of deviation maps over the hold-out band: within the even half, within the odd half, even by odd
    subroutine dev_gram( self, d, gee, goo, geo )
        class(state_mergeup),  intent(in)  :: self
        complex,               intent(in)  :: d(:,:,:,:,:)
        real(dp), allocatable, intent(out) :: gee(:,:), goo(:,:), geo(:,:)
        integer :: a, b, m
        m = size(d,5)
        allocate(gee(m,m), goo(m,m), geo(m,m), source=0._dp)
        do a = 1, m
            do b = 1, m
                geo(a,b) = self%shell_dot(d(:,:,:,1,a), d(:,:,:,2,b), self%kb + RB_MARGIN, self%khi)
                if( b < a ) cycle
                gee(a,b) = self%shell_dot(d(:,:,:,1,a), d(:,:,:,1,b), self%kb + RB_MARGIN, self%khi)
                goo(a,b) = self%shell_dot(d(:,:,:,2,a), d(:,:,:,2,b), self%kb + RB_MARGIN, self%khi)
                gee(b,a) = gee(a,b)
                goo(b,a) = goo(a,b)
            end do
        end do
    end subroutine dev_gram

    real(dp) function shell_dot( self, x, y, shlo, shhi )
        class(state_mergeup), intent(in) :: self
        complex,              intent(in) :: x(:,:,:), y(:,:,:)
        integer,              intent(in) :: shlo, shhi
        integer :: h, k, l, sh
        shell_dot = 0._dp
        do l = 1, self%box
            do k = 1, self%box
                do h = 1, self%nh
                    sh = self%fsh(h,k,l)
                    if( sh < shlo .or. sh > shhi ) cycle
                    shell_dot = shell_dot + self%fw(h,k,l) * (real(real(x(h,k,l)),dp)*real(real(y(h,k,l)),dp) + &
                        &real(aimag(x(h,k,l)),dp)*real(aimag(y(h,k,l)),dp))
                end do
            end do
        end do
    end function shell_dot

    real function shell_corr( self, x, y, shlo, shhi )
        class(state_mergeup), intent(in) :: self
        complex,              intent(in) :: x(:,:,:), y(:,:,:)
        integer,              intent(in) :: shlo, shhi
        real(dp) :: sxy, sxx, syy
        integer  :: h, k, l, sh
        sxy = 0._dp
        sxx = 0._dp
        syy = 0._dp
        do l = 1, self%box
            do k = 1, self%box
                do h = 1, self%nh
                    sh = self%fsh(h,k,l)
                    if( sh < shlo .or. sh > shhi ) cycle
                    sxy = sxy + self%fw(h,k,l) * (real(real(x(h,k,l)),dp)*real(real(y(h,k,l)),dp) + &
                        &real(aimag(x(h,k,l)),dp)*real(aimag(y(h,k,l)),dp))
                    sxx = sxx + self%fw(h,k,l) * real(abs(x(h,k,l)),dp)**2
                    syy = syy + self%fw(h,k,l) * real(abs(y(h,k,l)),dp)**2
                end do
            end do
        end do
        shell_corr = real(sxy / max(sqrt(sxx*syy), 1.e-30_dp))
    end function shell_corr

    subroutine kill( self )
        class(state_mergeup), intent(inout) :: self
        integer :: i
        if( allocated(self%fpls) ) deallocate(self%fpls)
        if( allocated(self%recs) )then
            do i = 1, size(self%recs)
                call self%recs(i)%kill
                call self%os(i)%kill
                call self%syms(i)%kill
                call self%smps(i)%kill
                call self%wimgs(i)%kill
            end do
            deallocate(self%recs, self%os, self%syms, self%smps, self%pouts, self%wimgs, self%bpls)
        endif
        if( allocated(self%rows) ) deallocate(self%rows, self%eo, self%vbin, self%eul, self%dfm)
        if( allocated(self%emhk) ) deallocate(self%emhk, self%emsh, self%emw, self%f2s, self%cw)
        if( allocated(self%jhk) ) deallocate(self%jhk, self%jsh, self%jw, self%f2j)
        if( allocated(self%fbhk) ) deallocate(self%fbhk, self%fbsh, self%fbw, self%f2b, self%cwb)
        self%nfb     = 0
        self%nshfb   = 0
        if( allocated(self%nv)  ) deallocate(self%nv)
        if( allocated(self%con) ) deallocate(self%con, self%cam, self%cse, self%cpj)
        if( allocated(self%cb)  ) deallocate(self%cb, self%cs)
        self%cfit = 0._dp
        self%nfit = 0
        if( allocated(self%msk3) ) deallocate(self%msk3, self%fsh, self%fw)
        if( allocated(self%fqh) ) deallocate(self%fqh, self%fqk)
        if( allocated(self%lmsk) ) deallocate(self%lmsk, self%lslot, self%msk0)
        if( allocated(self%cnl)  ) deallocate(self%cnl)
        self%l_equiv = .false.
        self%p => null()
        self%b => null()
        self%n = 0
        self%l_exists = .false.
        self%audit_calls = 0
    end subroutine kill

    integer function omp_thread()
        omp_thread = 1
        !$ omp_thread = omp_get_thread_num() + 1
    end function omp_thread

    subroutine pair_index( ii, k, g, q )
        integer, intent(in)  :: ii, k
        integer, intent(out) :: g, q
        integer :: c
        c = ii
        do g = 1, k-1
            if( c <= k - g )then
                q = g + c
                return
            endif
            c = c - (k - g)
        end do
        g = k - 1
        q = k
    end subroutine pair_index

    pure subroutine plane_cell( r, h, k, nyq, box, c1, c2, cnj )
        real,    intent(in)  :: r(3,3)
        integer, intent(in)  :: h, k, nyq, box
        integer, intent(out) :: c1(3), c2(3)
        logical, intent(out) :: cnj
        integer :: p(3)
        p   = nint(real(h) * r(1,:) + real(k) * r(2,:))
        cnj = p(1) < 0
        if( cnj ) p = -p
        c1  = 0
        c2  = 0
        if( any(abs(p) >= nyq) ) return
        c1 = [p(1) + 1, modulo(p(2), box) + 1, modulo(p(3), box) + 1]
        if( p(1) == 0 ) c2 = [1, modulo(-p(2), box) + 1, modulo(-p(3), box) + 1]
    end subroutine plane_cell

    pure subroutine zl_post( x, v, p, tq, lgq, r, et, et2, ea )
        real(dp), intent(in)  :: x, v, p(3), tq(:), lgq(:)
        real(dp), intent(out) :: r, et, et2, ea
        real(dp) :: lg(size(tq)), aq(size(tq)), q(size(tq)), l0, m, sq
        aq = exp(p(2) + p(3) * tq)
        lg = log(max(1._dp - p(1), 1.e-300_dp)) + lgq - 0.5_dp * (x - aq)**2 / v
        l0 = log(max(p(1), 1.e-300_dp)) - 0.5_dp * x * x / v
        m  = max(l0, maxval(lg))
        q  = exp(lg - m)
        sq = sum(q)
        r  = exp(l0 - m) / (exp(l0 - m) + sq)
        if( sq > 0._dp )then
            q   = q / sq
            et  = sum(q * tq)
            et2 = sum(q * tq * tq)
            ea  = sum(q * aq)
        else
            et  = 0._dp
            et2 = 1._dp
            ea  = exp(p(2) + 0.5_dp * p(3) * p(3))
        endif
    end subroutine zl_post

    function sshuf( lab, st ) result( out )
        use simple_rnd, only: shuffle
        integer, intent(in)  :: lab(:), st(:)
        integer, allocatable :: out(:), ix(:), sub(:)
        integer :: s, i
        out = lab
        do s = minval(st), maxval(st)
            ix = pack([(i, i=1,size(st))], st == s)
            if( size(ix) < 2 ) cycle
            sub = lab(ix)
            call shuffle(sub)
            out(ix) = sub
        end do
    end function sshuf

    pure integer function minnz( lab, k )
        integer, intent(in) :: lab(:), k
        integer :: g, c
        minnz = huge(minnz)
        do g = 1, k
            c = count(lab == g)
            if( c > 0 ) minnz = min(minnz, c)
        end do
    end function minnz

    pure function counts( lab, k ) result( c )
        integer, intent(in)  :: lab(:), k
        integer, allocatable :: c(:)
        integer :: g
        c = [(count(lab == g), g=1,k)]
    end function counts

    pure integer function nuniq( v )
        integer, intent(in) :: v(:)
        integer :: i
        nuniq = 0
        do i = 1, size(v)
            if( .not. any(v(1:i-1) == v(i)) ) nuniq = nuniq + 1
        end do
    end function nuniq

    subroutine isort( v )
        integer, intent(inout) :: v(:)
        integer :: i, j, t
        do i = 2, size(v)
            t = v(i)
            j = i - 1
            do while( j >= 1 )
                if( v(j) <= t ) exit
                v(j+1) = v(j)
                j = j - 1
            end do
            v(j+1) = t
        end do
    end subroutine isort

    ! the equivalence merge-up's inputs from a tree's own tables: every node whose split wrote K* >= 2 to split/split_stats.txt
    ! gives the effect r_hat of its weakest kept pair (the margin is the smallest); its stage-B leaves <node>_s1 / _s2 are
    ! pairs shown distinct; ntests = the largest split test count the nodes used. Leaf indices follow hetero_tree_leaves.txt
    subroutine equiv_tree_inputs( tdir, delta, dpairs, nrs, ntests )
        character(len=*),     intent(in)  :: tdir
        real,                 intent(out) :: delta
        integer, allocatable, intent(out) :: dpairs(:,:)
        integer,              intent(out) :: nrs, ntests
        character(len=256), allocatable :: lnames(:), nids(:)
        character(len=:),   allocatable :: statf
        integer :: i, a, b, k, kst, nt
        real    :: r, st(6)
        logical :: l_ex
        allocate(dpairs(2,0))
        delta  = 0.
        nrs    = 0
        ntests = 1
        call first_tokens(tdir//'/hetero_tree_leaves.txt', lnames)
        call first_tokens(tdir//'/hetero_tree_nodes.txt', nids)
        do i = 1, size(nids)
            statf = tdir//'/node_'//trim(nids(i))//'/split/split_stats.txt'
            inquire(file=statf, exist=l_ex)
            if( .not. l_ex ) cycle
            call read_split_stats(statf, st, kst, nt)
            ntests = max(ntests, nt)
            if( st(4) <= 0. .or. kst < 2 ) cycle
            r = split_effect(st(1), st(3), nint(st(6)))
            nrs = nrs + 1
            if( nrs == 1 )then
                delta = r
            else
                delta = min(delta, r)
            endif
            a = name_index(lnames, trim(nids(i))//'_s1')
            b = name_index(lnames, trim(nids(i))//'_s2')
            if( a > 0 .and. b > 0 ) dpairs = reshape([dpairs, a, b], [2, size(dpairs,2) + 1])
        end do
    end subroutine equiv_tree_inputs

    ! a split's effect per particle from its weakest kept pair statistic t at equal N n against the null mean m:
    ! y = (t - m) / (1 - m), r_hat = 2 y / ((1 - y) n)
    real function split_effect( t, m, n )
        real,    intent(in) :: t, m
        integer, intent(in) :: n
        real :: y
        y = min(max((t - m) / max(1. - m, 1.e-6), 1.e-4), 0.995)
        split_effect = 2. * y / ((1. - y) * real(max(n, 1)))
    end function split_effect

    ! split/split_stats.txt: one line of numbers written by state_split (weakest kept pair, weakest pair at the start, null mean,
    ! null sd, z_crit, equal N, K*, split tests)
    subroutine write_split_stats( fname, st, kst, ntests, kpre )
        character(len=*),  intent(in) :: fname
        real,              intent(in) :: st(6)
        integer,           intent(in) :: kst, ntests
        integer, optional, intent(in) :: kpre
        integer :: fu, io
        call fopen(fu, string(fname), status='replace', action='write', iostat=io)
        call fileiochk('state_split: cannot write '//fname, io)
        write(fu,'(A)') '# weakest kept pair statistic, weakest pair at the start, null mean, null sd, z_crit, equal N, K* (after the '//&
            &'guard and the floor), split tests, K* before them'
        if( present(kpre) )then
            write(fu,'(5(ES14.6,1X),I0,1X,I0,1X,I0,1X,I0)') st(1:5), nint(st(6)), kst, ntests, kpre
        else
            write(fu,'(5(ES14.6,1X),I0,1X,I0,1X,I0)') st(1:5), nint(st(6)), kst, ntests
        endif
        call fclose(fu)
    end subroutine write_split_stats

    ! per-particle quality terms for the node's reconstructions (the controller applies them): project row, 1 - p_empty, contrast
    subroutine write_weights( self, fname )
        class(state_mergeup), intent(inout) :: self
        character(len=*),     intent(in)    :: fname
        integer :: fu, io, i
        if( .not. allocated(self%cpj) ) call self%fit_contrast
        call fopen(fu, string(fname), status='replace', action='write', iostat=io)
        call fileiochk('state_split: cannot write '//fname, io)
        write(fu,'(A)') '# project row, 1 - zero-contrast posterior, posterior-mean contrast, log noise variance, stage-B corr'
        do i = 1, self%n
            write(fu,'(I0,1X,F8.5,1X,ES12.5,1X,F9.4,1X,F8.5)') self%rows(i), 1. - self%cpj(i), max(self%con(i), 1.e-6), &
                &log(max(self%nv(i), 1.e-12)), self%qcorr(i)
        end do
        call fclose(fu)
    end subroutine write_weights

    subroutine read_split_stats( fname, st, kst, ntests )
        character(len=*), intent(in)  :: fname
        real,             intent(out) :: st(6)
        integer,          intent(out) :: kst, ntests
        character(len=1024) :: ln
        integer :: fu, io, n6
        st     = 0.
        kst    = 0
        ntests = 1
        call fopen(fu, string(fname), status='old', action='read', iostat=io)
        if( io /= 0 ) return
        do
            read(fu,'(A)',iostat=io) ln
            if( io /= 0 ) exit
            if( len_trim(ln) == 0 .or. ln(1:1) == '#' ) cycle
            read(ln,*,iostat=io) st(1:5), n6, kst, ntests
            if( io == 0 ) st(6) = real(n6)
            exit
        end do
        call fclose(fu)
    end subroutine read_split_stats

    subroutine first_tokens( fname, toks )
        character(len=*),                intent(in)  :: fname
        character(len=256), allocatable, intent(out) :: toks(:)
        character(len=8192) :: line
        character(len=256)  :: tk
        integer :: funit, ios
        allocate(toks(0))
        open(newunit=funit, file=fname, status='old', action='read', iostat=ios)
        if( ios /= 0 ) return
        do
            read(funit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            if( len_trim(line) == 0 ) cycle
            read(line,*,iostat=ios) tk
            if( ios == 0 ) toks = [toks, tk]
        end do
        close(funit)
    end subroutine first_tokens

    integer function name_index( names, nm )
        character(len=*), intent(in) :: names(:), nm
        integer :: i
        name_index = 0
        do i = 1, size(names)
            if( trim(names(i)) == trim(nm) )then
                name_index = i
                return
            endif
        end do
    end function name_index

    function parent_dir( path ) result( d )
        character(len=*), intent(in)  :: path
        character(len=:), allocatable :: d
        integer :: p
        p = index(trim(path), '/', back=.true.)
        if( p > 1 )then
            d = path(1:p-1)
        else if( p == 1 )then
            d = '/'
        else
            d = '.'
        endif
    end function parent_dir

    function pair_list( pr ) result( s )
        integer, intent(in) :: pr(:,:)
        character(len=:), allocatable :: s
        integer :: k
        s = ''
        do k = 1, size(pr,2)
            s = s//' '//int2str(pr(1,k))//'-'//int2str(pr(2,k))
        end do
    end function pair_list

    real function medianr( v )
        use simple_stat, only: median_nocopy
        real, intent(in) :: v(:)
        real, allocatable :: w(:)
        medianr = 0.
        if( size(v) == 0 ) return
        w = v
        medianr = median_nocopy(w)
    end function medianr

    real function quantr( v, q )
        use simple_srch_sort_loc, only: hpsort
        real, intent(in) :: v(:), q
        real, allocatable :: w(:)
        real    :: pos
        integer :: lo
        quantr = 0.
        if( size(v) == 0 ) return
        w = v
        call hpsort(w)
        pos = real(size(w)-1) * q
        lo  = int(pos)
        quantr = w(lo+1) + (pos - real(lo)) * (w(min(lo+2,size(w))) - w(lo+1))
    end function quantr

    subroutine null_thr( nl, mu, sd, thr )
        type(shuffle_null), intent(in)  :: nl
        real,               intent(out) :: mu, sd, thr
        mu  = sum(nl%nul2) / real(size(nl%nul2))
        sd  = sqrt(sum((nl%nul2 - mu)**2) / real(max(1, size(nl%nul2) - 1)))
        thr = maxval(nl%nul2)
        if( nl%l_z ) thr = mu + nl%zc * sd
    end subroutine null_thr

    real function null_spread( n, nmax, zc )
        integer, intent(in) :: n, nmax
        real,    intent(in) :: zc
        real :: m
        null_spread = 0.
        if( n >= nmax ) return
        m = real(nmax - n)
        null_spread = sqrt(m / real(nmax * n) + zc * zc * (m - 1.) * real(nmax - 2) / &
            &(2. * real(nmax - 1)**2 * real(max(1, n - 1))))
    end function null_spread

    ! the split re-test's Bonferroni threshold over m split tests: the one-sided normal quantile at ALPHA_SPLIT / m
    real function split_zcrit( m )
        integer, intent(in) :: m
        real(dp) :: a, lo, hi, mid
        integer  :: it
        a  = ALPHA_SPLIT / real(max(1, m), dp)
        lo = 0._dp
        hi = 12._dp
        do it = 1, 80
            mid = 0.5_dp * (lo + hi)
            if( 0.5_dp * erfc(mid / sqrt(2._dp)) > a )then
                lo = mid
            else
                hi = mid
            endif
        end do
        split_zcrit = real(0.5_dp * (lo + hi))
    end function split_zcrit

    function sorted_upper( r, n ) result( str )
        use simple_srch_sort_loc, only: hpsort
        real,    intent(in) :: r(:,:)
        integer, intent(in) :: n
        character(len=:), allocatable :: str
        real, allocatable :: v(:)
        integer :: p, q
        allocate(v(0))
        do p = 1, n
            do q = p+1, n
                v = [v, r(p,q)]
            end do
        end do
        if( size(v) > 1 ) call hpsort(v)
        str = joinedr(v)
    end function sorted_upper

    function joined( vals, sep ) result( str )
        integer,                    intent(in) :: vals(:)
        character(len=*), optional, intent(in) :: sep
        character(len=:), allocatable :: str
        integer :: i
        str = ''
        do i = 1, size(vals)
            if( i > 1 )then
                if( present(sep) )then
                    str = str//sep
                else
                    str = str//' '
                endif
            endif
            str = str//int2str(vals(i))
        end do
    end function joined

    function joinedr( vals ) result( str )
        real, intent(in) :: vals(:)
        character(len=:), allocatable :: str
        integer :: i
        str = ''
        do i = 1, size(vals)
            if( i > 1 ) str = str//' '
            str = str//fmtr(real(vals(i),dp),2)
        end do
    end function joinedr

    function fmtr( x, ndec ) result( str )
        real(dp), intent(in)          :: x
        integer,  intent(in)          :: ndec
        character(len=:), allocatable :: str
        character(len=40) :: buf
        character(len=16) :: fmt
        write(fmt,'(A,I0,A)') '(F40.', ndec, ')'
        write(buf,fmt) x
        str = trim(adjustl(buf))
    end function fmtr

    ! diffusely differing split groups merge; pairs no stage-B state pair spans are judged on the children's own half maps (mu)
    subroutine locality_guard( lab, child, kstar, mskdiam, mu )
        integer, intent(in)    :: lab(:)
        integer, intent(inout) :: child(:), kstar
        real,    intent(in)    :: mskdiam
        class(state_mergeup), optional, intent(inout) :: mu
        real, parameter :: LOC_MIN = 0.5
        character(len=*), parameter :: CHD = 'locality_children'
        character(len=:), allocatable :: stb
        real,    allocatable :: rat(:,:)
        integer, allocatable :: smap(:), cnt(:), sz(:), rk(:)
        logical, allocatable :: have(:)
        real    :: r, rmin
        integer :: nst, ng, s, t, g, h, gm, hm, kold
        logical :: l_ex, l_kids
        if( kstar < 2 ) return
        stb = '../stageB'
        nst = maxval(lab)
        ng  = maxval(child)
        allocate(have(nst), source=.false.)
        allocate(smap(nst), source=0)
        allocate(cnt(ng))
        do s = 1, nst
            if( count(lab == s .and. child > 0) == 0 ) cycle
            cnt     = [(count(lab == s .and. child == g), g=1,ng)]
            smap(s) = maxloc(cnt, dim=1)
            inquire(file=stb//'/recvol_state'//int2str_pad(s,2)//'_even_unfil.mrc', exist=have(s))
            inquire(file=stb//'/recvol_state'//int2str_pad(s,2)//'_odd_unfil.mrc',  exist=l_ex)
            have(s) = have(s) .and. l_ex
        end do
        if( count(have) < 2 .and. .not. present(mu) )then
            write(logfhandle,'(A)') 'locality guard: stage-B half maps not found under '//stb//'; skipped'
            return
        endif
        if( count(have) < 2 ) write(logfhandle,'(A)') 'locality guard: stage-B half maps not found under '//stb//&
            &'; every group pair is judged on the delivered children''s half maps'
        allocate(rat(nst,nst), source=-1.)
        do s = 1, nst - 1
            do t = s + 1, nst
                if( .not. (have(s) .and. have(t)) .or. smap(s) == smap(t) ) cycle
                rat(s,t) = pair_locality(stb, s, t, mskdiam)
                rat(t,s) = rat(s,t)
                write(logfhandle,'(A)') 'locality guard: stage-B states '//int2str(s)//' / '//int2str(t)//' L / L_ref '//&
                    &fmtr(real(rat(s,t),dp),3)
            end do
        end do
        kold   = kstar
        l_kids = .false.
        do
            if( count([(any(child == g), g=1,ng)]) < 2 ) exit
            rmin = huge(rmin)
            gm   = 0
            hm   = 0
            do g = 1, ng - 1
                do h = g + 1, ng
                    if( .not. (any(child == g) .and. any(child == h)) ) cycle
                    r = -1.
                    do s = 1, nst
                        do t = 1, nst
                            if( smap(s) == g .and. smap(t) == h ) r = max(r, rat(s,t))
                        end do
                    end do
                    if( r < 0. )then
                        if( .not. present(mu) )then
                            write(logfhandle,'(A)') 'locality guard: groups '//int2str(g)//' / '//int2str(h)//&
                                &': no stage-B state pair spans them and no children''s maps; kept unjudged'
                            cycle
                        endif
                        if( .not. l_kids ) call mu%group_halves(child, CHD, stb)
                        l_kids = .true.
                        r = pair_locality(CHD, g, h, mskdiam)
                        write(logfhandle,'(A)') 'locality guard: groups '//int2str(g)//' / '//int2str(h)//&
                            &': no stage-B state pair spans them; L / L_ref from the delivered children''s half maps '//&
                            &fmtr(real(r,dp),3)
                    endif
                    if( r < rmin )then
                        rmin = r
                        gm   = g
                        hm   = h
                    endif
                end do
            end do
            if( rmin >= LOC_MIN ) exit
            write(logfhandle,'(A)') 'locality guard: groups '//int2str(gm)//' / '//int2str(hm)//' differ diffusely (L / L_ref '//&
                &fmtr(real(rmin,dp),3)//' < '//fmtr(real(LOC_MIN,dp),2)//') -> merged'
            where( child == hm ) child = gm
            where( smap  == hm ) smap  = gm
            l_kids = .false.
        end do
        sz = [(count(child == g), g=1,ng)]
        allocate(rk(0:ng), source=0)
        do g = 1, ng
            if( sz(g) > 0 ) rk(g) = count(sz > sz(g)) + count(sz(1:g-1) == sz(g)) + 1
        end do
        child = rk(child)
        kstar = count(sz > 0)
        write(logfhandle,'(A)') 'locality guard: K* '//int2str(kold)//' -> '//int2str(kstar)
    end subroutine locality_guard

    subroutine group_halves( self, child, dir, stb )
        use simple_imghead, only: get_mrcfile_info
        use simple_syslib,  only: dir_exists, simple_mkdir
        class(state_mergeup), intent(inout) :: self
        integer,              intent(in)    :: child(:)
        character(len=*),     intent(in)    :: dir, stb
        type(image)          :: vt
        complex, allocatable :: cm(:,:,:)
        real,    allocatable :: e(:)
        integer, allocatable :: sel(:), js(:)
        character(len=:), allocatable :: fn
        integer :: g, h, j, s, ldim(3), boxt
        real    :: smpdt
        logical :: l_ex
        ! LOC_MIN was set on stage-B halves, so the children's maps go on their grid (box LOC_BOX when none is on disk)
        boxt  = min(self%p%box, LOC_BOX)
        smpdt = real(self%p%box) * self%p%smpd / real(boxt)
        do s = 1, 99
            fn = stb//'/recvol_state'//int2str_pad(s,2)//'_even_unfil.mrc'
            inquire(file=fn, exist=l_ex)
            if( .not. l_ex ) cycle
            call get_mrcfile_info(string(fn), ldim, 'M', smpdt, .false.)
            boxt = ldim(1)
            exit
        end do
        if( boxt < self%box )then
            smpdt = smpdt * real(boxt) / real(self%box)
            boxt  = self%box
        endif
        if( .not. dir_exists(string(dir)) ) call simple_mkdir(string(dir))
        allocate(cm(self%nh,self%box,self%box))
        do g = 1, maxval(child)
            if( count(child == g) == 0 ) cycle
            do h = 0, 1
                sel = pack([(j, j=1,self%n)], child == g .and. self%eo == h)
                if( size(sel) == 0 ) THROW_HARD('a delivered group has no particles in one half; group_halves')
                call self%accumulate_jobs(sel, [1, size(sel)+1], .false., js)
                call self%recs(js(1))%get_rho_shell_medians(e)
                call self%finish(js(1), EPS_REG * e, .false., cm)
                call vt%new([boxt,boxt,boxt], smpdt)
                call self%wimgs(js(1))%pad(vt, antialiasing=.false.)
                call vt%ifft
                call vt%write(string(dir//'/recvol_state'//int2str_pad(g,2)//'_'//trim(merge('even', 'odd ', h == 0))//'_unfil.mrc'))
                call vt%kill
            end do
        end do
        write(logfhandle,'(A)') 'locality guard: delivered children''s half maps written to '//dir//' (box '//int2str(boxt)//&
            &', '//fmtr(real(smpdt,dp),3)//' A)'
    end subroutine group_halves

    real function pair_locality( stb, s, t, mskdiam ) result( ratio )
        use simple_imghead, only: get_mrcfile_info
        character(len=*), intent(in) :: stb
        integer,          intent(in) :: s, t
        real,             intent(in) :: mskdiam
        type(image) :: img(4)
        character(len=:), allocatable :: fn(:)
        real    :: smpd, p
        integer :: ldim(3), i
        allocate(character(len=len(stb)+40) :: fn(4))
        fn(1) = stb//'/recvol_state'//int2str_pad(s,2)//'_even_unfil.mrc'
        fn(2) = stb//'/recvol_state'//int2str_pad(s,2)//'_odd_unfil.mrc'
        fn(3) = stb//'/recvol_state'//int2str_pad(t,2)//'_even_unfil.mrc'
        fn(4) = stb//'/recvol_state'//int2str_pad(t,2)//'_odd_unfil.mrc'
        call get_mrcfile_info(string(trim(fn(1))), ldim, 'M', smpd, .false.)
        do i = 1, 4
            call img(i)%new(ldim, smpd)
            call img(i)%read(string(trim(fn(i))))
        end do
        call pair_maps(img, mskdiam, ratio, p)
        do i = 1, 4
            call img(i)%kill
        end do
    end function pair_locality

    ! the effect floor: a child pair whose reproducible difference power is under EFF_MIN x the node's half-map noise power merges
    subroutine effect_floor( lab, child, kstar, mskdiam )
        use simple_imghead, only: get_mrcfile_info
        integer, intent(in)    :: lab(:)
        integer, intent(inout) :: child(:), kstar
        real,    intent(in)    :: mskdiam
        real, parameter :: EFF_MIN = 2.
        character(len=*), parameter :: STB = '../stageB'
        type(image)       :: img(4)
        real, allocatable :: em(:,:,:,:), om(:,:,:,:), eu(:,:,:), ou(:,:,:), ge(:,:,:), go(:,:,:), he(:,:,:), ho(:,:,:)
        integer, allocatable :: cnt(:,:), sz(:), rk(:)
        real    :: smpd, nnode, p, r, rmin, dum
        integer :: ldim(3), nst, ng, s, g, h, gm, hm, kold, i
        logical :: l_ex
        if( kstar < 2 ) return
        nst = maxval(lab)
        ng  = maxval(child)
        allocate(cnt(nst,ng))
        do s = 1, nst
            cnt(s,:) = [(count(lab == s .and. child == g), g=1,ng)]
        end do
        do s = 1, nst
            if( sum(cnt(s,:)) == 0 ) cycle
            inquire(file=STB//'/recvol_state'//int2str_pad(s,2)//'_even_unfil.mrc', exist=l_ex)
            if( .not. l_ex )then
                write(logfhandle,'(A)') 'effect floor: stage-B half maps not found under '//STB//'; skipped'
                return
            endif
            call get_mrcfile_info(string(STB//'/recvol_state'//int2str_pad(s,2)//'_even_unfil.mrc'), ldim, 'M', smpd, .false.)
        end do
        allocate(em(ldim(1),ldim(2),ldim(3),nst), om(ldim(1),ldim(2),ldim(3),nst), source=0.)
        allocate(eu(ldim(1),ldim(2),ldim(3)), ou(ldim(1),ldim(2),ldim(3)), source=0.)
        do i = 1, 4
            call img(i)%new(ldim, smpd)
        end do
        do s = 1, nst
            if( sum(cnt(s,:)) == 0 ) cycle
            call img(1)%read(string(STB//'/recvol_state'//int2str_pad(s,2)//'_even_unfil.mrc'))
            call img(2)%read(string(STB//'/recvol_state'//int2str_pad(s,2)//'_odd_unfil.mrc'))
            em(:,:,:,s) = img(1)%get_rmat()
            om(:,:,:,s) = img(2)%get_rmat()
            eu = eu + real(sum(cnt(s,:))) * em(:,:,:,s)
            ou = ou + real(sum(cnt(s,:))) * om(:,:,:,s)
        end do
        eu = eu / real(max(sum(cnt), 1))
        ou = ou / real(max(sum(cnt), 1))
        call img(1)%set_rmat(eu, .false.)
        call img(2)%set_rmat(ou, .false.)
        call img(3)%set_rmat(eu, .false.)
        call img(4)%set_rmat(ou, .false.)
        call pair_maps(img, mskdiam, dum, p, nnode)
        kold = kstar
        do
            if( count([(any(child == g), g=1,ng)]) < 2 ) exit
            rmin = huge(rmin)
            gm   = 0
            hm   = 0
            do g = 1, ng - 1
                if( .not. any(child == g) ) cycle
                call group_maps(g, ge, go)
                do h = g + 1, ng
                    if( .not. any(child == h) ) cycle
                    call group_maps(h, he, ho)
                    call img(1)%set_rmat(ge, .false.)
                    call img(2)%set_rmat(go, .false.)
                    call img(3)%set_rmat(he, .false.)
                    call img(4)%set_rmat(ho, .false.)
                    call pair_maps(img, mskdiam, dum, p)
                    r = p / max(nnode, 1.e-30)
                    write(logfhandle,'(A)') 'effect floor: groups '//int2str(g)//' / '//int2str(h)//' difference power over the '//&
                        &'node''s half-map noise power P / N '//fmtr(real(r,dp),3)
                    if( r < rmin )then
                        rmin = r
                        gm   = g
                        hm   = h
                    endif
                end do
            end do
            if( rmin >= EFF_MIN ) exit
            write(logfhandle,'(A)') 'effect floor: groups '//int2str(gm)//' / '//int2str(hm)//' differ by too little (P / N '//&
                &fmtr(real(rmin,dp),3)//' < '//fmtr(real(EFF_MIN,dp),2)//') -> merged'
            where( child == hm ) child = gm
            cnt(:,gm) = cnt(:,gm) + cnt(:,hm)
            cnt(:,hm) = 0
        end do
        sz = [(count(child == g), g=1,ng)]
        allocate(rk(0:ng), source=0)
        do g = 1, ng
            if( sz(g) > 0 ) rk(g) = count(sz > sz(g)) + count(sz(1:g-1) == sz(g)) + 1
        end do
        child = rk(child)
        kstar = count(sz > 0)
        write(logfhandle,'(A)') 'effect floor: K* '//int2str(kold)//' -> '//int2str(kstar)
        do i = 1, 4
            call img(i)%kill
        end do

    contains

        subroutine group_maps( g, e, o )
            integer,           intent(in)  :: g
            real, allocatable, intent(out) :: e(:,:,:), o(:,:,:)
            integer :: s
            allocate(e(ldim(1),ldim(2),ldim(3)), o(ldim(1),ldim(2),ldim(3)), source=0.)
            do s = 1, nst
                if( cnt(s,g) == 0 ) cycle
                e = e + real(cnt(s,g)) * em(:,:,:,s)
                o = o + real(cnt(s,g)) * om(:,:,:,s)
            end do
            e = e / real(max(sum(cnt(:,g)), 1))
            o = o / real(max(sum(cnt(:,g)), 1))
        end subroutine group_maps

    end subroutine effect_floor

    ! the guard's pair statistic on four half maps (s: img 1, 2; t: img 3, 4): L / L_ref, the difference power P, s's noise power
    subroutine pair_maps( img, mskdiam, ratio, p, pnoise )
        use simple_srch_sort_loc, only: hpsort
        type(image),    intent(inout) :: img(4)
        real,           intent(in)    :: mskdiam
        real,           intent(out)   :: ratio, p
        real, optional, intent(out)   :: pnoise
        real, parameter :: SIG = 3., TOPQ = 0.05, EDGE = 6.
        type(image)       :: w
        real, allocatable :: pw(:,:), fb(:), c14(:), c23(:), g(:), gk(:), bf(:), m(:,:,:), sp(:)
        real, allocatable :: r1(:,:,:), r2(:,:,:), r3(:,:,:), r4(:,:,:), res(:,:,:), nse(:,:,:), con(:,:,:), ex(:,:,:)
        real    :: smpd, rad, d, den, lsh, lref
        integer :: ldim(3), i, j, k, n, sh
        ldim = img(1)%get_ldim()
        smpd = img(1)%get_smpd()
        rad  = 0.5 * mskdiam / smpd
        allocate(m(ldim(1),ldim(2),ldim(3)))
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    d = sqrt(real((i - ldim(1)/2 - 1)**2 + (j - ldim(2)/2 - 1)**2 + (k - ldim(3)/2 - 1)**2))
                    m(i,j,k) = 0.
                    if( d <= rad + EDGE ) m(i,j,k) = 0.5 * (1. + cos(PI * min(max(d - rad, 0.), EDGE) / EDGE))
                end do
            end do
        end do
        do i = 1, 4
            call img(i)%set_rmat(img(i)%get_rmat() * m, .false.)
            call img(i)%fft()
        end do
        n = img(1)%get_filtsz()
        allocate(pw(n,4), source=0.)
        allocate(fb(n), c14(n), c23(n), g(n), gk(n), bf(n))
        do i = 1, 4
            call img(i)%spectrum('power', sp, norm=.false.)
            pw(1:min(n,size(sp)),i) = sp(1:min(n,size(sp)))
        end do
        call img(3)%fsc(img(4), fb)
        call img(1)%fsc(img(4), c14)
        call img(2)%fsc(img(3), c23)
        do sh = 1, n
            den   = fb(sh) * sqrt(pw(sh,3) * pw(sh,4))
            g(sh) = 0.
            if( den > 0. ) g(sh) = 0.5 * (c14(sh) * sqrt(pw(sh,1) * pw(sh,4)) + c23(sh) * sqrt(pw(sh,2) * pw(sh,3))) / den
            gk(sh) = exp(-2. * (PI * SIG * real(sh) / real(ldim(1)))**2)
            d      = real(ldim(1)) * smpd / real(sh)
            bf(sh) = merge(1., 0., d <= 40. .and. d > 10.)
        end do
        do i = 1, 4
            call img(i)%ifft()
        end do
        r1 = img(1)%get_rmat()
        r2 = img(2)%get_rmat()
        r3 = img(3)%get_rmat()
        r4 = img(4)%get_rmat()
        call w%new(ldim, smpd)
        res = 0.5 * (r1 + r2) - filtered(0.5 * (r3 + r4))
        nse = 0.5 * (r1 - r2) - filtered(0.5 * (r3 - r4))
        res = banded(res)
        nse = banded(nse)
        con = banded(0.25 * (r1 + r2 + r3 + r4))
        ex  = max(smoothed(res * res) - smoothed(nse * nse), 0.)
        lsh  = top_share(pack(ex, m > 0.5))
        lref = top_share(pack(smoothed(con * con), m > 0.5))
        ratio = 0.
        if( lref > 0. ) ratio = lsh / lref
        p = sum(pack(ex, m > 0.5))
        if( present(pnoise) )then
            nse    = banded(0.5 * (r1 - r2))
            pnoise = sum(pack(smoothed(nse * nse), m > 0.5))
        endif
        do i = 1, 4
            call img(i)%ifft()
        end do
        call w%kill

    contains

        function filtered( x ) result( y )
            real, intent(in)  :: x(:,:,:)
            real, allocatable :: y(:,:,:)
            call w%set_rmat(x, .false.)
            call w%apply_filter(g)
            y = w%get_rmat()
        end function filtered

        function banded( x ) result( y )
            real, intent(in)  :: x(:,:,:)
            real, allocatable :: y(:,:,:)
            call w%set_rmat(x, .false.)
            call w%apply_filter(bf)
            y = w%get_rmat()
            y = y - sum(y) / real(size(y))
        end function banded

        function smoothed( x ) result( y )
            real, intent(in)  :: x(:,:,:)
            real, allocatable :: y(:,:,:)
            call w%set_rmat(x, .false.)
            call w%apply_filter(gk)
            y = w%get_rmat()
        end function smoothed

        real function top_share( v )
            real, intent(in)  :: v(:)
            real, allocatable :: vs(:)
            integer :: nt
            top_share = 0.
            if( size(v) == 0 ) return
            if( sum(v) <= 0. ) return
            vs = v
            call hpsort(vs)
            nt = max(1, int(TOPQ * real(size(vs))))
            top_share = sum(vs(size(vs)-nt+1:)) / sum(vs)
        end function top_share

    end subroutine pair_maps

end module simple_state_mergeup
