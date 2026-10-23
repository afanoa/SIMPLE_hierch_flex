!@descr: cross-half split test of the states of a node (SSNR junk screen, docked equal-N pair maps against a stratified shuffle null, cross-half EM, re-test beyond the EM band) and the merge-up of states along a state tree
module simple_state_mergeup
!$ use omp_lib
use simple_defs,          only: dp, logfhandle, OSMPL_PAD_FAC, MAXIMGBATCHSZ, PI, nthr_glob
use simple_type_defs,     only: fplane_type
use simple_error,         only: simple_exception
use simple_string,        only: string
use simple_string_utils,  only: int2str
use simple_fileio,        only: fopen, fclose, fileiochk, nlines
use simple_ori,           only: ori
use simple_sym,           only: sym
use simple_image,         only: image
use simple_reconstructor, only: reconstructor, exp_samples
use simple_parameters,    only: parameters
use simple_builder,       only: builder
use simple_state_tree,    only: state_tree
use simple_stat,          only: avg_ranks
implicit none

public :: state_mergeup, merge_poles, merge_poles_halves, halves_merge_margin, halves_accept_margin, halves_se, EM_LIM
public :: HALVES_KSIG, halves_floor, halves_cand, mask_neff, seed_poles
private
#include "simple_local_flags.inc"

integer, parameter :: NSHUF = 20, NVIEW = 40, NEM_ITERS = 3, EM_MIN = 1000, MINSEL_JUNK = 50, MAXJOBS = 1000, RND_SEED = 0
integer, parameter :: NDOCK_ITS = 8, RB_MARGIN = 2, RB_MARGIN_FOCUS = 4, SITE_ITS = 50, SITE_ITS_MIN = 10
real,    parameter :: SITE_CONV = 0.005, HALVES_KSIG = 3.
integer, parameter :: SITE_GRID(12) = [500, 1000, 1500, 2500, 4000, 6000, 10000, 15000, 25000, 40000, 60000, 90000]
real,    parameter :: SITE_PFLOOR = 1.e-3
real,    parameter :: EPS_REG = 0.05, EM_LIM = 20., JUNK_CORR = 0.8, JUNK_SSNR = 0.35, JUNK_PROJ = 0.75, DOCK_MAX = 3.
real,    parameter :: JUNK_LIM = 12., JUNK_GAP = 2., JUNK_MINW = 0.005
integer, parameter :: JUNK_ITS = 500

type :: ivec
    integer, allocatable :: v(:)
end type ivec

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
    real,    allocatable :: eul(:,:), dfm(:), emw(:), f2s(:,:), cw(:), msk3(:,:,:), fw(:,:,:), msk0(:,:,:), lmsk(:,:,:,:)
    real(dp), allocatable :: fqh(:), fqk(:)
    integer, allocatable :: fsh(:,:,:)
    integer :: n = 0, nslots = 1, box = 0, nh = 0, nyq = 0, klo = 0, khi = 0, kb = 0, kt = 0, nem = 0, nshem = 0
    logical :: l_exists = .false.
    integer, allocatable :: jhk(:,:), jsh(:)
    real,    allocatable :: jw(:), f2j(:,:)
    integer :: nj = 0, kj = 0
  contains
    procedure :: new
    procedure :: focus
    procedure :: leaf_focus
    procedure :: run
    procedure :: test
    procedure :: site_em
    procedure :: junk_screen
    procedure :: get_rows
    procedure :: holdout_band
    procedure :: kill
    procedure, private :: consensus
    procedure, private :: screen
    procedure, private :: decide
    procedure, private :: eqn_stats
    procedure, private :: cross_em
    procedure, private :: estep
    procedure, private :: maps
    procedure, private :: accumulate
    procedure, private :: finish
    procedure, private :: dock
    procedure, private :: strata
    procedure, private :: pair_stat
    procedure, private :: shell_fsc
    procedure, private :: shell_corr
    procedure, private :: shell_dot
    procedure          :: half_rel
    procedure          :: half_corr
end type state_mergeup

contains

    subroutine new( self, params, build, em_lim )
        use simple_matcher_3Drec,   only: init_rec, prep_imgs4rec, cleanup_rec_buffers
        use simple_matcher_ptcl_io, only: prepimgbatch, ptcl_read_session
        use simple_imgarr_utils,    only: alloc_imgarr, dealloc_imgarr
        use simple_rnd,             only: seed_rnd_fixed
        class(state_mergeup),      intent(inout) :: self
        class(parameters), target, intent(inout) :: params
        class(builder),    target, intent(inout) :: build
        real,                      intent(in)    :: em_lim
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
        allocate(self%eo(self%n), self%vbin(self%n), self%eul(3,self%n), self%dfm(self%n))
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
                &crop_imgs=crop_imgs, kfromto=[1,self%nyq], observation_model=.true., whiten=.false.)
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
            self%fpls(i)%cmplx_plane = real(self%fpls(i)%transfer_plane) * self%fpls(i)%cmplx_plane
            deallocate(self%fpls(i)%transfer_plane)
        end do
        !$omp end parallel do
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
        integer,    allocatable :: grp(:), names(:), l2u(:), uptcl(:), cand(:), child(:), chrow(:), newu(:)
        integer,    allocatable :: to(:), cnt(:), merged(:), junk(:), tst(:), pop(:)
        logical,    allocatable :: isj(:), newj(:)
        integer :: g, lv, gi, il, u, nid, i, j, kstar, c, nmerge, ntest, nbefore, big
        units = max(labels, 0)
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
                call self%test(cand, size(tst), child, kstar, junk)
                if( allocated(self%lslot) ) self%msk3 = self%msk0
                do j = 1, size(junk)
                    call pass(tst(junk(j)), .true.)
                end do
                allocate(chrow(size(units)), source=0)
                chrow(self%rows) = child
                allocate(to(size(tst)), source=0)
                allocate(pop(max(1,kstar)), source=0)
                do c = 1, kstar
                    pop(c) = count(chrow == c)
                end do
                do j = 1, size(tst)
                    if( any(junk == j) )then
                        to(j) = -1
                        cycle
                    endif
                    u = tst(j)
                    allocate(cnt(max(1,kstar)), source=0)
                    do i = 1, size(units)
                        if( units(i) == u .and. chrow(i) > 0 ) cnt(chrow(i)) = cnt(chrow(i)) + 1
                    end do
                    if( sum(cnt) > 0 ) to(j) = maxloc(cnt, dim=1)
                    deallocate(cnt)
                end do
                big = 1
                if( any(to > 0) ) big = maxloc(pop, dim=1, mask=[(any(to == c), c=1,size(pop))])
                where( to == 0 ) to = big
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
                deallocate(cand, chrow, to, pop, child, junk)
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

    subroutine test( self, cand, k0, child, kstar, junk )
        class(state_mergeup), intent(inout) :: self
        integer,              intent(in)    :: cand(:), k0
        integer, allocatable, intent(out)   :: child(:), junk(:)
        integer,              intent(out)   :: kstar
        complex, allocatable :: cons(:,:,:,:), nm(:,:,:,:)
        integer, allocatable :: ia(:), ib(:), la(:), lb(:), sta(:), stb(:), pa(:), pb(:), jidx(:), joff(:), gcom(:)
        integer, allocatable :: plist(:), gof(:), exa(:), exb(:), starta(:), startb(:), sizes(:), order(:), rank(:)
        real,    allocatable :: eps(:), nul(:), nul2(:), jscale(:), vals(:), vals2(:)
        real    :: thr, thr2
        integer :: neq, s, s1, s2, nbatch, g, q, m, npair, ii, jj, nb, npr, ipr, kem, ng
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
        write(logfhandle,'(A)') 'equal N '//int2str(neq)//' per group; band '//int2str(self%klo)//'..'//int2str(self%khi)//&
            &', re-test band '//int2str(self%kt)//'..'//int2str(self%khi)//'; '//int2str(NSHUF)//&
            &' shuffles; stratified null: '//int2str(nuniq(sta))//' / '//int2str(nuniq(stb))//&
            &' nuisance strata (view x defocus)'
        allocate(nul(0), nul2(0))
        if( ng >= 2 )then
            nbatch = max(1, min(NSHUF, MAXJOBS / max(1, 2*ng)))
            npr    = ng * (ng - 1) / 2
            do s1 = 1, NSHUF, nbatch
                s2 = min(NSHUF, s1 + nbatch - 1)
                nb = s2 - s1 + 1
                allocate(jidx(0), joff(1), jscale(0))
                joff(1) = 1
                do s = s1, s2
                    pa = sshuf(la, sta)
                    pb = sshuf(lb, stb)
                    call add_subset_jobs(ia, pa)
                    call add_subset_jobs(ib, pb)
                end do
                allocate(nm(self%nh,self%box,self%box,2*ng*nb))
                call self%maps(jidx, joff, jscale, eps, .true., .false., nm)
                !$omp parallel do default(shared) private(m) schedule(dynamic) proc_bind(close)
                do m = 1, 2*ng*nb
                    call self%dock(nm(:,:,:,m), cons(:,:,:,2-mod((m-1)/ng,2)), self%kt, self%khi)
                end do
                !$omp end parallel do
                allocate(vals(npr*nb), vals2(npr*nb))
                !$omp parallel do default(shared) private(ipr,s,g,q,ii) schedule(dynamic) proc_bind(close)
                do ipr = 1, npr*nb
                    s  = (ipr - 1) / npr
                    ii = ipr - s * npr
                    call pair_index(ii, ng, g, q)
                    vals(ipr)  = self%pair_stat(nm(:,:,:,s*2*ng+g), nm(:,:,:,s*2*ng+q), &
                        &nm(:,:,:,s*2*ng+ng+g), nm(:,:,:,s*2*ng+ng+q), self%klo)
                    vals2(ipr) = self%pair_stat(nm(:,:,:,s*2*ng+g), nm(:,:,:,s*2*ng+q), &
                        &nm(:,:,:,s*2*ng+ng+g), nm(:,:,:,s*2*ng+ng+q), self%kt)
                end do
                !$omp end parallel do
                nul  = [nul, vals]
                nul2 = [nul2, vals2]
                deallocate(nm, jidx, joff, jscale, vals, vals2)
            end do
        endif
        thr  = huge(thr)
        thr2 = huge(thr2)
        if( size(nul) > 0 )then
            thr  = maxval(nul)
            thr2 = maxval(nul2)
            write(logfhandle,'(A)') 'shuffle null ('//int2str(size(nul))//' pair statistics): median '//&
                &fmtr(real(medianr(nul),dp),3)//', 95 % '//fmtr(real(quantr(nul,0.95),dp),3)//', max '//&
                &fmtr(real(thr,dp),3)//' -> threshold '//fmtr(real(thr,dp),3)
            write(logfhandle,'(A)') 'beyond-EM band '//int2str(self%kt)//'..'//int2str(self%khi)//': shuffle null median '//&
                &fmtr(real(medianr(nul2),dp),3)//', max '//fmtr(real(thr2,dp),3)//' -> threshold '//fmtr(real(thr2,dp),3)
        endif
        call self%decide(ia, ib, la, lb, k0, eps, neq, cons, self%klo, thr, 'board groups', plist, gof, exa, exb, kstar)
        if( kstar >= 2 )then
            allocate(starta(size(ia)), startb(size(ib)), source=0)
            do m = 1, size(plist)
                where( la == plist(m) ) starta = gof(m)
                where( lb == plist(m) ) startb = gof(m)
            end do
            kem = kstar
            call self%cross_em(ia, ib, starta, startb, kem, eps, neq, la, lb)
            write(logfhandle,'(A)') 'cross-half EM half A: MAP changed vs start '//fmtr(real(count(la /= starta),dp)/&
                &real(size(la),dp),3)//'; half B '//fmtr(real(count(lb /= startb),dp)/real(size(lb),dp),3)
            call self%decide(ia, ib, la, lb, kem, eps, neq, cons, self%kt, thr2, 'cross-half EM groups', plist, gof, &
                &exa, exb, kstar)
        endif
        npair = size(plist)
        allocate(sizes(kstar), order(kstar), rank(kstar), source=0)
        do m = 1, npair
            sizes(gof(m)) = sizes(gof(m)) + count(la == plist(m)) + count(lb == plist(m))
        end do
        order = [(ii, ii=1,kstar)]
        do ii = 2, kstar
            jj = ii
            do while( jj > 1 )
                if( sizes(order(jj-1)) >= sizes(order(jj)) ) exit
                order(jj-1:jj) = [order(jj), order(jj-1)]
                jj = jj - 1
            end do
        end do
        do ii = 1, kstar
            rank(order(ii)) = ii
        end do
        allocate(child(self%n), source=0)
        do ii = 1, size(ia)
            child(ia(ii)) = group_of(la(ii), exa)
        end do
        do ii = 1, size(ib)
            child(ib(ii)) = group_of(lb(ii), exb)
        end do
        call junk_report(ia, 1, 'A')
        call junk_report(ib, 2, 'B')
        write(logfhandle,'(A)') 'K* = '//int2str(kstar)//'; sizes '//joined([(count(child == ii), ii=1,kstar)])//&
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

        subroutine add_subset_jobs( idx, lab )
            integer, intent(in) :: idx(:), lab(:)
            integer, allocatable :: sel(:)
            integer :: gg
            do gg = 1, ng
                sel = pack(idx, lab == gcom(gg))
                if( size(sel) < neq ) THROW_HARD('a shuffled group is smaller than the equal N')
                call rnd_subset(sel, neq)
                jidx   = [jidx, sel(1:neq)]
                joff   = [joff, size(jidx) + 1]
                jscale = [jscale, 1.]
            end do
        end subroutine add_subset_jobs

        subroutine junk_report( idx, hh, side )
            integer,          intent(in) :: idx(:), hh
            character(len=*), intent(in) :: side
            complex, allocatable :: jm(:,:,:,:)
            integer, allocatable :: sel(:), jo(:), ji(:)
            real,    allocatable :: js(:)
            character(len=:), allocatable :: line
            real    :: cc
            integer :: gg, nsel, mm
            allocate(ji(0), jo(1), js(0))
            jo(1) = 1
            do gg = 1, kstar
                sel  = pack(idx, [(child(idx(mm)) == gg, mm=1,size(idx))])
                nsel = 0
                if( size(sel) >= MINSEL_JUNK )then
                    nsel = min(size(sel), neq)
                    call rnd_subset(sel, nsel)
                endif
                ji = [ji, sel(1:nsel)]
                jo = [jo, size(ji) + 1]
                js = [js, real(max(nsel,1)) / real(neq)]
            end do
            allocate(jm(self%nh,self%box,self%box,kstar))
            call self%maps(ji, jo, js, eps, .true., .false., jm)
            line = ' half '//side//' states:'
            do gg = 1, kstar
                nsel = count([(child(idx(mm)) == gg, mm=1,size(idx))])
                if( nsel < MINSEL_JUNK )then
                    line = line//' '//int2str(gg)//':(n='//int2str(nsel)//')'
                else
                    cc   = self%shell_corr(jm(:,:,:,gg), cons(:,:,:,hh), 2, self%klo - 1)
                    line = line//' '//int2str(gg)//':(n='//int2str(nsel)//', consensus corr '//fmtr(real(cc,dp),2)
                    if( cc < JUNK_CORR ) line = line//', junk-like'
                    line = line//')'
                endif
            end do
            write(logfhandle,'(A)') line
            deallocate(jm)
        end subroutine junk_report

    end subroutine test

    subroutine consensus( self, ia, ib, la, lb, k0, neq, eps, cons )
        class(state_mergeup), intent(inout) :: self
        integer,              intent(in)    :: ia(:), ib(:), la(:), lb(:), k0
        integer,              intent(out)   :: neq
        real,    allocatable, intent(out)   :: eps(:)
        complex, allocatable, intent(inout) :: cons(:,:,:,:)
        real, allocatable :: epsa(:), epsb(:)
        neq = min(minnz(la, k0), minnz(lb, k0))
        call self%accumulate(ia, 1, .false.)
        call self%recs(1)%get_rho_shell_medians(epsa)
        call self%accumulate(ib, 2, .false.)
        call self%recs(2)%get_rho_shell_medians(epsb)
        eps = EPS_REG * 0.5 * (epsa * real(neq) / real(size(ia)) + epsb * real(neq) / real(size(ib)))
        if( .not. allocated(cons) ) allocate(cons(self%nh,self%box,self%box,2))
        call self%finish(1, eps * real(size(ia)) / real(neq), .true., cons(:,:,:,1))
        call self%finish(2, eps * real(size(ib)) / real(neq), .true., cons(:,:,:,2))
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
        integer, allocatable :: units(:), jidx(:), joff(:), sa(:), sb(:), part(:)
        real,    allocatable :: jscale(:), met(:,:,:)
        logical, allocatable :: flag(:), keepa(:), keepb(:)
        integer :: nu, u, nfl, i
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
            call self%accumulate(sa, 1, .false.)
            call self%finish(1, eps * real(size(sa)) / real(neq), .true., sc(:,:,:,1))
            call self%accumulate(sb, 2, .false.)
            call self%finish(2, eps * real(size(sb)) / real(neq), .true., sc(:,:,:,2))
            call flags(sc, 0.5 * real(size(sa) + size(sb)))
        endif
        junk = pack(units, flag)
        part = pack(units, .not. flag .and. min(met(2,2,:), met(2,3,:)) < JUNK_PROJ)
        do u = 1, nu
            write(logfhandle,'(A)') '  junk screen unit '//int2str(units(u))//': projection on the screened consensus '//&
                &fmtr(real(met(2,2,u),dp),2)//' / '//fmtr(real(met(2,3,u),dp),2)//', rel. SSNR '//&
                &fmtr(real(met(3,2,u),dp),2)//' / '//fmtr(real(met(3,3,u),dp),2)//' -> '//merge('JUNK','ok  ',flag(u))
        end do
        if( size(part) > 0 ) write(logfhandle,'(A)') '  partial units (projection < '//fmtr(real(JUNK_PROJ,dp),2)//&
            &' with the node signal; tested as composition, not junk): '//joined(part)
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

    subroutine decide( self, ia, ib, la, lb, kmax, eps, neq, cons, kl, thr, tag, plist, gof, exa, exb, kout )
        class(state_mergeup), intent(inout) :: self
        integer,              intent(in)    :: ia(:), ib(:), la(:), lb(:), kmax, neq, kl
        real,                 intent(in)    :: eps(0:), thr
        complex,              intent(in)    :: cons(:,:,:,:)
        character(len=*),     intent(in)    :: tag
        integer, allocatable, intent(out)   :: plist(:), gof(:), exa(:), exb(:)
        integer,              intent(out)   :: kout
        complex, allocatable :: allm(:,:,:,:)
        integer, allocatable :: presa(:), presb(:), jidx(:), joff(:), mapa(:), mapb(:), groupof(:), pp(:), qq(:)
        real,    allocatable :: jscale(:), rmat(:,:), merged(:), lowc(:), vals(:)
        logical, allocatable :: l_alive(:)
        real    :: best
        integer :: neq_, g, npair, p, q, bp, bq, mm
        presa = pack([(g, g=1,kmax)], [(count(la == g) > 0, g=1,kmax)])
        presb = pack([(g, g=1,kmax)], [(count(lb == g) > 0, g=1,kmax)])
        neq_  = min(minval([(count(la == presa(g)), g=1,size(presa))]), minval([(count(lb == presb(g)), g=1,size(presb))]))
        allocate(jidx(0), joff(1), jscale(0))
        joff(1) = 1
        call add_jobs(ia, la, presa)
        call add_jobs(ib, lb, presb)
        allocate(allm(self%nh,self%box,self%box,size(jscale)))
        call self%maps(jidx, joff, jscale, eps, .true., .false., allm)
        !$omp parallel do default(shared) private(mm) schedule(dynamic) proc_bind(close)
        do mm = 1, size(jscale)
            call self%dock(allm(:,:,:,mm), cons(:,:,:,merge(2,1,mm <= size(presa))), self%kt, self%khi)
        end do
        !$omp end parallel do
        plist = pack(presa, [(any(presb == presa(g)), g=1,size(presa))])
        npair = size(plist)
        allocate(mapa(npair), mapb(npair))
        do mm = 1, npair
            mapa(mm) = findloc(presa, plist(mm), dim=1)
            mapb(mm) = size(presa) + findloc(presb, plist(mm), dim=1)
        end do
        allocate(exa(0:kmax), exb(0:kmax), source=0)
        call extras(presa, 0, 1, exa)
        call extras(presb, size(presa), 2, exb)
        allocate(lowc(npair))
        do mm = 1, npair
            lowc(mm) = self%shell_corr(allm(:,:,:,mapa(mm)) - cons(:,:,:,1), allm(:,:,:,mapb(mm)) - cons(:,:,:,2), &
                &2, self%klo - 2)
        end do
        write(logfhandle,'(A)') '['//tag//'] equal N '//int2str(neq_)//'; test band '//int2str(kl)//'..'//&
            &int2str(self%khi)//'; paired '//int2str(npair)//' groups (coarse-shell correlation median '//&
            &fmtr(real(medianr(lowc),dp),3)//')'
        deallocate(allm)
        groupof = [(p, p=1,npair)]
        allocate(rmat(npair,npair), source=0.)
        allocate(pp(0), qq(0))
        do p = 1, npair
            do q = p+1, npair
                pp = [pp, p]
                qq = [qq, q]
            end do
        end do
        if( size(pp) > 0 )then
            allocate(vals(size(pp)))
            call self%eqn_stats(ia, ib, la, lb, plist, groupof, pp, qq, eps, neq, cons, kl, vals)
            do mm = 1, size(pp)
                rmat(pp(mm),qq(mm)) = vals(mm)
            end do
            deallocate(vals)
        endif
        write(logfhandle,'(A)') 'pair statistics at the start (sorted): '//sorted_upper(rmat, npair)
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
            if( bp == 0 .or. best >= thr ) exit
            merged = [merged, best]
            l_alive(bq) = .false.
            where( groupof == bq ) groupof = bp
            pp = pack([(q, q=1,npair)], l_alive .and. [(q /= bp, q=1,npair)])
            if( size(pp) == 0 ) cycle
            qq = pp
            pp = [(min(bp,qq(mm)), mm=1,size(qq))]
            qq = [(max(bp,qq(mm)), mm=1,size(qq))]
            allocate(vals(size(pp)))
            call self%eqn_stats(ia, ib, la, lb, plist, groupof, pp, qq, eps, neq, cons, kl, vals)
            do mm = 1, size(pp)
                rmat(pp(mm),qq(mm)) = vals(mm)
            end do
            deallocate(vals)
        end do
        allocate(gof(npair))
        do mm = 1, npair
            gof(mm) = count(l_alive(1:groupof(mm)))
        end do
        kout = count(l_alive)
        write(logfhandle,'(A)') int2str(size(merged))//' merges at '//joinedr(merged)//' -> '//int2str(kout)//' groups'
        deallocate(rmat, l_alive, groupof, merged, lowc, mapa, mapb, jidx, joff, jscale)

    contains

        subroutine add_jobs( idx, lab, pres )
            integer, intent(in) :: idx(:), lab(:), pres(:)
            integer, allocatable :: sel(:)
            integer :: gg
            do gg = 1, size(pres)
                sel = pack(idx, lab == pres(gg))
                call rnd_subset(sel, neq_)
                jidx   = [jidx, sel(1:neq_)]
                joff   = [joff, size(jidx) + 1]
                jscale = [jscale, real(neq_) / real(neq)]
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

    end subroutine decide

    subroutine eqn_stats( self, ia, ib, la, lb, plist, groupof, pp, qq, eps, neq, cons, kl, vals )
        class(state_mergeup), intent(inout) :: self
        integer,              intent(in)    :: ia(:), ib(:), la(:), lb(:), plist(:), groupof(:), pp(:), qq(:), neq, kl
        real,                 intent(in)    :: eps(0:)
        complex,              intent(in)    :: cons(:,:,:,:)
        real,                 intent(out)   :: vals(:)
        complex, allocatable :: pm(:,:,:,:)
        integer, allocatable :: jidx(:), joff(:), ga(:), gb(:)
        real,    allocatable :: jscale(:)
        integer :: j, m
        ga = [(group_of_label(la(m)), m=1,size(la))]
        gb = [(group_of_label(lb(m)), m=1,size(lb))]
        allocate(jidx(0), joff(1), jscale(0))
        joff(1) = 1
        do j = 1, size(pp)
            call pair_jobs(ia, ga, pp(j), qq(j))
            call pair_jobs(ib, gb, pp(j), qq(j))
        end do
        allocate(pm(self%nh,self%box,self%box,size(jscale)))
        call self%maps(jidx, joff, jscale, eps, .true., .false., pm)
        !$omp parallel do default(shared) private(m) schedule(dynamic) proc_bind(close)
        do m = 1, size(jscale)
            call self%dock(pm(:,:,:,m), cons(:,:,:,merge(2,1,mod(m-1,4) < 2)), self%kt, self%khi)
        end do
        !$omp end parallel do
        !$omp parallel do default(shared) private(j) schedule(dynamic) proc_bind(close)
        do j = 1, size(pp)
            vals(j) = self%pair_stat(pm(:,:,:,4*j-3), pm(:,:,:,4*j-2), pm(:,:,:,4*j-1), pm(:,:,:,4*j), kl)
        end do
        !$omp end parallel do
        deallocate(pm, jidx, joff, jscale)

    contains

        integer function group_of_label( lab )
            integer, intent(in) :: lab
            integer :: mm
            group_of_label = 0
            do mm = 1, size(plist)
                if( plist(mm) == lab ) group_of_label = groupof(mm)
            end do
        end function group_of_label

        subroutine pair_jobs( idx, grp, p, q )
            integer, intent(in) :: idx(:), grp(:), p, q
            integer, allocatable :: sp(:), sq(:)
            integer :: nn
            allocate(sp, source=pack(idx, grp == p))
            allocate(sq, source=pack(idx, grp == q))
            nn = max(1, min(size(sp), size(sq)))
            call rnd_subset(sp, nn)
            call rnd_subset(sq, nn)
            jidx   = [jidx, sp(1:nn)]
            joff   = [joff, size(jidx) + 1]
            jidx   = [jidx, sq(1:nn)]
            joff   = [joff, size(jidx) + 1]
            jscale = [jscale, real(nn) / real(neq), real(nn) / real(neq)]
        end subroutine pair_jobs

    end subroutine eqn_stats

    subroutine cross_em( self, ia, ib, starta, startb, k, eps, neq, la, lb )
        class(state_mergeup), intent(inout) :: self
        integer,              intent(in)    :: ia(:), ib(:), starta(:), startb(:), k, neq
        real,                 intent(in)    :: eps(0:)
        integer,              intent(inout) :: la(:), lb(:)
        complex,             allocatable :: ma(:,:,:,:), mb(:,:,:,:)
        type(reconstructor), allocatable :: refa(:), refb(:)
        real,    allocatable :: wa(:,:), wb(:,:), wan(:,:), wbn(:,:), fl(:,:)
        integer, allocatable :: mpa(:), mpb(:), ca(:), cb(:)
        logical, allocatable :: l_alive(:)
        character(len=:), allocatable :: line
        integer :: it, g, i
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
        write(logfhandle,'(A)') 'cross-half EM: '//int2str(k)//' groups, '//int2str(NEM_ITERS)//' iterations, shells 2..'//&
            &int2str(self%kb)//', '//int2str(self%nem)//' plane samples per particle'
        do it = 1, NEM_ITERS
            call em_maps(ia, wa, ma)
            call em_maps(ib, wb, mb)
            !$omp parallel do default(shared) private(g) schedule(dynamic) proc_bind(close)
            do g = 1, 2*k
                if( g <= k )then
                    call self%dock(ma(:,:,:,g), mb(:,:,:,k+1), 2, self%kb)
                else
                    call self%dock(mb(:,:,:,g-k), ma(:,:,:,k+1), 2, self%kb)
                endif
            end do
            !$omp end parallel do
            !$omp parallel do default(shared) private(g) schedule(dynamic) proc_bind(close)
            do g = 1, k
                call filtered_refs(g)
            end do
            !$omp end parallel do
            call self%estep(ia, wa, refb, l_alive, wan)
            call self%estep(ib, wb, refa, l_alive, wbn)
            mpa  = maxloc(wa, dim=2)
            mpb  = maxloc(wb, dim=2)
            line = '  iteration '//int2str(it)//': MAP changed A '//&
                &fmtr(real(count(maxloc(wan,dim=2) /= mpa),dp)/real(size(ia),dp),3)//' B '//&
                &fmtr(real(count(maxloc(wbn,dim=2) /= mpb),dp)/real(size(ib),dp),3)
            wa  = wan
            wb  = wbn
            mpa = maxloc(wa, dim=2)
            mpb = maxloc(wb, dim=2)
            ca  = [(count(mpa == g), g=1,k)]
            cb  = [(count(mpb == g), g=1,k)]
            do g = 1, k
                if( .not. l_alive(g) ) cycle
                if( min(ca(g), cb(g)) < EM_MIN )then
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
        end do
        la = maxloc(wa, dim=2, mask=spread(l_alive, 1, size(ia)))
        lb = maxloc(wb, dim=2, mask=spread(l_alive, 1, size(ib)))
        do g = 1, k
            call refa(g)%kill
            call refb(g)%kill
        end do
        deallocate(ma, mb, refa, refb, wa, wb, fl, l_alive)

    contains

        subroutine em_maps( idx, w, m )
            integer, intent(in)    :: idx(:)
            real,    intent(in)    :: w(:,:)
            complex, intent(inout) :: m(:,:,:,:)
            integer, allocatable :: jidx(:), joff(:)
            real,    allocatable :: jscale(:), jw(:)
            integer :: gg, nn
            nn     = size(idx)
            jidx   = [[(idx, gg=1,k)], idx]
            joff   = [(1 + gg*nn, gg=0,k+1)]
            jw     = [reshape(w, [nn*k]), [(1., gg=1,nn)]]
            jscale = [[(max(sum(w(:,gg)), 1.)/real(neq), gg=1,k)], real(nn)/real(neq)]
            call self%maps(jidx, joff, jscale, eps, .false., .true., m, jw)
        end subroutine em_maps

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

    end subroutine cross_em

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
        real,    allocatable :: sm(:,:,:), wa(:,:), wb(:,:), rot(:,:,:), fw2(:), key(:)
        real(dp),allocatable :: epsb(:), epsf(:)
        integer, allocatable :: ia(:), ib(:), oa(:), ob(:), fhk(:,:)
        real(dp) :: bta(self%nshem), btb(self%nshem), s2a(self%nshem), s2b(self%nshem), pra(2), prb(2)
        real    :: f(0:self%nyq), cbest, chga, chgb
        integer :: it, g, ig, nh2, i, ngrid, nfull, h, k, r2
        call vol%read_and_crop(mfile, self%p%smpd, self%box, self%p%smpd_crop)
        sm = min(max(vol%get_rmat(), 0.), 1.) * self%msk3
        call vol%kill
        allocate(rot(3,3,self%n))
        do i = 1, self%n
            call self%os(1)%set_euler(self%eul(:,i))
            rot(:,:,i) = self%os(1)%get_mat()
        end do
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
        do it = 1, merge(0, SITE_ITS, size(score) > 0)
            call nn_map(ia, wa(:,1), self%emhk, self%emw, epsb, mp(:,:,:,1,1))
            call nn_map(ia, wa(:,2), self%emhk, self%emw, epsb, mp(:,:,:,2,1))
            call nn_map(ia, [(1., i=1,size(ia))], self%emhk, self%emw, epsb, mp(:,:,:,3,1))
            call nn_map(ib, wb(:,1), self%emhk, self%emw, epsb, mp(:,:,:,1,2))
            call nn_map(ib, wb(:,2), self%emhk, self%emw, epsb, mp(:,:,:,2,2))
            call nn_map(ib, [(1., i=1,size(ib))], self%emhk, self%emw, epsb, mp(:,:,:,3,2))
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
        do ig = 1, ngrid
            nh2 = SITE_GRID(ig) / 2
            call nn_map(ia(oa(1:min(nh2,size(ia)))), [(1., i=1,min(nh2,size(ia)))], fhk, fw2, epsf, tm(:,:,:,3))
            call nn_map(ib(ob(1:min(nh2,size(ib)))), [(1., i=1,min(nh2,size(ib)))], fhk, fw2, epsf, tm(:,:,:,4))
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

        subroutine cell_of( i, h, k, c1, c2, cnj )
            integer, intent(in)  :: i, h, k
            integer, intent(out) :: c1(3), c2(3)
            logical, intent(out) :: cnj
            integer :: p(3)
            p   = nint(real(h) * rot(1,:,i) + real(k) * rot(2,:,i))
            cnj = p(1) < 0
            if( cnj ) p = -p
            c1  = 0
            c2  = 0
            if( any(abs(p) >= self%nyq) ) return
            c1 = [p(1) + 1, modulo(p(2), self%box) + 1, modulo(p(3), self%box) + 1]
            if( p(1) == 0 ) c2 = [1, modulo(-p(2), self%box) + 1, modulo(-p(3), self%box) + 1]
        end subroutine cell_of

        function density( pidx, hk, sw ) result( epsv )
            integer, intent(in)   :: pidx(:), hk(:,:)
            real,    intent(in)   :: sw(:)
            real(dp), allocatable :: epsv(:)
            real(dp), allocatable :: den(:,:,:)
            real,     allocatable :: vals(:)
            integer :: ii, j, c1(3), c2(3), s, cnt
            logical :: cnj
            allocate(den(self%nh,self%box,self%box), source=0._dp)
            !$omp parallel do default(shared) private(ii,j,c1,c2,cnj) schedule(static) reduction(+:den) proc_bind(close)
            do ii = 1, size(pidx)
                do j = 1, size(sw)
                    call cell_of(pidx(ii), hk(1,j), hk(2,j), c1, c2, cnj)
                    if( c1(1) == 0 ) cycle
                    if( c2(1) > 0 )then
                        den(c1(1),c1(2),c1(3)) = den(c1(1),c1(2),c1(3)) + 0.5_dp * sw(j) * &
                            &self%fpls(pidx(ii))%ctfsq_plane(hk(1,j),hk(2,j))
                        den(c2(1),c2(2),c2(3)) = den(c2(1),c2(2),c2(3)) + 0.5_dp * sw(j) * &
                            &self%fpls(pidx(ii))%ctfsq_plane(hk(1,j),hk(2,j))
                    else
                        den(c1(1),c1(2),c1(3)) = den(c1(1),c1(2),c1(3)) + sw(j) * self%fpls(pidx(ii))%ctfsq_plane(hk(1,j),hk(2,j))
                    endif
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
            complex(dp) :: v
            real(dp)    :: c2v, wsum
            integer     :: ii, j, c1(3), c2(3), hh, kk, ll
            logical     :: cnj
            allocate(num(self%nh,self%box,self%box), source=cmplx(0._dp,0._dp,dp))
            allocate(den(self%nh,self%box,self%box), source=0._dp)
            !$omp parallel do default(shared) private(ii,j,c1,c2,cnj,v,c2v) schedule(static) reduction(+:num,den) proc_bind(close)
            do ii = 1, size(pidx)
                if( pw(ii) <= 0. ) cycle
                do j = 1, size(sw)
                    call cell_of(pidx(ii), hk(1,j), hk(2,j), c1, c2, cnj)
                    if( c1(1) == 0 ) cycle
                    v   = real(pw(ii) * sw(j), dp) * cmplx(self%fpls(pidx(ii))%cmplx_plane(hk(1,j),hk(2,j)), kind=dp)
                    c2v = real(pw(ii) * sw(j), dp) * real(self%fpls(pidx(ii))%ctfsq_plane(hk(1,j),hk(2,j)), dp)
                    if( cnj ) v = conjg(v)
                    if( c2(1) > 0 )then
                        num(c1(1),c1(2),c1(3)) = num(c1(1),c1(2),c1(3)) + 0.5_dp * v
                        den(c1(1),c1(2),c1(3)) = den(c1(1),c1(2),c1(3)) + 0.5_dp * c2v
                        num(c2(1),c2(2),c2(3)) = num(c2(1),c2(2),c2(3)) + 0.5_dp * conjg(v)
                        den(c2(1),c2(2),c2(3)) = den(c2(1),c2(2),c2(3)) + 0.5_dp * c2v
                    else
                        num(c1(1),c1(2),c1(3)) = num(c1(1),c1(2),c1(3)) + v
                        den(c1(1),c1(2),c1(3)) = den(c1(1),c1(2),c1(3)) + c2v
                    endif
                end do
            end do
            !$omp end parallel do
            wsum = max(sum(real(pw,dp)), 1._dp)
            do ll = 1, self%box
                do kk = 1, self%box
                    do hh = 1, self%nh
                        out(hh,kk,ll) = cmplx(num(hh,kk,ll) / (den(hh,kk,ll) + epsv(self%fsh(hh,kk,ll)) * wsum), kind=kind(out))
                    end do
                end do
            end do
        end subroutine nn_map

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
            integer  :: ii, gg, j, c1(3), c2(3), sh
            logical  :: cnj
            allocate(aa(self%nshem,2,size(pidx)), bb(self%nshem,2,size(pidx)), source=0.)
            !$omp parallel do default(shared) private(ii,gg,j,c1,c2,cnj,sh,pv,dv) schedule(static) proc_bind(close)
            do ii = 1, size(pidx)
                do j = 1, self%nem
                    call cell_of(pidx(ii), self%emhk(1,j), self%emhk(2,j), c1, c2, cnj)
                    if( c1(1) == 0 ) cycle
                    sh = self%emsh(j)
                    dv = self%fpls(pidx(ii))%cmplx_plane(self%emhk(1,j),self%emhk(2,j))
                    do gg = 1, 2
                        pv = refs(c1(1),c1(2),c1(3),gg)
                        if( cnj ) pv = conjg(pv)
                        aa(sh,gg,ii) = aa(sh,gg,ii) + self%emw(j) * (real(pv)*real(dv) + aimag(pv)*aimag(dv))
                        bb(sh,gg,ii) = bb(sh,gg,ii) + self%emw(j) * &
                            &self%fpls(pidx(ii))%ctfsq_plane(self%emhk(1,j),self%emhk(2,j)) * abs(pv)**2
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

    ! per-node junk screen at the seeding labels: amplitudes against every state's other-half reference, 3-component mixture cut
    subroutine junk_screen( self, seeds )
        use simple_srch_sort_loc, only: hpsort
        class(state_mergeup), intent(inout) :: self
        integer,              intent(inout) :: seeds(:)
        type(ivec)            :: hx(2)
        complex,  allocatable :: mp(:,:,:,:,:), rf(:,:,:,:)
        real,     allocatable :: rot(:,:,:), fs(:,:), amp(:), var(:), post(:), x(:), v(:), key(:)
        real(dp), allocatable :: epsv(:)
        integer,  allocatable :: lab(:), best(:), cnt(:,:), ord(:)
        logical,  allocatable :: ljunk(:), lj(:)
        character(len=:), allocatable :: line
        real(dp) :: mu(3), s2(3), pw(3), gap(3), xm
        real     :: f(0:self%nyq), mv
        integer  :: nk, i, g, e, c, cg, h, k, l, funit, rep, ios, nf
        nk = maxval([0, seeds])
        if( nk == 0 )then
            write(logfhandle,'(A)') '>>> JUNK_SCREEN no seed state, nothing screened'
            return
        endif
        lab = seeds
        call fopen(rep, string('site_junk_report.txt'), status='replace', action='write', iostat=ios)
        call fileiochk('junk_screen: cannot write site_junk_report.txt', ios)
        allocate(rot(3,3,self%n))
        do i = 1, self%n
            call self%os(1)%set_euler(self%eul(:,i))
            rot(:,:,i) = self%os(1)%get_mat()
        end do
        allocate(cnt(nk,2))
        do e = 1, 2
            hx(e)%v  = pack([(i, i=1,self%n)], self%eo == e - 1)
            cnt(:,e) = [(count(lab(hx(e)%v) == g), g=1,nk)]
        end do
        call say(int2str(nk)//' seed states, '//int2str(self%n)//' particles, band shells 2..'//int2str(self%kj)//' ('//&
            &fmtr(real(real(self%box)*self%p%smpd_crop/real(self%kj),dp),1)//' A, '//int2str(self%nj)//' plane samples), '//&
            &'every particle against the other half''s state references')
        allocate(epsv(0:self%nyq))
        epsv = 0.5_dp * (dens(hx(1)%v) + dens(hx(2)%v))
        allocate(mp(self%nh,self%box,self%box,0:nk,2))
        do e = 1, 2
            call wiener(hx(e)%v, mp(:,:,:,0,e))
            do g = 1, nk
                call wiener(pack(hx(e)%v, lab(hx(e)%v) == g), mp(:,:,:,g,e))
            end do
        end do
        allocate(fs(0:self%nyq,nk))
        do g = 1, nk
            call self%shell_fsc(mp(:,:,:,g,1) - mp(:,:,:,0,1), mp(:,:,:,g,2) - mp(:,:,:,0,2), f)
            fs(:,g) = min(max(f, 0.), 1.)
        end do
        allocate(rf(self%nh,self%box,self%box,nk))
        allocate(amp(self%n), var(self%n), post(self%n), source=0.)
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
        deallocate(mp, rf, rot)
        x = pack(amp, best > 0)
        v = pack(var, best > 0)
        allocate(ljunk(self%n), source=.false.)
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
            if( any(lj) )then
                do i = 1, self%n
                    if( best(i) > 0 ) post(i) = jpost(amp(i), var(i))
                end do
            endif
            ljunk = best > 0 .and. post >= 0.5
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
            call say('junk components (weight >= '//fmtr(real(JUNK_MINW,dp),3)//', >= '//int2str(nint(JUNK_GAP))//&
                &' sd below the genuine one):'//line)
        endif
        nf = count(ljunk)
        call say(int2str(nf)//' of '//int2str(self%n)//' particles flagged junk ('//fmtr(real(nf,dp)/real(self%n,dp),4)//&
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
        write(funit,'(A)') '# project row, junk (1) or not (0), amplitude, noise variance, best state, junk posterior'
        do i = 1, self%n
            write(funit,'(I0,1X,I0,1X,F9.4,1X,ES11.4,1X,I0,1X,F7.4)') self%rows(i), merge(1, 0, ljunk(i)), amp(i), var(i), &
                &best(i), post(i)
        end do
        call fclose(funit)
        where( ljunk ) seeds = 0

    contains

        subroutine say( str )
            character(len=*), intent(in) :: str
            write(logfhandle,'(A)') '>>> JUNK_SCREEN '//str
            write(rep,'(A)') str
        end subroutine say

        subroutine cell( i, h, k, c1, c2, cnj )
            integer, intent(in)  :: i, h, k
            integer, intent(out) :: c1(3), c2(3)
            logical, intent(out) :: cnj
            integer :: p(3)
            p   = nint(real(h) * rot(1,:,i) + real(k) * rot(2,:,i))
            cnj = p(1) < 0
            if( cnj ) p = -p
            c1  = 0
            c2  = 0
            if( any(abs(p) >= self%nyq) ) return
            c1 = [p(1) + 1, modulo(p(2), self%box) + 1, modulo(p(3), self%box) + 1]
            if( p(1) == 0 ) c2 = [1, modulo(-p(2), self%box) + 1, modulo(-p(3), self%box) + 1]
        end subroutine cell

        function dens( pidx ) result( ev )
            integer, intent(in)   :: pidx(:)
            real(dp), allocatable :: ev(:)
            real(dp), allocatable :: den(:,:,:)
            real,     allocatable :: vals(:)
            real(dp) :: c2v
            integer  :: ii, j, c1(3), c2(3), s, nv
            logical  :: cnj
            allocate(den(self%nh,self%box,self%box), source=0._dp)
            !$omp parallel do default(shared) private(ii,j,c1,c2,cnj,c2v) schedule(static) reduction(+:den) proc_bind(close)
            do ii = 1, size(pidx)
                do j = 1, self%nj
                    call cell(pidx(ii), self%jhk(1,j), self%jhk(2,j), c1, c2, cnj)
                    if( c1(1) == 0 ) cycle
                    c2v = real(self%jw(j) * self%fpls(pidx(ii))%ctfsq_plane(self%jhk(1,j),self%jhk(2,j)), dp)
                    if( c2(1) > 0 )then
                        den(c1(1),c1(2),c1(3)) = den(c1(1),c1(2),c1(3)) + 0.5_dp * c2v
                        den(c2(1),c2(2),c2(3)) = den(c2(1),c2(2),c2(3)) + 0.5_dp * c2v
                    else
                        den(c1(1),c1(2),c1(3)) = den(c1(1),c1(2),c1(3)) + c2v
                    endif
                end do
            end do
            !$omp end parallel do
            allocate(ev(0:self%nyq), source=1._dp)
            do s = 0, self%nyq
                vals = real(pack(den, self%fsh == s .and. den > 0._dp))
                nv   = size(vals)
                if( nv == 0 ) cycle
                call hpsort(vals)
                ev(s) = real(vals((nv + 1) / 2), dp) / real(max(1, size(pidx)), dp)
            end do
        end function dens

        subroutine wiener( pidx, out )
            integer, intent(in)  :: pidx(:)
            complex, intent(out) :: out(:,:,:)
            complex(dp), allocatable :: num(:,:,:)
            real(dp),    allocatable :: den(:,:,:)
            complex(dp) :: cv
            real(dp)    :: c2v
            integer     :: ii, j, c1(3), c2(3), hh, kk, ll
            logical     :: cnj
            allocate(num(self%nh,self%box,self%box), source=cmplx(0._dp,0._dp,dp))
            allocate(den(self%nh,self%box,self%box), source=0._dp)
            !$omp parallel do default(shared) private(ii,j,c1,c2,cnj,cv,c2v) schedule(static) reduction(+:num,den) proc_bind(close)
            do ii = 1, size(pidx)
                do j = 1, self%nj
                    call cell(pidx(ii), self%jhk(1,j), self%jhk(2,j), c1, c2, cnj)
                    if( c1(1) == 0 ) cycle
                    cv  = real(self%jw(j),dp) * cmplx(self%fpls(pidx(ii))%cmplx_plane(self%jhk(1,j),self%jhk(2,j)), kind=dp)
                    c2v = real(self%jw(j) * self%fpls(pidx(ii))%ctfsq_plane(self%jhk(1,j),self%jhk(2,j)), dp)
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
            !$omp end parallel do
            do ll = 1, self%box
                do kk = 1, self%box
                    do hh = 1, self%nh
                        out(hh,kk,ll) = cmplx(num(hh,kk,ll) / (den(hh,kk,ll) + epsv(self%fsh(hh,kk,ll)) * &
                            &real(max(1, size(pidx)),dp)), kind=kind(out))
                    end do
                end do
            end do
        end subroutine wiener

        subroutine score( pidx, loff, e )
            integer, intent(in) :: pidx(:), e
            logical, intent(in) :: loff(:)
            real(dp) :: sa(self%kj), sb(self%kj), sf(self%kj), sc(self%kj), beta(self%kj), wsh(self%kj)
            real(dp) :: a(self%kj,nk), b(self%kj,nk), num, den, s, smax
            real     :: cw(self%kj), c2v
            complex  :: pv, dv
            integer  :: ii, i, j, g, c1(3), c2(3), sh, nd
            logical  :: cnj
            cw = 0.
            do j = 1, self%nj
                cw(self%jsh(j)) = cw(self%jsh(j)) + self%jw(j)
            end do
            sa = 0._dp
            sb = 0._dp
            sf = 0._dp
            sc = 0._dp
            !$omp parallel do default(shared) private(ii,i,j,g,c1,c2,cnj,sh,pv,dv,c2v) schedule(static) &
            !$omp& reduction(+:sa,sb,sf,sc) proc_bind(close)
            do ii = 1, size(pidx)
                i = pidx(ii)
                g = lab(i)
                if( g < 1 ) cycle
                if( .not. loff(g) ) cycle
                do j = 1, self%nj
                    call cell(i, self%jhk(1,j), self%jhk(2,j), c1, c2, cnj)
                    if( c1(1) == 0 ) cycle
                    sh  = self%jsh(j)
                    dv  = self%fpls(i)%cmplx_plane(self%jhk(1,j),self%jhk(2,j))
                    c2v = self%fpls(i)%ctfsq_plane(self%jhk(1,j),self%jhk(2,j))
                    pv  = rf(c1(1),c1(2),c1(3),g)
                    if( cnj ) pv = conjg(pv)
                    sa(sh) = sa(sh) + self%jw(j) * (real(pv)*real(dv) + aimag(pv)*aimag(dv))
                    sb(sh) = sb(sh) + self%jw(j) * c2v * abs(pv)**2
                end do
                sf = sf + real(self%f2j(:,i),dp)
                sc = sc + real(cw,dp)
            end do
            !$omp end parallel do
            beta = min(max(sa / max(sb, 1.e-30_dp), 0._dp), 2._dp)
            wsh  = 1._dp / max((sf - 2._dp*beta*sa + beta*beta*sb) / max(sc, 1._dp), 1.e-12_dp)
            !$omp parallel do default(shared) private(ii,i,j,g,c1,c2,cnj,sh,pv,dv,c2v,a,b,num,den,s,smax) schedule(static) &
            !$omp& proc_bind(close)
            do ii = 1, size(pidx)
                i = pidx(ii)
                a = 0._dp
                b = 0._dp
                do j = 1, self%nj
                    call cell(i, self%jhk(1,j), self%jhk(2,j), c1, c2, cnj)
                    if( c1(1) == 0 ) cycle
                    sh  = self%jsh(j)
                    dv  = self%fpls(i)%cmplx_plane(self%jhk(1,j),self%jhk(2,j))
                    c2v = self%fpls(i)%ctfsq_plane(self%jhk(1,j),self%jhk(2,j))
                    do g = 1, nk
                        pv = rf(c1(1),c1(2),c1(3),g)
                        if( cnj ) pv = conjg(pv)
                        a(sh,g) = a(sh,g) + self%jw(j) * (real(pv)*real(dv) + aimag(pv)*aimag(dv))
                        b(sh,g) = b(sh,g) + self%jw(j) * c2v * abs(pv)**2
                    end do
                end do
                smax = -huge(smax)
                do g = 1, nk
                    if( .not. loff(g) ) cycle
                    num = sum(wsh * beta * a(:,g))
                    den = sum(wsh * beta * beta * b(:,g))
                    s   = num - 0.5_dp * den
                    if( s <= smax ) cycle
                    smax    = s
                    best(i) = g
                    amp(i)  = real(num / max(den, 1.e-12_dp))
                    var(i)  = real(1._dp / max(den, 1.e-12_dp))
                end do
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

    real function halves_merge_margin( cm, ra, rb, a, b )
        real,    intent(in) :: cm(:,:), ra(:), rb(:)
        integer, intent(in) :: a, b
        halves_merge_margin = min(nm(cm(a,b), ra(a), rb(b)), nm(cm(b,a), ra(b), rb(a))) - &
            &min(nm(cm(a,a), ra(a), rb(a)), nm(cm(b,b), ra(b), rb(b)))
    contains
        real function nm( c, r1, r2 )
            real, intent(in) :: c, r1, r2
            nm = c / sqrt(max(r1, 1.e-3) * max(r2, 1.e-3))
        end function nm
    end function halves_merge_margin

    real function halves_accept_margin( cm, gp, cand, a, l_row )
        real,    intent(in) :: cm(:,:)
        integer, intent(in) :: gp(:), a
        logical, intent(in) :: cand(:), l_row
        real    :: best
        integer :: f
        best = -huge(best)
        do f = 1, size(gp)
            if( gp(f) == gp(a) .or. .not. cand(f) ) cycle
            if( l_row )then
                best = max(best, cm(a,f))
            else
                best = max(best, cm(f,a))
            endif
        end do
        halves_accept_margin = cm(a,a) - best
        if( best == -huge(best) ) halves_accept_margin = cm(a,a)
    end function halves_accept_margin

    real function halves_se( reps )
        real, intent(in) :: reps(:)
        halves_se = sqrt(sum((reps - sum(reps) / real(size(reps)))**2) / real(max(1, size(reps) - 1))) / sqrt(2.)
    end function halves_se

    real function halves_floor( cq, r1q, r2q, a, b )
        real,    intent(in) :: cq(:,:,:), r1q(:,:), r2q(:,:)
        integer, intent(in) :: a, b
        integer :: k
        halves_floor = max(halves_se([(cq(a,a,k) / sqrt(max(r1q(a,k), 1.e-3) * max(r2q(a,k), 1.e-3)), k=1,size(cq,3))]), &
            &halves_se([(cq(b,b,k) / sqrt(max(r1q(b,k), 1.e-3) * max(r2q(b,k), 1.e-3)), k=1,size(cq,3))]))
    end function halves_floor

    function halves_cand( craw, rel1, rel2, cq ) result( cand )
        real,    intent(in)  :: craw(:,:), rel1(:), rel2(:), cq(:,:,:)
        logical, allocatable :: cand(:)
        integer :: a
        cand = [(craw(a,a) > 0. .and. rel1(a) > 0. .and. rel2(a) > 0. .and. craw(a,a) >= HALVES_KSIG * halves_se(cq(a,a,:)), &
            &a=1,size(craw,1))]
    end function halves_cand

    ! groups merge, closest pair first, when every cross pair of candidates is within HALVES_KSIG SEs (floored at self-match SEs); a
    ! group counts when each member beats its best outside partner by HALVES_KSIG SEs; solo ones never merge, tested last vs poles
    subroutine merge_poles_halves( posts, craw, rel1, rel2, cq, r1q, r2q, minseed, l_rest, solo, grp, seeds, ok, twin, gm )
        use simple_srch_sort_loc, only: hpsort
        real,                 intent(in)  :: posts(:,:), craw(:,:), rel1(:), rel2(:), cq(:,:,:), r1q(:,:), r2q(:,:)
        integer,              intent(in)  :: minseed
        logical,              intent(in)  :: l_rest, solo(:)
        integer, allocatable, intent(out) :: grp(:), seeds(:), twin(:), gm(:)
        logical, allocatable, intent(out) :: ok(:)
        real(dp), allocatable :: gpost(:,:)
        real,     allocatable :: rat(:,:), key(:)
        integer,  allocatable :: lab(:), plab(:), pa(:), pb(:), ord(:)
        logical,  allocatable :: cand(:), cons(:), other(:), part(:)
        integer :: ne, n, a, b, ga, gb, npole, ipass, br, bc, k
        ne = size(craw,1)
        n  = size(posts,1)
        cand = halves_cand(craw, rel1, rel2, cq)
        part = cand .and. .not. solo
        allocate(rat(ne,ne), source=-huge(1.))
        allocate(pa(0), pb(0), key(0))
        do a = 1, ne
            do b = a + 1, ne
                if( .not. (cand(a) .and. cand(b)) .or. solo(a) .or. solo(b) ) cycle
                rat(a,b) = halves_merge_margin(craw, rel1, rel2, a, b) / max(halves_se([(halves_merge_margin(cq(:,:,k), &
                    &r1q(:,k), r2q(:,k), a, b), k=1,size(cq,3))]), halves_floor(cq, r1q, r2q, a, b), 1.e-6)
                rat(b,a) = rat(a,b)
                pa  = [pa, a]
                pb  = [pb, b]
                key = [key, -rat(a,b)]
            end do
        end do
        ord = [(k, k=1,size(key))]
        if( size(key) > 1 ) call hpsort(key, ord)
        grp = [(a, a=1,ne)]
        do k = 1, size(ord)
            a = pa(ord(k))
            b = pb(ord(k))
            if( rat(a,b) < -HALVES_KSIG .or. grp(a) == grp(b) ) cycle
            if( any(spread(grp == grp(a), 2, ne) .and. spread(grp == grp(b), 1, ne) .and. rat < -HALVES_KSIG) ) cycle
            ga = grp(a)
            gb = grp(b)
            where( grp == gb ) grp = ga
        end do
        allocate(cons(ne), ok(ne))
        allocate(twin(ne), source=0)
        do ipass = 1, 2
            do a = 1, ne
                cons(a) = cand(a) .and. strong(a, .true.) .and. strong(a, .false.)
            end do
            if( ipass == 2 ) exit
            do a = 1, ne
                if( cons(a) .or. .not. cand(a) .or. solo(a) ) cycle
                other = [(b /= a .and. cand(b) .and. .not. solo(b), b=1,ne)]
                if( .not. any(other) ) cycle
                br = maxloc(craw(a,:), dim=1, mask=other)
                bc = maxloc(craw(:,a), dim=1, mask=other)
                if( br /= bc ) cycle
                if( any(spread(grp == grp(a), 2, ne) .and. spread(grp == grp(br), 1, ne) .and. rat < -HALVES_KSIG) ) cycle
                twin(a) = br
                ga = grp(a)
                gb = grp(br)
                where( grp == ga ) grp = gb
            end do
        end do
        gm = grp
        do a = 1, ne
            ok(a) = all(pack(cons, grp == grp(a)))
        end do
        part = ok .and. .not. solo
        do a = 1, ne
            if( solo(a) ) ok(a) = cand(a) .and. strong(a, .true.) .and. strong(a, .false.)
        end do
        allocate(lab(ne), source=0)
        npole = 0
        do a = 1, ne
            if( .not. ok(a) .or. lab(a) > 0 ) cycle
            npole = npole + 1
            where( grp == grp(a) ) lab = npole
        end do
        grp = lab
        allocate(seeds(n), source=0)
        if( npole == 0 ) return
        allocate(gpost(n,npole), source=0._dp)
        do a = 1, ne
            if( grp(a) == 0 ) cycle
            gpost(:,grp(a)) = gpost(:,grp(a)) + real(posts(:,a),dp) / real(count(grp == grp(a)),dp)
        end do
        call seed_poles(gpost, minseed, l_rest, seeds, plab)
        do a = 1, ne
            if( grp(a) > 0 ) grp(a) = plab(grp(a))
        end do

    contains

        logical function strong( e, l_row )
            integer, intent(in) :: e
            logical, intent(in) :: l_row
            real :: am
            integer :: kk
            am     = halves_accept_margin(craw, grp, part, e, l_row)
            strong = am > 0. .and. am >= HALVES_KSIG * &
                &halves_se([(halves_accept_margin(cq(:,:,kk), grp, part, e, l_row), kk=1,size(cq,3))])
        end function strong

    end subroutine merge_poles_halves

    subroutine estep( self, idx, w, refs, l_alive, wnew )
        class(state_mergeup), intent(inout) :: self
        integer,              intent(in)    :: idx(:)
        real,                 intent(in)    :: w(:,:)
        type(reconstructor),  intent(in)    :: refs(:)
        logical,              intent(in)    :: l_alive(:)
        real, allocatable,    intent(out)   :: wnew(:,:)
        real,    allocatable :: aa(:,:,:), bb(:,:,:)
        integer, allocatable :: mp(:)
        real(dp) :: sa(self%nshem), sb(self%nshem), sf(self%nshem), sc(self%nshem), scl(self%nshem), sig2(self%nshem)
        real(dp) :: prior(size(refs)), s(size(refs)), smax
        complex  :: pv, wv
        integer  :: n, k, i, ip, g, j, h, kk, sh, it, nyq_save
        n  = size(idx)
        k  = size(refs)
        mp = maxloc(w, dim=2)
        do g = 1, k
            prior(g) = log(max(real(sum(w(:,g)),dp) / real(n,dp), 1.e-6_dp))
        end do
        allocate(aa(self%nshem,k,n), bb(self%nshem,k,n), source=0.)
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
            self%fpls(ip)%nyq = OSMPL_PAD_FAC * self%kb
            do g = 1, k
                if( .not. l_alive(g) .and. g /= mp(i) ) cycle
                call refs(g)%project_fplane(self%os(it), self%fpls(ip), self%pouts(it), apply_ctf_amp=.false., &
                    &samples=self%smps(it))
                do j = 1, self%nem
                    h  = self%emhk(1,j)
                    kk = self%emhk(2,j)
                    sh = self%emsh(j)
                    pv = self%pouts(it)%cmplx_plane(h,kk)
                    wv = self%fpls(ip)%cmplx_plane(h,kk)
                    aa(sh,g,i) = aa(sh,g,i) + self%emw(j) * (real(pv)*real(wv) + aimag(pv)*aimag(wv))
                    bb(sh,g,i) = bb(sh,g,i) + self%emw(j) * self%fpls(ip)%ctfsq_plane(h,kk) * abs(pv)**2
                end do
            end do
            self%fpls(ip)%nyq = nyq_save
            sa = sa + real(aa(:,mp(i),i),dp)
            sb = sb + real(bb(:,mp(i),i),dp)
            sf = sf + real(self%f2s(:,ip),dp)
            sc = sc + real(self%cw,dp)
        end do
        !$omp end parallel do
        scl  = sa / max(sb, 1.e-30_dp)
        sig2 = max((sf - 2._dp*scl*sa + scl*scl*sb) / max(sc, 1._dp), 1.e-12_dp)
        allocate(wnew(n,k))
        !$omp parallel do default(shared) private(i,g,s,smax) schedule(static) proc_bind(close)
        do i = 1, n
            do g = 1, k
                if( l_alive(g) )then
                    s(g) = 0.5_dp * sum((2._dp*scl*real(aa(:,g,i),dp) - scl*scl*real(bb(:,g,i),dp)) / sig2) + prior(g)
                else
                    s(g) = -huge(1._dp)
                endif
            end do
            smax = maxval(s)
            do g = 1, k
                wnew(i,g) = 0.
                if( l_alive(g) ) wnew(i,g) = real(exp(s(g) - smax))
            end do
            if( sum(wnew(i,:)) > 0. )then
                wnew(i,:) = wnew(i,:) / sum(wnew(i,:))
            else
                wnew(i,:) = w(i,:)
            endif
        end do
        !$omp end parallel do
        write(logfhandle,'(A)') '  cross-half E-step: model scale at shells 2 / '//int2str(self%kb)//' '//&
            &fmtr(scl(1),4)//' / '//fmtr(scl(self%nshem),4)
        deallocate(aa, bb, mp)
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

    subroutine finish( self, it, epsv, l_mask, cm )
        class(state_mergeup), intent(inout) :: self
        integer,              intent(in)    :: it
        real,                 intent(in)    :: epsv(0:)
        logical,              intent(in)    :: l_mask
        complex,              intent(out)   :: cm(:,:,:)
        real, pointer :: rmat(:,:,:)
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

    subroutine maps( self, jidx, joff, jscale, eps, l_mask, l_band, cms, jw )
        class(state_mergeup), intent(inout) :: self
        integer,              intent(in)    :: jidx(:), joff(:)
        real,                 intent(in)    :: jscale(:), eps(0:)
        logical,              intent(in)    :: l_mask, l_band
        complex,              intent(inout) :: cms(:,:,:,:)
        real, optional,       intent(in)    :: jw(:)
        integer :: j, it
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

    real function pair_stat( self, a1, a2, b1, b2, kl )
        class(state_mergeup), intent(in) :: self
        complex,              intent(in) :: a1(:,:,:), a2(:,:,:), b1(:,:,:), b2(:,:,:)
        integer,              intent(in) :: kl
        real    :: f(0:self%nyq)
        integer :: k
        call self%shell_fsc(a1 - a2, b1 - b2, f)
        pair_stat = -huge(1.)
        do k = kl, max(self%khi - 2, kl)
            pair_stat = max(pair_stat, sum(f(k:k+2)) / 3.)
        end do
    end function pair_stat

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
        if( allocated(self%msk3) ) deallocate(self%msk3, self%fsh, self%fw)
        if( allocated(self%fqh) ) deallocate(self%fqh, self%fqk)
        if( allocated(self%lmsk) ) deallocate(self%lmsk, self%lslot, self%msk0)
        self%p => null()
        self%b => null()
        self%n = 0
        self%l_exists = .false.
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

    subroutine rnd_subset( sel, m )
        use simple_rnd, only: partial_shuffle
        integer, intent(inout) :: sel(:)
        integer, intent(in)    :: m
        call partial_shuffle(sel, m)
    end subroutine rnd_subset

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

end module simple_state_mergeup
