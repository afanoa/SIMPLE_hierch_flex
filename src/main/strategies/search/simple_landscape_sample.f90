!@descr: landscape sampling along the pole pairs of a node: pair log-likelihood-ratio coordinates per Fourier pixel-parity set, pixel-split cross-validated sample counts, held-out and final sample maps
module simple_landscape_sample
!$ use omp_lib
use simple_defs,         only: dp, logfhandle, MAXIMGBATCHSZ, PI
use simple_type_defs,    only: fplane_type
use simple_error,        only: simple_exception
use simple_string,       only: string
use simple_string_utils, only: int2str
use simple_fileio,       only: fopen, fclose, fileiochk, nlines
use simple_ori,          only: ori
use simple_sym,          only: sym
use simple_image,        only: image
use simple_parameters,   only: parameters
use simple_builder,      only: builder
implicit none

public :: landscape_sample
private
#include "simple_local_flags.inc"

integer, parameter :: NVIEW = 40, NDOCK_ITS = 8, MINHALF = 150, NQB = 40, NCAND = 11
integer, parameter :: CANDS(NCAND) = [4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 24]
real,    parameter :: EM_LIM_LS = 12., QLO = 0.01, ALPHA = 0.01, MINFSC = 0.143, FSC_RES = 20., EPS_REG = 0.05, DOCK_MAX = 3.

type :: fold_res
    logical :: ok = .false.
    integer :: npr = 0
    real    :: zc = 0.
    integer, allocatable :: cell(:), grid(:), prs(:,:)
    real,    allocatable :: u(:), st(:,:), fs(:)
end type fold_res

type :: landscape_sample
    private
    class(parameters), pointer :: p => null()
    class(builder),    pointer :: b => null()
    type(fplane_type), allocatable :: fpls(:)
    type(image) :: wimg
    integer, allocatable :: rows(:), eo(:), st(:), vbin(:), hk(:,:), psh(:), ppar(:), fsh(:,:,:)
    real,    allocatable :: rot(:,:,:,:), dfm(:), pw(:), msk3(:,:,:), fw(:,:,:), sc(:,:,:), ppn(:), f2(:,:,:), cw(:,:), nw(:,:,:)
    real(dp),allocatable :: fqh(:), fqk(:), epsv(:)
    complex, allocatable :: refs(:,:,:,:,:)
    real    :: pfrac(-1:1) = 1.
    integer :: n = 0, box = 0, nh = 0, nyq = 0, kb = 0, kres = 0, np = 0, nst = 0, nrows = 0, nsym = 1
    logical :: l_exists = .false.
  contains
    procedure :: new
    procedure :: run
    procedure :: kill
    procedure, private :: read_refs
    procedure, private :: score
    procedure, private :: sample_pair
    procedure, private :: cell_of
    procedure, private :: nnmap
    procedure, private :: ustat
    procedure, private :: struct_centres
    procedure, private :: evaluate
    procedure, private :: merge_sparse
    procedure, private :: cellfsc
    procedure, private :: dock
    procedure, private :: shell_fsc
    procedure, private :: write_map
end type landscape_sample

contains

    subroutine new( self, params, build )
        use simple_matcher_3Drec,   only: init_rec, prep_imgs4rec, cleanup_rec_buffers
        use simple_matcher_ptcl_io, only: prepimgbatch, ptcl_read_session
        use simple_imgarr_utils,    only: alloc_imgarr, dealloc_imgarr
        use simple_defs,            only: nthr_glob
        use simple_srch_sort_loc,   only: hpsort
        class(landscape_sample),   intent(inout) :: self
        class(parameters), target, intent(inout) :: params
        class(builder),    target, intent(inout) :: build
        type(fplane_type), allocatable :: fdummy(:)
        type(image),       allocatable :: crop_imgs(:)
        real,              allocatable :: ones(:,:,:), sdn(:), vals(:)
        real(dp),          allocatable :: den(:,:,:)
        type(ptcl_read_session) :: reads
        type(image) :: mskimg
        type(ori)   :: o
        type(sym)   :: se
        real    :: dirs(3,NVIEW), v(3), t, z, rz, g, e1, e2, pws(-1:1), es(3)
        integer :: i, j, ip, h, k, l, r2, b1, b2, ie, sh, q, c1(3), c2(3), cnt, is
        logical :: cnj
        call self%kill
        self%p => params
        self%b => build
        !$ call omp_set_max_active_levels(1)
        self%nrows = build%spproj_field%get_noris()
        self%rows  = pack([(i, i=1,self%nrows)], [(build%spproj_field%get_state(i) > 0 .and. &
            &(build%spproj_field%get_eo(i) == 0 .or. build%spproj_field%get_eo(i) == 1), i=1,self%nrows)])
        self%n = size(self%rows)
        if( self%n < 2 ) THROW_HARD('no particles with a state and an even/odd half')
        self%box  = params%box_crop
        self%nh   = self%box / 2 + 1
        self%nyq  = self%box / 2
        self%kb   = min(self%nyq - 1, nint(real(self%box) * params%smpd_crop / EM_LIM_LS))
        self%kres = nint(real(self%box) * params%smpd_crop / FSC_RES)
        call se%new(trim(params%pgrp))
        self%nsym = se%get_nsym()
        allocate(self%eo(self%n), self%st(self%n), self%vbin(self%n), self%rot(3,3,self%nsym,self%n), self%dfm(self%n))
        do i = 1, NVIEW
            t  = real(i-1) + 0.5
            z  = 1. - 2. * t / real(NVIEW)
            rz = sqrt(max(1. - z*z, 0.))
            g  = PI * (3. - sqrt(5.)) * t
            dirs(:,i) = [rz*cos(g), rz*sin(g), z]
        end do
        call o%new(is_ptcl=.true.)
        do i = 1, self%n
            ip = self%rows(i)
            self%eo(i)  = build%spproj_field%get_eo(ip)
            self%st(i)  = build%spproj_field%get_state(ip)
            self%dfm(i) = 0.5 * (build%spproj_field%get_dfx(ip) + build%spproj_field%get_dfy(ip))
            call o%set_euler(build%spproj_field%get_euler(ip))
            self%rot(:,:,1,i) = o%get_mat()
            do is = 2, self%nsym
                call se%apply(build%spproj_field%get_euler(ip), is, es)
                call o%set_euler(es)
                self%rot(:,:,is,i) = o%get_mat()
            end do
            e1 = build%spproj_field%e1get(ip) * PI / 180.
            e2 = build%spproj_field%e2get(ip) * PI / 180.
            v  = [sin(e2)*cos(e1), sin(e2)*sin(e1), cos(e2)]
            self%vbin(i) = maxloc(matmul(v, dirs), dim=1)
        end do
        call o%kill
        call se%kill
        self%nst = maxval(self%st)
        write(logfhandle,'(A)') '>>> LANDSCAPE_SAMPLE loading '//int2str(self%n)//' particles of '//int2str(self%nst)//&
            &' states at box_crop '//int2str(self%box)//' (Nyquist '//fmtr(real(2.*params%smpd_crop,dp),1)//' A)'
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
            self%np = 0
            do k = lbound(self%fpls(1)%cmplx_plane,2), 0
                do h = lbound(self%fpls(1)%cmplx_plane,1), ubound(self%fpls(1)%cmplx_plane,1)
                    r2 = h*h + k*k
                    if( r2 > self%nyq*(self%nyq-1) ) cycle
                    self%np = self%np + 1
                    if( ie == 1 ) cycle
                    self%hk(:,self%np) = [h, k]
                    self%psh(self%np)  = nint(sqrt(real(r2)))
                    self%ppar(self%np) = modulo(h + k, 2)
                    self%pw(self%np)   = merge(1., 2., k == 0)
                end do
            end do
            if( ie == 1 ) allocate(self%hk(2,self%np), self%psh(self%np), self%ppar(self%np), self%pw(self%np))
        end do
        allocate(self%cw(0:self%nyq,0:1), self%f2(0:self%nyq,0:1,self%n), self%nw(0:self%nyq,0:1,0:1), source=0.)
        do j = 1, self%np
            self%cw(self%psh(j),self%ppar(j)) = self%cw(self%psh(j),self%ppar(j)) + self%pw(j)
        end do
        pws = [sum(self%pw), sum(self%pw, mask=self%ppar == 0), sum(self%pw, mask=self%ppar == 1)]
        self%pfrac = pws / pws(-1)
        !$omp parallel do default(shared) private(i,j,h,k) schedule(static) proc_bind(close)
        do i = 1, self%n
            if( sdn(i) > 0. ) self%fpls(i)%cmplx_plane = sdn(i) * self%fpls(i)%cmplx_plane
            do j = 1, self%np
                h = self%hk(1,j)
                k = self%hk(2,j)
                self%f2(self%psh(j),self%ppar(j),i) = self%f2(self%psh(j),self%ppar(j),i) + &
                    &self%pw(j) * abs(self%fpls(i)%cmplx_plane(h,k))**2
            end do
            self%fpls(i)%cmplx_plane = real(self%fpls(i)%transfer_plane) * self%fpls(i)%cmplx_plane
            deallocate(self%fpls(i)%transfer_plane)
        end do
        !$omp end parallel do
        do ie = 0, 1
            do q = 0, 1
                do sh = 0, self%nyq
                    t = sum(self%f2(sh,q,:), mask=self%eo == ie)
                    if( t > 0. ) self%nw(sh,q,ie) = self%cw(sh,q) * real(count(self%eo == ie)) / t
                end do
            end do
        end do
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
        allocate(den(self%nh,self%box,self%box), source=0._dp)
        !$omp parallel do default(shared) private(i,j,is,c1,c2,cnj) schedule(static) reduction(+:den) proc_bind(close)
        do i = 1, self%n
            do j = 1, self%np
                do is = 1, self%nsym
                    call self%cell_of(i, is, self%hk(1,j), self%hk(2,j), c1, c2, cnj)
                    if( c1(1) == 0 ) cycle
                    if( c2(1) > 0 )then
                        den(c1(1),c1(2),c1(3)) = den(c1(1),c1(2),c1(3)) + 0.5_dp * self%pw(j) * &
                            &self%fpls(i)%ctfsq_plane(self%hk(1,j),self%hk(2,j))
                        den(c2(1),c2(2),c2(3)) = den(c2(1),c2(2),c2(3)) + 0.5_dp * self%pw(j) * &
                            &self%fpls(i)%ctfsq_plane(self%hk(1,j),self%hk(2,j))
                    else
                        den(c1(1),c1(2),c1(3)) = den(c1(1),c1(2),c1(3)) + self%pw(j) * &
                            &self%fpls(i)%ctfsq_plane(self%hk(1,j),self%hk(2,j))
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        allocate(self%epsv(0:self%nyq), source=0._dp)
        do sh = 0, self%nyq
            vals = real(pack(den, self%fsh == sh .and. den > 0._dp))
            cnt  = size(vals)
            if( cnt == 0 ) cycle
            call hpsort(vals)
            self%epsv(sh) = real(EPS_REG,dp) * real(vals((cnt + 1) / 2), dp) / real(self%n, dp)
        end do
        deallocate(den)
        allocate(ones(self%box,self%box,self%box), source=1.)
        call mskimg%new([self%box,self%box,self%box], params%smpd_crop, wthreads=.false.)
        call mskimg%set_rmat(ones, .false.)
        call mskimg%mask3D_soft(params%msk_crop, backgr=0.)
        self%msk3 = mskimg%get_rmat()
        call mskimg%kill
        call self%wimg%new([self%box,self%box,self%box], params%smpd_crop, wthreads=.false.)
        write(logfhandle,'(A)') '>>> LANDSCAPE_SAMPLE score band shells 2..'//int2str(self%kb)//' ('//&
            &fmtr(real(real(self%box)*params%smpd_crop/real(self%kb),dp),1)//' A), maps to shell '//int2str(self%nyq-1)//&
            &', cell FSC at shell '//int2str(self%kres)//', '//int2str(self%np)//' plane samples'
        self%l_exists = .true.
    end subroutine new

    subroutine run( self, polefile, splitfile, mapdir )
        class(landscape_sample), intent(inout) :: self
        class(string),           intent(in)    :: polefile, splitfile, mapdir
        character(len=1024)  :: line
        integer, allocatable :: child(:), cnts(:), grp(:), rep(:), pa(:), pb(:), mem(:)
        logical, allocatable :: is_rest(:), keep(:), intree(:)
        real,    allocatable :: dist(:,:), dmin(:)
        integer :: funit, ios, nl, j, s, i, ng, g, a, b, npair, ip, jbest, iline
        logical :: l_split
        call self%read_refs(mapdir)
        nl = nlines(polefile)
        allocate(is_rest(self%nst), source=.false.)
        call fopen(funit, polefile, status='old', action='read', iostat=ios)
        call fileiochk('landscape_sample: cannot read '//polefile%to_char(), ios)
        do iline = 1, nl
            read(funit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            read(line,*,iostat=ios) j
            if( ios /= 0 .or. j < 1 .or. j > self%nst ) cycle
            is_rest(j) = index(' '//trim(line)//' ', ' rest ') > 0
        end do
        call fclose(funit)
        if( nl /= self%nst ) write(logfhandle,'(A)') '>>> LANDSCAPE_SAMPLE warning: '//int2str(nl)//' pole lines for '//&
            &int2str(self%nst)//' states'
        allocate(keep(self%n), source=.true.)
        allocate(grp(self%nst), source=0)
        l_split = splitfile%strlen() > 0
        if( l_split )then
            if( nlines(splitfile) /= self%nrows ) THROW_HARD('split labels do not match the project rows')
            allocate(child(self%nrows))
            call fopen(funit, splitfile, status='old', action='read', iostat=ios)
            call fileiochk('landscape_sample: cannot read '//splitfile%to_char(), ios)
            do i = 1, self%nrows
                read(funit,*) child(i)
            end do
            call fclose(funit)
            keep = child(self%rows) > 0
            allocate(cnts(maxval(child)+1))
            do s = 1, self%nst
                if( is_rest(s) ) cycle
                cnts = 0
                do i = 1, self%n
                    if( self%st(i) == s .and. keep(i) ) cnts(child(self%rows(i))) = cnts(child(self%rows(i))) + 1
                end do
                if( maxval(cnts) > 0 ) grp(s) = maxloc(cnts, dim=1)
            end do
        else
            do s = 1, self%nst
                if( .not. is_rest(s) ) grp(s) = s
            end do
        endif
        ng  = 0
        allocate(rep(self%nst), source=0)
        do s = 1, self%nst
            if( grp(s) == 0 ) cycle
            if( any(grp(1:s-1) == grp(s)) ) cycle
            ng = ng + 1
            rep(ng) = s
            do j = s + 1, self%nst
                if( grp(j) == grp(s) .and. count(self%st == j) > count(self%st == rep(ng)) ) rep(ng) = j
            end do
        end do
        write(logfhandle,'(A)') '>>> LANDSCAPE_SAMPLE '//int2str(ng)//' accepted poles (states '//joined(rep(1:ng))//&
            &'); rest states '//joined(pack([(s, s=1,self%nst)], is_rest))
        call fopen(funit, string('landscape_pairs.txt'), status='replace', action='write', iostat=ios)
        call fileiochk('landscape_sample: cannot write landscape_pairs.txt', ios)
        if( ng < 2 )then
            write(logfhandle,'(A)') '>>> LANDSCAPE_SAMPLE fewer than two accepted poles: no pairs'
            call fclose(funit)
            return
        endif
        call self%score
        if( ng <= 3 )then
            npair = ng * (ng - 1) / 2
            allocate(pa(npair), pb(npair))
            ip = 0
            do a = 1, ng
                do b = a + 1, ng
                    ip = ip + 1
                    pa(ip) = a
                    pb(ip) = b
                end do
            end do
        else
            allocate(dist(ng,ng), source=0.)
            do a = 1, ng
                do b = a + 1, ng
                    dist(a,b) = map_dist(rep(a), rep(b))
                    dist(b,a) = dist(a,b)
                end do
            end do
            npair = ng - 1
            allocate(pa(npair), pb(npair), intree(ng), dmin(ng), mem(ng))
            intree = .false.
            intree(1) = .true.
            dmin = dist(1,:)
            mem  = 1
            do ip = 1, npair
                jbest = minloc(dmin, mask=.not. intree, dim=1)
                pa(ip) = mem(jbest)
                pb(ip) = jbest
                intree(jbest) = .true.
                do g = 1, ng
                    if( .not. intree(g) .and. dist(jbest,g) < dmin(g) )then
                        dmin(g) = dist(jbest,g)
                        mem(g)  = jbest
                    endif
                end do
            end do
        endif
        do ip = 1, npair
            call self%sample_pair(rep(pa(ip)), rep(pb(ip)), members(pa(ip)), members(pb(ip)), is_rest, keep, funit)
        end do
        call fclose(funit)

    contains

        function members( ig ) result( ms )
            integer, intent(in)  :: ig
            integer, allocatable :: ms(:)
            ms = pack([(s, s=1,self%nst)], grp == grp(rep(ig)))
        end function members

        real function map_dist( s1, s2 )
            integer, intent(in) :: s1, s2
            integer :: hh, kk, ll
            map_dist = 0.
            do ll = 1, self%box
                do kk = 1, self%box
                    do hh = 1, self%nh
                        if( self%fsh(hh,kk,ll) < 2 .or. self%fsh(hh,kk,ll) > self%kb ) cycle
                        map_dist = map_dist + self%fw(hh,kk,ll) * abs(0.5 * (self%refs(hh,kk,ll,1,s1) + &
                            &self%refs(hh,kk,ll,2,s1) - self%refs(hh,kk,ll,1,s2) - self%refs(hh,kk,ll,2,s2)))**2
                    end do
                end do
            end do
        end function map_dist

    end subroutine run

    subroutine read_refs( self, mapdir )
        use simple_refine3D_fnames, only: refine3D_state_halfvol_fname
        class(landscape_sample), intent(inout) :: self
        class(string),           intent(in)    :: mapdir
        character(len=4), parameter :: HALVES(2) = ['even', 'odd ']
        type(image)  :: vol
        type(string) :: fname
        integer :: s, ih
        allocate(self%refs(self%nh,self%box,self%box,2,self%nst))
        do s = 1, self%nst
            do ih = 1, 2
                fname = mapdir//'/'//refine3D_state_halfvol_fname(s, trim(HALVES(ih)), unfil=.true.)
                call vol%read_and_crop(fname, self%p%smpd, self%box, self%p%smpd_crop)
                call vol%fft
                self%refs(:,:,:,ih,s) = vol%get_cmat()
                call vol%kill
            end do
        end do
    end subroutine read_refs

    subroutine score( self )
        class(landscape_sample), intent(inout) :: self
        real(dp) :: sa(0:self%nyq), sb(0:self%nyq), sf(0:self%nyq), scn(0:self%nyq), beta(0:self%nyq), sig2(0:self%nyq)
        real(dp) :: ao(0:self%nyq), bo(0:self%nyq), s(0:1), x
        complex  :: pv, dv
        real     :: ws
        integer  :: i, j, g, ih, sh, q, c1(3), c2(3), is
        logical  :: cnj
        sa  = 0._dp
        sb  = 0._dp
        sf  = 0._dp
        scn = 0._dp
        !$omp parallel do default(shared) private(i,j,g,ih,sh,c1,c2,cnj,pv,dv,ws,is) schedule(static) &
        !$omp& reduction(+:sa,sb,sf,scn) proc_bind(close)
        do i = 1, self%n
            g  = self%st(i)
            ih = 2 - self%eo(i)
            do j = 1, self%np
                sh = self%psh(j)
                if( sh < 2 .or. sh > self%kb ) cycle
                dv = self%fpls(i)%cmplx_plane(self%hk(1,j),self%hk(2,j))
                ws = self%pw(j) / real(self%nsym)
                do is = 1, self%nsym
                    call self%cell_of(i, is, self%hk(1,j), self%hk(2,j), c1, c2, cnj)
                    if( c1(1) == 0 ) cycle
                    pv = self%refs(c1(1),c1(2),c1(3),ih,g)
                    if( cnj ) pv = conjg(pv)
                    sa(sh) = sa(sh) + ws * (real(pv)*real(dv) + aimag(pv)*aimag(dv))
                    sb(sh) = sb(sh) + ws * self%fpls(i)%ctfsq_plane(self%hk(1,j),self%hk(2,j)) * abs(pv)**2
                end do
            end do
            sf(2:self%kb)  = sf(2:self%kb)  + real(self%f2(2:self%kb,0,i) + self%f2(2:self%kb,1,i), dp)
            scn(2:self%kb) = scn(2:self%kb) + real(self%cw(2:self%kb,0) + self%cw(2:self%kb,1), dp)
        end do
        !$omp end parallel do
        beta = min(max(sa / max(sb, 1.e-30_dp), 0._dp), 4._dp)
        sig2 = max((sf - 2._dp*beta*sa + beta*beta*sb) / max(scn, 1._dp), 1.e-12_dp)
        allocate(self%sc(0:1,self%nst,self%n), self%ppn(self%n), source=0.)
        !$omp parallel do default(shared) private(i,j,g,ih,sh,q,c1,c2,cnj,pv,dv,s,ao,bo,x,ws,is) schedule(static) &
        !$omp& proc_bind(close)
        do i = 1, self%n
            ih = 2 - self%eo(i)
            do g = 1, self%nst
                s  = 0._dp
                ao = 0._dp
                bo = 0._dp
                do j = 1, self%np
                    sh = self%psh(j)
                    if( sh < 2 .or. sh > self%kb ) cycle
                    dv = self%fpls(i)%cmplx_plane(self%hk(1,j),self%hk(2,j))
                    q  = self%ppar(j)
                    ws = self%pw(j) / real(self%nsym)
                    do is = 1, self%nsym
                        call self%cell_of(i, is, self%hk(1,j), self%hk(2,j), c1, c2, cnj)
                        if( c1(1) == 0 ) cycle
                        pv = self%refs(c1(1),c1(2),c1(3),ih,g)
                        if( cnj ) pv = conjg(pv)
                        x  = real(ws * (real(pv)*real(dv) + aimag(pv)*aimag(dv)), dp)
                        ao(sh) = ao(sh) + x
                        bo(sh) = bo(sh) + real(ws * self%fpls(i)%ctfsq_plane(self%hk(1,j),self%hk(2,j)) * abs(pv)**2, dp)
                        s(q)   = s(q) + 0.5_dp * (2._dp * beta(sh) * x - beta(sh)**2 * real(ws * &
                            &self%fpls(i)%ctfsq_plane(self%hk(1,j),self%hk(2,j)) * abs(pv)**2, dp)) / sig2(sh)
                    end do
                end do
                self%sc(:,g,i) = real(s)
                if( g == self%st(i) ) self%ppn(i) = real(sum((real(self%f2(2:self%kb,0,i) + self%f2(2:self%kb,1,i), dp) - &
                    &2._dp*beta(2:self%kb)*ao(2:self%kb) + beta(2:self%kb)**2*bo(2:self%kb)) / sig2(2:self%kb)) / &
                    &sum(real(self%cw(2:self%kb,0) + self%cw(2:self%kb,1), dp)))
            end do
        end do
        !$omp end parallel do
        where( self%ppn <= 0. ) self%ppn = 1.
        write(logfhandle,'(A)') '>>> LANDSCAPE_SAMPLE scores against '//int2str(self%nst)//' states; model scale at shells 2 / '//&
            &int2str(self%kb)//' '//fmtr(beta(2),2)//' / '//fmtr(beta(self%kb),2)//'; particle noise levels 5/50/95 % '//&
            &fmtr(real(quantr(self%ppn,0.05),dp),2)//'/'//fmtr(real(quantr(self%ppn,0.5),dp),2)//'/'//&
            &fmtr(real(quantr(self%ppn,0.95),dp),2)
    end subroutine score

    subroutine sample_pair( self, sa, sb, mema, memb, is_rest, keep, funit )
        use simple_syslib, only: simple_mkdir
        class(landscape_sample), intent(inout) :: self
        integer,                 intent(in)    :: sa, sb, mema(:), memb(:), funit
        logical,                 intent(in)    :: is_rest(:), keep(:)
        character(len=1), parameter :: FOLDS(2) = ['A', 'B']
        type(fold_res), allocatable :: res(:,:)
        type(fold_res)       :: fin
        complex, allocatable :: cons(:,:,:,:,:), m0(:,:,:), m1(:,:,:)
        real,    allocatable :: t(:,:), cfull(:), cent(:)
        integer, allocatable :: sel(:), ord(:), sid(:), out(:), lab(:)
        type(string) :: tag
        real    :: f(0:self%nyq), fctr(2,2), fcf(NQB,2), ftt(NQB,2), lv, res143, mu, sd
        integer :: best(2), i, k, ifold, ic, nf, nn, nbest, q, ih, m, funit2, ios, nsel, kk
        sel  = pack([(i, i=1,self%n)], [(keep(i) .and. (any(mema == self%st(i)) .or. any(memb == self%st(i)) .or. &
            &is_rest(self%st(i))), i=1,self%n)])
        nsel = size(sel)
        tag  = string('landscape_'//int2str(sa)//'_'//int2str(sb))
        allocate(t(nsel,3))
        do q = 0, 1
            t(:,q+1) = (self%sc(q,sb,sel) - self%sc(q,sa,sel)) / self%ppn(sel)
        end do
        t(:,3) = t(:,1) + t(:,2)
        do k = 1, 3
            mu = sum(t(:,k)) / real(nsel)
            sd = sqrt(sum((t(:,k) - mu)**2) / real(nsel))
            t(:,k) = (t(:,k) - mu) / max(sd, 1.e-12)
        end do
        write(logfhandle,'(A)') '>>> LANDSCAPE_SAMPLE pair '//int2str(sa)//' - '//int2str(sb)//': '//int2str(nsel)//&
            &' particles (states '//joined(pack([(k, k=1,self%nst)], [(any(mema == k) .or. any(memb == k) .or. is_rest(k), &
            &k=1,self%nst)]))//'); fold A vs B coordinate corr '//fmtr(real(pcorr(t(:,1), t(:,2)),dp),3)
        allocate(cons(self%nh,self%box,self%box,0:1,0:1))
        do q = 0, 1
            do ih = 0, 1
                call self%nnmap(pack(sel, self%eo(sel) == ih), q, cons(:,:,:,ih,q))
            end do
        end do
        allocate(res(NCAND,2))
        best = 0
        do ifold = 1, 2
            q  = ifold - 1
            nf = 0
            write(logfhandle,'(A)') '  fold '//FOLDS(ifold)//' (coordinate from pixel set '//FOLDS(ifold)//', tests on set '//&
                &FOLDS(3-ifold)//'):'
            call self%struct_centres(sel, t(:,ifold), 1 - q, cons(:,:,:,:,1-q), fctr(:,ifold), fcf(:,ifold), ftt(:,ifold))
            do ic = 1, NCAND
                nn   = CANDS(ic)
                cent = centres(nn, fctr(:,ifold), fcf(:,ifold), ftt(:,ifold))
                call self%evaluate(sel, t(:,ifold), nn, cent, 1 - q, cons(:,:,:,:,1-q), res(ic,ifold))
                if( res(ic,ifold)%ok )then
                    best(ifold) = nn
                    nf = 0
                else
                    nf = nf + 1
                endif
                if( best(ifold) > 0 .and. nf >= 2 ) exit
            end do
        end do
        write(logfhandle,'(A)') '  resolvable samples: fold A '//int2str(best(1))//', fold B '//int2str(best(2))
        if( minval(best) == 0 )then
            write(logfhandle,'(A)') '  a fold resolves no candidate: no landscape for this pair'
            write(funit,'(I0,1X,I0,1X,I0,1X,I0,1X,I0)') sa, sb, 0, best(1), best(2)
            return
        endif
        nbest = minval(best)
        allocate(cfull(nbest))
        do k = 1, nbest
            lv = 0.
            do ifold = 1, 2
                cent = centres(nbest, fctr(:,ifold), fcf(:,ifold), ftt(:,ifold))
                lv   = lv + 0.5 * real(count(t(:,ifold) < cent(k))) / real(nsel)
            end do
            cfull(k) = quantr(t(:,3), min(max(lv, 0.), 1.))
        end do
        fin%grid = make_cells(t(:,3), nbest, cfull)
        fin%cell = fin%grid
        call self%merge_sparse(sel, t(:,3), -1, fin%cell, fin%fs)
        call ordered_cells(fin%cell, t(:,3), ord)
        m = size(ord)
        allocate(sid(nbest), source=0)
        do k = 1, m
            sid(ord(k)) = k
        end do
        allocate(out(self%nrows), lab(nsel), source=0)
        lab = sid(fin%cell)
        out(self%rows(sel)) = lab
        call write_ints(tag//'_samples.txt', out)
        call fopen(funit2, tag//'_coords.txt', status='replace', action='write', iostat=ios)
        call fileiochk('landscape_sample: cannot write '//tag%to_char()//'_coords.txt', ios)
        do i = 1, nsel
            write(funit2,'(I0,3(1X,F9.4),1X,I0)') self%rows(sel(i)), t(i,3), t(i,1), t(i,2), lab(i)
        end do
        call fclose(funit2)
        write(logfhandle,'(A)') '  CHOSEN '//int2str(nbest)//' -> '//int2str(m)//' samples ('//int2str(nbest-m)//&
            &' merged as not reproducible); sizes '//joined([(count(lab == k), k=1,m)])
        call simple_mkdir(tag)
        call simple_mkdir(tag//'/final')
        allocate(m0(self%nh,self%box,self%box), m1(self%nh,self%box,self%box))
        do k = 1, m
            call self%nnmap(pack(sel, lab == k .and. self%eo(sel) == 0), -1, m0)
            call self%nnmap(pack(sel, lab == k .and. self%eo(sel) == 1), -1, m1)
            call self%shell_fsc(m0, m1, f)
            res143 = real(self%box) * self%p%smpd_crop / real(self%nyq - 1)
            do kk = 2, self%nyq - 1
                if( f(kk) < 0.143 )then
                    res143 = real(self%box) * self%p%smpd_crop / real(kk)
                    exit
                endif
            end do
            write(logfhandle,'(A)') '  sample '//int2str(k)//': '//int2str(count(lab == k))//' particles, centre '//&
                &fmtr(real(meanm(t(:,3), lab == k),dp),2)//' sd, FSC 0.143 at '//fmtr(real(res143,dp),1)//' A'
            call self%write_map(0.5 * (m0 + m1), tag//'/final/state_'//pad2(k)//'.mrc')
        end do
        ic = findloc(CANDS, nbest, dim=1)
        do ifold = 1, 2
            associate( r => res(ic,ifold) )
                write(logfhandle,'(A)') '  fold '//FOLDS(ifold)//' at n '//int2str(nbest)//': steps (end-to-end units) '//&
                    &joinedr(r%st(1,1:r%npr))//'; z '//joinedr(r%st(2,1:r%npr))//' (threshold '//fmtr(real(r%zc,dp),2)//&
                    &'); steps CV '//fmtr(real(cvof(r%st(1,1:r%npr)),dp),2)
                write(logfhandle,'(A)') '  check 3, fold '//FOLDS(ifold)//' (test-pixel positions, end cells at 0 / 1): end_hi '//&
                    &fmtr(real(meanm(r%u, r%grid == nbest) - meanm(r%u, r%grid == 1),dp),2)//', tail_lo '//&
                    &fmtr(real(meanm(r%u, t(:,ifold) < quantr(t(:,ifold), QLO)) - meanm(r%u, r%grid == 1),dp),2)//&
                    &', tail_hi '//fmtr(real(meanm(r%u, t(:,ifold) > quantr(t(:,ifold), 1. - QLO)) - &
                    &meanm(r%u, r%grid == 1),dp),2)
                call simple_mkdir(tag//'/heldout_'//FOLDS(ifold))
                call ordered_cells(r%cell, t(:,ifold), ord)
                do k = 1, size(ord)
                    call self%nnmap(pack(sel, r%cell == ord(k) .and. self%eo(sel) == 0), 2 - ifold, m0)
                    call self%nnmap(pack(sel, r%cell == ord(k) .and. self%eo(sel) == 1), 2 - ifold, m1)
                    call self%write_map(0.5 * (m0 + m1), tag//'/heldout_'//FOLDS(ifold)//'/state_'//pad2(k)//'.mrc')
                end do
            end associate
        end do
        write(funit,'(I0,1X,I0,1X,I0,1X,I0,1X,I0,2(1X,F6.3))') sa, sb, m, best(1), best(2), &
            &cvof(res(ic,1)%st(1,1:res(ic,1)%npr)), cvof(res(ic,2)%st(1,1:res(ic,2)%npr))
    end subroutine sample_pair

    subroutine struct_centres( self, sel, tc, q, cons, ctr, cf, tt )
        class(landscape_sample), intent(inout) :: self
        integer,                 intent(in)    :: sel(:), q
        real,                    intent(in)    :: tc(:)
        complex,                 intent(in)    :: cons(:,:,:,0:)
        real,                    intent(out)   :: ctr(2), cf(NQB), tt(NQB)
        integer, allocatable :: grid(:), b(:)
        real,    allocatable :: u(:)
        real    :: um(NQB), w(NQB), qb(0:NQB)
        integer :: k, i
        grid = make_cells(tc, 4)
        call self%ustat(sel, grid, 4, q, cons, u)
        do k = 0, NQB
            qb(k) = quantr(tc, real(k) / real(NQB))
        end do
        allocate(b(size(tc)))
        do i = 1, size(tc)
            b(i) = min(max(count(qb(0:NQB-1) <= tc(i)), 1), NQB)
        end do
        do k = 1, NQB
            w(k)  = real(count(b == k))
            tt(k) = meanm(tc, b == k)
            um(k) = meanm(u, b == k)
        end do
        cf  = pava(um, max(w, 1.)) + 1.e-6 * real([(k, k=0,NQB-1)])
        ctr = [quantr(tc, QLO), quantr(tc, 1. - QLO)]
        write(logfhandle,'(A)') '    structural position E[u | t], 10 quantile steps: '//&
            &joinedr([(interp(quantr(tc, 0.05 + 0.1 * real(k)), tt, cf), k=0,9)])
    end subroutine struct_centres

    subroutine evaluate( self, sel, tc, nn, cent, q, cons, r )
        class(landscape_sample), intent(inout) :: self
        integer,                 intent(in)    :: sel(:), nn, q
        real,                    intent(in)    :: tc(:), cent(:)
        complex,                 intent(in)    :: cons(:,:,:,0:)
        type(fold_res),          intent(inout) :: r
        integer, allocatable :: occ(:), prs(:,:)
        integer :: k, ip, p1, p2, ih, np_
        real(dp) :: ds(0:1), vs(0:1), wts(0:1), d
        real, allocatable :: up(:), uq(:)
        r%grid = make_cells(tc, nn, cent)
        r%cell = r%grid
        call self%merge_sparse(sel, tc, q, r%cell, r%fs)
        allocate(occ(nn), source=0)
        do k = 1, nn
            if( any(r%grid == k) ) occ(k) = maxval(r%cell, mask=r%grid == k)
        end do
        allocate(prs(2,0))
        do k = 1, nn - 1
            if( occ(k) == 0 .or. occ(k+1) == 0 .or. occ(k) == occ(k+1) ) cycle
            if( size(prs,2) > 0 )then
                if( any(prs(1,:) == occ(k) .and. prs(2,:) == occ(k+1)) ) cycle
            endif
            prs = reshape([prs, occ(k), occ(k+1)], [2, size(prs,2)+1])
        end do
        np_   = size(prs,2)
        r%npr = np_
        r%prs = prs
        r%ok  = .false.
        if( np_ < 1 ) return
        call self%ustat(sel, r%grid, nn, q, cons, r%u)
        allocate(r%st(4,np_), source=0.)
        r%zc = norm_ppf(1. - ALPHA / real(np_))
        r%ok = .true.
        do ip = 1, np_
            p1 = prs(1,ip)
            p2 = prs(2,ip)
            do ih = 0, 1
                up = pack(r%u, r%cell == p1 .and. self%eo(sel) == ih)
                uq = pack(r%u, r%cell == p2 .and. self%eo(sel) == ih)
                if( size(up) < 2 .or. size(uq) < 2 )then
                    ds(ih) = 0._dp
                    vs(ih) = 1.e30_dp
                    cycle
                endif
                ds(ih) = real(sum(uq) / size(uq) - sum(up) / size(up), dp)
                vs(ih) = real(varp(up) / size(up) + varp(uq) / size(uq), dp)
            end do
            wts = 1._dp / vs
            d   = sum(wts * ds) / sum(wts)
            r%st(:,ip) = real([d, d * sqrt(sum(wts)), ds(0), ds(1)])
            r%ok = r%ok .and. r%st(2,ip) > r%zc .and. ds(0) > 0._dp .and. ds(1) > 0._dp
        end do
        write(logfhandle,'(A)') '    n '//int2str(nn)//': '//int2str(count([(any(r%cell == k), k=1,nn)]))//' samples, steps '//&
            &joinedr(r%st(1,:))//'; z '//joinedr(r%st(2,:))//' (threshold '//fmtr(real(r%zc,dp),2)//') -> '//&
            &merge('ALL RESOLVED    ', 'not all resolved', r%ok)
    end subroutine evaluate

    subroutine merge_sparse( self, sel, tc, q, cell, fs )
        class(landscape_sample), intent(inout) :: self
        integer,                 intent(in)    :: sel(:), q
        real,                    intent(in)    :: tc(:)
        integer,                 intent(inout) :: cell(:)
        real, allocatable,       intent(out)   :: fs(:)
        integer :: nc, k, g, tgt
        real    :: cg, dmin
        nc = maxval(cell)
        allocate(fs(nc), source=huge(1.))
        do k = 1, nc
            if( any(cell == k) ) fs(k) = self%cellfsc(pack(sel, cell == k), q)
        end do
        do
            if( count([(any(cell == k), k=1,nc)]) < 2 ) exit
            g = 0
            do k = 1, nc
                if( .not. any(cell == k) .or. fs(k) >= MINFSC ) cycle
                if( g == 0 )then
                    g = k
                else if( count(cell == k) < count(cell == g) )then
                    g = k
                endif
            end do
            if( g == 0 ) exit
            cg   = meanm(tc, cell == g)
            tgt  = 0
            dmin = huge(1.)
            do k = 1, nc
                if( k == g .or. .not. any(cell == k) ) cycle
                if( abs(meanm(tc, cell == k) - cg) < dmin )then
                    dmin = abs(meanm(tc, cell == k) - cg)
                    tgt  = k
                endif
            end do
            where( cell == g ) cell = tgt
            fs(g)   = huge(1.)
            fs(tgt) = self%cellfsc(pack(sel, cell == tgt), q)
        end do
    end subroutine merge_sparse

    real function cellfsc( self, idx, q )
        class(landscape_sample), intent(inout) :: self
        integer,                 intent(in)    :: idx(:), q
        complex, allocatable :: m0(:,:,:), m1(:,:,:)
        real :: f(0:self%nyq)
        cellfsc = -1.
        if( min(count(self%eo(idx) == 0), count(self%eo(idx) == 1)) < MINHALF ) return
        allocate(m0(self%nh,self%box,self%box), m1(self%nh,self%box,self%box))
        call self%nnmap(pack(idx, self%eo(idx) == 0), q, m0)
        call self%nnmap(pack(idx, self%eo(idx) == 1), q, m1)
        call self%shell_fsc(m0, m1, f)
        cellfsc = sum(f(max(self%kres-1,2):self%kres+1)) / real(self%kres + 2 - max(self%kres-1,2))
    end function cellfsc

    subroutine ustat( self, sel, grid, nn, q, cons, u )
        class(landscape_sample), intent(inout) :: self
        integer,                 intent(in)    :: sel(:), grid(:), nn, q
        complex,                 intent(in)    :: cons(:,:,:,0:)
        real, allocatable,       intent(out)   :: u(:)
        complex, allocatable :: dm(:,:,:,:), e(:,:,:)
        integer, allocatable :: strat(:), idx(:)
        real     :: qd(2), scl, ws
        real(dp) :: acc
        complex  :: pv, dv
        integer  :: ih, i, j, sh, c1(3), c2(3), k, ii, ns, is
        logical  :: cnj
        allocate(dm(self%nh,self%box,self%box,0:1), e(self%nh,self%box,self%box))
        do ih = 0, 1
            call self%nnmap(pack(sel, grid == nn .and. self%eo(sel) == ih), q, e)
            call self%dock(e, cons(:,:,:,1-ih), 2, self%nyq-1)
            dm(:,:,:,ih) = e
            call self%nnmap(pack(sel, grid == 1 .and. self%eo(sel) == ih), q, e)
            call self%dock(e, cons(:,:,:,1-ih), 2, self%nyq-1)
            dm(:,:,:,ih) = dm(:,:,:,ih) - e
        end do
        allocate(u(size(sel)), source=0.)
        !$omp parallel do default(shared) private(ii,i,ih,j,sh,c1,c2,cnj,pv,dv,acc,ws,is) schedule(static) proc_bind(close)
        do ii = 1, size(sel)
            i   = sel(ii)
            ih  = self%eo(i)
            acc = 0._dp
            do j = 1, self%np
                sh = self%psh(j)
                if( self%ppar(j) /= q .or. sh < 2 ) cycle
                dv = self%fpls(i)%cmplx_plane(self%hk(1,j),self%hk(2,j))
                ws = self%pw(j) / real(self%nsym)
                do is = 1, self%nsym
                    call self%cell_of(i, is, self%hk(1,j), self%hk(2,j), c1, c2, cnj)
                    if( c1(1) == 0 ) cycle
                    pv = dm(c1(1),c1(2),c1(3),1-ih)
                    if( cnj ) pv = conjg(pv)
                    acc = acc + real(ws * self%nw(sh,q,ih) * (real(pv)*real(dv) + aimag(pv)*aimag(dv)), dp)
                end do
            end do
            u(ii) = real(acc)
        end do
        !$omp end parallel do
        do ih = 0, 1
            idx = pack([(ii, ii=1,size(sel))], self%eo(sel) == ih)
            if( size(idx) == 0 ) cycle
            qd(1) = quantr(self%dfm(sel(idx)), 1. / 3.)
            qd(2) = quantr(self%dfm(sel(idx)), 2. / 3.)
            allocate(strat(size(idx)))
            do k = 1, size(idx)
                strat(k) = (self%vbin(sel(idx(k))) - 1) * 3 + count(qd < self%dfm(sel(idx(k))))
            end do
            do ns = 0, 3 * NVIEW - 1
                if( .not. any(strat == ns) ) cycle
                u(idx) = merge(u(idx) - meanm(u(idx), strat == ns), u(idx), strat == ns)
            end do
            deallocate(strat)
            scl = meanm(u(idx), grid(idx) == nn) - meanm(u(idx), grid(idx) == 1)
            if( abs(scl) > 1.e-30 ) u(idx) = u(idx) / scl
        end do
    end subroutine ustat

    subroutine nnmap( self, pidx, q, out )
        class(landscape_sample), intent(inout) :: self
        integer,                 intent(in)    :: pidx(:), q
        complex,                 intent(out)   :: out(:,:,:)
        complex(dp), allocatable :: num(:,:,:)
        real(dp),    allocatable :: den(:,:,:)
        real, pointer :: rmat(:,:,:)
        complex(dp) :: v, vj
        real(dp)    :: c2v, reg
        integer     :: ii, j, c1(3), c2(3), hh, kk, ll, is
        logical     :: cnj
        allocate(num(self%nh,self%box,self%box), source=cmplx(0._dp,0._dp,dp))
        allocate(den(self%nh,self%box,self%box), source=0._dp)
        !$omp parallel do default(shared) private(ii,j,c1,c2,cnj,v,vj,c2v,is) schedule(static) reduction(+:num,den) &
        !$omp& proc_bind(close)
        do ii = 1, size(pidx)
            do j = 1, self%np
                if( q >= 0 .and. self%ppar(j) /= q ) cycle
                vj  = real(self%pw(j), dp) * cmplx(self%fpls(pidx(ii))%cmplx_plane(self%hk(1,j),self%hk(2,j)), kind=dp)
                c2v = real(self%pw(j), dp) * real(self%fpls(pidx(ii))%ctfsq_plane(self%hk(1,j),self%hk(2,j)), dp)
                do is = 1, self%nsym
                    call self%cell_of(pidx(ii), is, self%hk(1,j), self%hk(2,j), c1, c2, cnj)
                    if( c1(1) == 0 ) cycle
                    v = vj
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
        end do
        !$omp end parallel do
        reg = real(self%pfrac(max(q,-1)), dp) * real(max(size(pidx),1), dp)
        do ll = 1, self%box
            do kk = 1, self%box
                do hh = 1, self%nh
                    out(hh,kk,ll) = cmplx(num(hh,kk,ll) / (den(hh,kk,ll) + self%epsv(self%fsh(hh,kk,ll)) * reg), kind=kind(out))
                end do
            end do
        end do
        call self%wimg%set_ft(.true.)
        call self%wimg%set_cmat(out)
        call self%wimg%ifft
        call self%wimg%get_rmat_ptr(rmat)
        rmat = rmat * self%msk3
        call self%wimg%fft
        out = self%wimg%get_cmat()
    end subroutine nnmap

    subroutine cell_of( self, i, is, h, k, c1, c2, cnj )
        class(landscape_sample), intent(in)  :: self
        integer,                 intent(in)  :: i, is, h, k
        integer,                 intent(out) :: c1(3), c2(3)
        logical,                 intent(out) :: cnj
        integer :: p(3)
        p   = nint(real(h) * self%rot(1,:,is,i) + real(k) * self%rot(2,:,is,i))
        cnj = p(1) < 0
        if( cnj ) p = -p
        c1  = 0
        c2  = 0
        if( any(abs(p) >= self%nyq) ) return
        c1 = [p(1) + 1, modulo(p(2), self%box) + 1, modulo(p(3), self%box) + 1]
        if( p(1) == 0 ) c2 = [1, modulo(-p(2), self%box) + 1, modulo(-p(3), self%box) + 1]
    end subroutine cell_of

    subroutine dock( self, cm, cref, b1, b2 )
        class(landscape_sample), intent(in)    :: self
        complex,                 intent(inout) :: cm(:,:,:)
        complex,                 intent(in)    :: cref(:,:,:)
        integer,                 intent(in)    :: b1, b2
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
        class(landscape_sample), intent(in)  :: self
        complex,                 intent(in)  :: x(:,:,:), y(:,:,:)
        real,                    intent(out) :: f(0:)
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

    subroutine write_map( self, cm, fname )
        class(landscape_sample), intent(inout) :: self
        complex,                 intent(in)    :: cm(:,:,:)
        class(string),           intent(in)    :: fname
        call self%wimg%set_ft(.true.)
        call self%wimg%set_cmat(cm)
        call self%wimg%ifft
        call self%wimg%write(fname)
    end subroutine write_map

    subroutine kill( self )
        class(landscape_sample), intent(inout) :: self
        if( allocated(self%fpls) ) deallocate(self%fpls)
        call self%wimg%kill
        if( allocated(self%rows) ) deallocate(self%rows, self%eo, self%st, self%vbin, self%rot, self%dfm)
        if( allocated(self%hk)   ) deallocate(self%hk, self%psh, self%ppar, self%pw)
        if( allocated(self%cw)   ) deallocate(self%cw, self%f2, self%nw)
        if( allocated(self%fsh)  ) deallocate(self%fsh, self%fw, self%fqh, self%fqk)
        if( allocated(self%epsv) ) deallocate(self%epsv)
        if( allocated(self%msk3) ) deallocate(self%msk3)
        if( allocated(self%refs) ) deallocate(self%refs)
        if( allocated(self%sc)   ) deallocate(self%sc, self%ppn)
        self%p => null()
        self%b => null()
        self%n = 0
        self%l_exists = .false.
    end subroutine kill

    function make_cells( tc, nn, cent ) result( cell )
        real,           intent(in) :: tc(:)
        integer,        intent(in) :: nn
        real, optional, intent(in) :: cent(:)
        integer, allocatable :: cell(:)
        real    :: lo, hi, w
        integer :: i
        allocate(cell(size(tc)))
        if( present(cent) )then
            do i = 1, size(tc)
                cell(i) = 1 + count(0.5 * (cent(1:nn-1) + cent(2:nn)) < tc(i))
            end do
        else
            lo = quantr(tc, QLO)
            hi = quantr(tc, 1. - QLO)
            w  = (hi - lo) / real(nn - 1)
            do i = 1, size(tc)
                cell(i) = 1 + min(max(nint((tc(i) - lo) / w), 0), nn - 1)
            end do
        endif
    end function make_cells

    function centres( nn, ctr, cf, tt ) result( c )
        integer, intent(in) :: nn
        real,    intent(in) :: ctr(2), cf(:), tt(:)
        real, allocatable :: c(:)
        real    :: flo, fhi
        integer :: k
        flo = interp(ctr(1), tt, cf)
        fhi = interp(ctr(2), tt, cf)
        allocate(c(nn))
        c(1)  = ctr(1)
        c(nn) = ctr(2)
        do k = 1, nn - 2
            c(k+1) = interp(flo + (fhi - flo) * real(k) / real(nn - 1), cf, tt)
        end do
    end function centres

    function pava( y, w ) result( f )
        real, intent(in) :: y(:), w(:)
        real, allocatable :: f(:), v(:), wt(:)
        integer, allocatable :: cnt(:)
        integer :: i, m, k, j
        allocate(v(size(y)), wt(size(y)), cnt(size(y)), f(size(y)))
        m = 0
        do i = 1, size(y)
            m = m + 1
            v(m)   = y(i)
            wt(m)  = w(i)
            cnt(m) = 1
            do while( m > 1 )
                if( v(m-1) <= v(m) ) exit
                v(m-1)   = (v(m-1) * wt(m-1) + v(m) * wt(m)) / (wt(m-1) + wt(m))
                wt(m-1)  = wt(m-1) + wt(m)
                cnt(m-1) = cnt(m-1) + cnt(m)
                m = m - 1
            end do
        end do
        j = 0
        do k = 1, m
            f(j+1:j+cnt(k)) = v(k)
            j = j + cnt(k)
        end do
    end function pava

    real function interp( x, xp, fp )
        real, intent(in) :: x, xp(:), fp(:)
        integer :: i
        if( x <= xp(1) )then
            interp = fp(1)
        else if( x >= xp(size(xp)) )then
            interp = fp(size(fp))
        else
            i = count(xp <= x)
            interp = fp(i) + (fp(i+1) - fp(i)) * (x - xp(i)) / max(xp(i+1) - xp(i), tiny(x))
        endif
    end function interp

    real function norm_ppf( p )
        real, intent(in) :: p
        real(dp), parameter :: A(6) = [-3.969683028665376e+01_dp, 2.209460984245205e+02_dp, -2.759285104469687e+02_dp, &
            &1.383577518672690e+02_dp, -3.066479806614716e+01_dp, 2.506628277459239e+00_dp]
        real(dp), parameter :: B(5) = [-5.447609879822406e+01_dp, 1.615858368580409e+02_dp, -1.556989798598866e+02_dp, &
            &6.680131188771972e+01_dp, -1.328068155288572e+01_dp]
        real(dp), parameter :: C(6) = [-7.784894002430293e-03_dp, -3.223964580411365e-01_dp, -2.400758277161838e+00_dp, &
            &-2.549732539343734e+00_dp, 4.374664141464968e+00_dp, 2.938163982698783e+00_dp]
        real(dp), parameter :: D(4) = [7.784695709041462e-03_dp, 3.224671290700398e-01_dp, 2.445134137142996e+00_dp, &
            &3.754408661907416e+00_dp]
        real(dp) :: q, r, x
        x = real(min(max(p, 1.e-12), 1. - 1.e-12), dp)
        if( x < 0.02425_dp )then
            q = sqrt(-2._dp * log(x))
            norm_ppf = real((((((C(1)*q+C(2))*q+C(3))*q+C(4))*q+C(5))*q+C(6)) / ((((D(1)*q+D(2))*q+D(3))*q+D(4))*q+1._dp))
        else if( x > 1._dp - 0.02425_dp )then
            q = sqrt(-2._dp * log(1._dp - x))
            norm_ppf = -real((((((C(1)*q+C(2))*q+C(3))*q+C(4))*q+C(5))*q+C(6)) / ((((D(1)*q+D(2))*q+D(3))*q+D(4))*q+1._dp))
        else
            q = x - 0.5_dp
            r = q * q
            norm_ppf = real((((((A(1)*r+A(2))*r+A(3))*r+A(4))*r+A(5))*r+A(6))*q / (((((B(1)*r+B(2))*r+B(3))*r+B(4))*r+B(5))*r+1._dp))
        endif
    end function norm_ppf

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

    real function meanm( v, m )
        real,    intent(in) :: v(:)
        logical, intent(in) :: m(:)
        meanm = 0.
        if( count(m) > 0 ) meanm = sum(v, mask=m) / real(count(m))
    end function meanm

    real function varp( v )
        real, intent(in) :: v(:)
        varp = sum((v - sum(v) / real(size(v)))**2) / real(size(v))
    end function varp

    real function pcorr( x, y )
        real, intent(in) :: x(:), y(:)
        real(dp) :: mx, my
        mx = sum(real(x,dp)) / size(x)
        my = sum(real(y,dp)) / size(y)
        pcorr = real(sum((x - mx) * (y - my)) / max(sqrt(sum((x - mx)**2) * sum((y - my)**2)), 1.e-30_dp))
    end function pcorr

    real function cvof( d )
        real, intent(in) :: d(:)
        cvof = 0.
        if( size(d) > 0 ) cvof = sqrt(varp(d)) / max(abs(sum(d) / real(size(d))), 1.e-12)
    end function cvof

    subroutine sort_by( ord, key )
        use simple_srch_sort_loc, only: hpsort
        integer,  intent(inout) :: ord(:)
        real(dp), intent(in)    :: key(:)
        real, allocatable :: k(:)
        k = real(key)
        call hpsort(k, ord)
    end subroutine sort_by

    subroutine ordered_cells( cell, tc, ord )
        integer,              intent(in)  :: cell(:)
        real,                 intent(in)  :: tc(:)
        integer, allocatable, intent(out) :: ord(:)
        real(dp), allocatable :: cv(:)
        integer :: k
        ord = pack([(k, k=1,maxval(cell))], [(any(cell == k), k=1,maxval(cell))])
        allocate(cv(size(ord)))
        do k = 1, size(ord)
            cv(k) = real(meanm(tc, cell == ord(k)), dp)
        end do
        call sort_by(ord, cv)
    end subroutine ordered_cells

    subroutine write_ints( fname, v )
        class(string), intent(in) :: fname
        integer,       intent(in) :: v(:)
        integer :: funit, ios, i
        call fopen(funit, fname, status='replace', action='write', iostat=ios)
        call fileiochk('landscape_sample: cannot write '//fname%to_char(), ios)
        do i = 1, size(v)
            write(funit,'(I0)') v(i)
        end do
        call fclose(funit)
    end subroutine write_ints

    function pad2( k ) result( str )
        integer, intent(in) :: k
        character(len=2) :: str
        write(str,'(I2.2)') k
    end function pad2

    function joined( vals ) result( str )
        integer, intent(in) :: vals(:)
        character(len=:), allocatable :: str
        integer :: i
        str = ''
        do i = 1, size(vals)
            if( i > 1 ) str = str//' '
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
            str = str//fmtr(real(vals(i),dp),3)
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

end module simple_landscape_sample
