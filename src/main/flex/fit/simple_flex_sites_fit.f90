!@descr: flex_pca sparse sites: the node's observation operator on the coarse lattice and the EM of block-sparse site modes
module simple_flex_sites_fit
!$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num
use, intrinsic :: iso_c_binding, only: c_float_complex
use simple_core_module_api,               only: dp, fplane_type, hpsort, logfhandle, maximgbatchsz, median, ori, simple_exception
use simple_rnd,                           only: ran3
use simple_builder,                       only: builder
use simple_parameters,                    only: parameters
use simple_image,                         only: image
use simple_flex_pca_planes,               only: flex_plane_store
use simple_flex_reconstructor_latent_ops, only: planes_batch_load, projected_model_kfromto
use simple_matcher_3Drec,                 only: init_rec, cleanup_rec_buffers
use simple_matcher_ptcl_io,               only: prepimgbatch, ptcl_read_session
implicit none
private
#include "simple_local_flags.inc"

public :: flex_sites_ops, flex_sites_grid, flex_sites_state, sites_em_iter, sites_gains, sites_posterior_all, group_norm2, sites_diag_precision
public :: sites_fit_nuisance, modes_ft, modes_real

real(dp), parameter :: PI_DP = 3.14159265358979323846d0

type :: flex_sites_ops
    integer  :: box = 0, nh = 0, nc = 0, n3 = 0, nyq = 0, kb = 0, n = 0, nthr = 1
    real     :: smpd = 0., msk = 0.
    integer,     allocatable :: cell_shell(:), cell_lg(:,:), pind(:), eo(:), col_p(:), pid_c(:)
    integer(8),  allocatable :: ip_p(:), ip_c(:)
    real(dp),    allocatable :: cell_fw(:), beta(:), sig2(:), amp(:), mu_real(:)
    real,        allocatable :: t_p(:), t_c(:)
    complex,     allocatable :: r_p(:), r_c(:)
    complex(dp), allocatable :: mu(:)
    logical,     allocatable :: hasdata(:)
    type(image), allocatable :: fimg(:)
  contains
    procedure :: new     => ops_new
    procedure :: to_ft   => ops_to_ft
    procedure :: to_real => ops_to_real
    procedure :: cell_of => ops_cell_of
    procedure :: kill    => ops_kill
end type flex_sites_ops

type :: flex_sites_grid
    integer  :: nsupp = 0, ng = 0, admm_its = 100, rot_sweeps = 8
    real(dp) :: kappa = 3.d0
    integer, allocatable :: sv(:), gid(:), gvox(:)
end type flex_sites_grid

type :: flex_sites_state
    integer :: qs = 0, nb = 0
    real(dp),    allocatable :: x(:,:), z(:,:), u(:,:), zr(:,:), ur(:,:), gn(:,:), alpha(:)
    integer,     allocatable :: blk(:)
    complex(dp), allocatable :: wn(:,:)
    real(dp) :: rho = -1.d0, rhor = -1.d0
end type flex_sites_state

contains


    subroutine ops_new( self, params, build, plane_store, pinds, nptcls )
        class(flex_sites_ops),   intent(inout) :: self
        class(parameters),       intent(in)    :: params
        class(builder),          intent(inout) :: build
        class(flex_plane_store), intent(inout) :: plane_store
        integer,                 intent(in)    :: pinds(:), nptcls
        type(ptcl_read_session)        :: reads
        type(fplane_type), allocatable :: fpls(:)
        type(ori),         allocatable :: os(:)
        type(image)                    :: vol
        integer,     allocatable :: hs(:), ks(:), rs(:), rc(:,:), cnt(:)
        real,        allocatable :: rt(:,:)
        complex,     allocatable :: ry(:,:)
        real(dp),    allocatable :: f2(:,:), aa(:,:), bb(:,:), cw(:), num(:), den(:)
        complex(dp) :: mus, rv
        real(dp)    :: tv, sw, a2b2
        real        :: rmat(3,3), loc(3)
        integer     :: h, k, i, j, c, cc, ib, bl(2), bsz, ns, p(3), r2, x, y, z, it, mate, lg(3), kft(2)
        integer(8)  :: e
        logical     :: l_pcache
        self%box  = params%box_crop
        self%smpd = params%smpd_crop
        self%nh   = self%box / 2 + 1
        self%nc   = self%nh * self%box * self%box
        self%n3   = self%box**3
        self%nyq  = self%box / 2
        kft       = projected_model_kfromto(params%box_crop, params%smpd_crop, params%lp)
        self%kb   = min(kft(2), self%nyq - 1)
        self%msk  = params%mskdiam / 2. / self%smpd
        self%n    = nptcls
        self%nthr = 1
        !$ self%nthr = omp_get_max_threads()
        allocate(self%fimg(self%nthr))
        do i = 1, self%nthr
            call self%fimg(i)%new([self%box,self%box,self%box], self%smpd, wthreads=.false.)
        end do
        allocate(self%cell_shell(self%nc), self%cell_fw(self%nc), self%cell_lg(3,self%nc))
        do z = 1, self%box
            do y = 1, self%box
                do x = 1, self%nh
                    c  = x + self%nh * ((y - 1) + self%box * (z - 1))
                    lg = self%fimg(1)%comp_addr_logi(x, y, z)
                    self%cell_lg(:,c)  = lg
                    self%cell_shell(c) = min(nint(sqrt(real(sum(lg**2)))), self%nyq)
                    self%cell_fw(c)    = merge(2.d0, 1.d0, lg(1) > 0 .and. lg(1) < self%nyq)
                    if( self%cell_of(lg) /= c ) THROW_HARD('sparse sites: cell layout differs from the Fourier layout of the image')
                end do
            end do
        end do
        call vol%new([self%box,self%box,self%box], self%smpd)
        call vol%read_and_crop(params%vols(1), params%smpd, self%box, self%smpd)
        self%mu_real = real(reshape(vol%get_rmat(), [self%n3]), dp)
        call vol%kill
        allocate(self%mu(self%nc))
        call self%to_ft(self%mu_real, self%mu)
        ns = 0
        allocate(hs((2*self%kb+1)**2), ks((2*self%kb+1)**2), rs((2*self%kb+1)**2))
        do k = -self%kb, 0
            do h = -self%kb, self%kb
                if( k == 0 .and. h <= 0 ) cycle
                r2 = h * h + k * k
                if( r2 < 4 .or. r2 > self%kb * (self%kb + 1) ) cycle
                ns = ns + 1
                hs(ns) = h; ks(ns) = k; rs(ns) = nint(sqrt(real(r2)))
            end do
        end do
        allocate(rc(ns,nptcls), rt(ns,nptcls), ry(ns,nptcls), self%pind(nptcls), self%eo(nptcls), os(MAXIMGBATCHSZ))
        self%pind = pinds(1:nptcls)
        l_pcache  = plane_store%cache_in_use()
        call init_rec(params, build, MAXIMGBATCHSZ, fpls, cropped=l_pcache)
        if( l_pcache )then
            call prepimgbatch(params, build, MAXIMGBATCHSZ, box=params%box_crop, smpd=params%smpd_crop)
        else
            call prepimgbatch(params, build, MAXIMGBATCHSZ)
        endif
        call reads%new(params)
        do ib = 1, nptcls, MAXIMGBATCHSZ
            bl  = [ib, min(nptcls, ib + MAXIMGBATCHSZ - 1)]
            bsz = bl(2) - bl(1) + 1
            call planes_batch_load(plane_store, params, build, nptcls, pinds, bl, fpls, l_pcache, reads)
            do i = 1, bsz
                call build%spproj_field%get_ori(pinds(bl(1)+i-1), os(i))
                self%eo(bl(1)+i-1) = build%spproj_field%get_eo(pinds(bl(1)+i-1))
            end do
            !$omp parallel do default(shared) private(i,j,h,k,rmat,loc,p,c) schedule(static) proc_bind(close)
            do i = 1, bsz
                rmat = os(i)%get_mat()
                do j = 1, ns
                    h   = hs(j)
                    k   = ks(j)
                    loc = real(h) * rmat(1,:) + real(k) * rmat(2,:)
                    p   = nint(loc)
                    c   = 0
                    if( all(abs(p) < self%nyq) )then
                        if( p(1) < 0 )then
                            c = -self%cell_of(-p)
                        else
                            c = self%cell_of(p)
                        endif
                    endif
                    rc(j,bl(1)+i-1) = c
                    rt(j,bl(1)+i-1) = real(fpls(i)%transfer_plane(h,k))
                    ry(j,bl(1)+i-1) = fpls(i)%cmplx_plane(h,k)
                end do
            end do
            !$omp end parallel do
        end do
        call reads%kill
        call cleanup_rec_buffers(build, fpls)
        allocate(f2(0:self%kb,nptcls), aa(0:self%kb,nptcls), bb(0:self%kb,nptcls), source=0.d0)
        !$omp parallel do default(shared) private(i,j,c,mus) schedule(static) proc_bind(close)
        do i = 1, nptcls
            do j = 1, ns
                f2(rs(j),i) = f2(rs(j),i) + 2.d0 * abs(cmplx(ry(j,i), kind=dp))**2
                c = rc(j,i)
                if( c == 0 ) cycle
                mus = self%mu(abs(c))
                if( c < 0 ) mus = conjg(mus)
                aa(rs(j),i) = aa(rs(j),i) + 2.d0 * real(conjg(mus) * real(rt(j,i),dp) * cmplx(ry(j,i), kind=dp), dp)
                bb(rs(j),i) = bb(rs(j),i) + 2.d0 * abs(mus)**2 * real(rt(j,i),dp)**2
            end do
        end do
        !$omp end parallel do
        allocate(cw(0:self%kb), source=0.d0)
        do j = 1, ns
            cw(rs(j)) = cw(rs(j)) + 2.d0
        end do
        allocate(self%beta(0:self%kb), self%sig2(0:self%kb), num(0:self%kb), den(0:self%kb))
        allocate(self%amp(nptcls), source=1.d0)
        do it = 1, 3
            do h = 0, self%kb
                self%beta(h) = max(sum(self%amp * aa(h,:)) / max(sum(self%amp**2 * bb(h,:)), 1.d-30), 1.d-12)
                self%sig2(h) = max(sum(f2(h,:) - 2.d0 * self%amp * self%beta(h) * aa(h,:) + (self%amp * self%beta(h))**2 * &
                    &bb(h,:)) / max(real(nptcls,dp) * cw(h), 1.d0), 1.d-30)
            end do
            do i = 1, nptcls
                self%amp(i) = sum(self%beta(2:) * aa(2:,i) / self%sig2(2:)) / &
                    &max(sum(self%beta(2:)**2 * bb(2:,i) / self%sig2(2:)), 1.d-30)
            end do
            self%amp = min(max(self%amp, 0.05d0), 20.d0)
            self%amp = self%amp / median(self%amp)
        end do
        write(logfhandle,'(A,I0,A,F7.3,A,I0,A,F6.1,A,I0,A,I0,A)') '>>> SPARSE_SITES lattice: box ', self%box, ' at ', self%smpd, &
            &' A, band to ring ', self%kb, ' (', real(self%box) * self%smpd / real(self%kb), ' A), ', nptcls, ' particles, ', ns, &
            &' samples per particle'
        allocate(cnt(nptcls), source=0)
        do i = 1, nptcls
            do j = 1, ns
                c = rc(j,i)
                if( c == 0 ) cycle
                cnt(i) = cnt(i) + merge(2, 1, self%cell_lg(1,abs(c)) == 0)
            end do
        end do
        allocate(self%ip_p(nptcls+1))
        self%ip_p(1) = 1
        do i = 1, nptcls
            self%ip_p(i+1) = self%ip_p(i) + cnt(i)
        end do
        allocate(self%col_p(self%ip_p(nptcls+1)-1), self%t_p(self%ip_p(nptcls+1)-1), self%r_p(self%ip_p(nptcls+1)-1))
        !$omp parallel do default(shared) private(i,j,c,cc,e,mus,tv,sw,a2b2,rv,mate) schedule(static) proc_bind(close)
        do i = 1, nptcls
            e = self%ip_p(i)
            do j = 1, ns
                c = rc(j,i)
                if( c == 0 ) cycle
                cc  = abs(c)
                mus = self%mu(cc)
                if( c < 0 ) mus = conjg(mus)
                tv   = real(rt(j,i), dp)
                sw   = 2.d0 / self%sig2(rs(j))
                a2b2 = (self%amp(i) * self%beta(rs(j)) * tv)**2 * sw
                rv   = self%amp(i) * self%beta(rs(j)) * tv * sw * cmplx(ry(j,i), kind=dp) - a2b2 * mus
                if( c < 0 ) rv = conjg(rv)
                if( self%cell_lg(1,cc) == 0 )then
                    mate = self%cell_of([0, -self%cell_lg(2,cc), -self%cell_lg(3,cc)])
                    self%col_p(e) = cc;   self%t_p(e) = real(0.5d0 * a2b2); self%r_p(e) = cmplx(0.5d0 * rv)
                    e = e + 1
                    self%col_p(e) = mate; self%t_p(e) = real(0.5d0 * a2b2); self%r_p(e) = cmplx(0.5d0 * conjg(rv))
                else
                    self%col_p(e) = cc;   self%t_p(e) = real(a2b2);         self%r_p(e) = cmplx(rv)
                endif
                e = e + 1
            end do
        end do
        !$omp end parallel do
        deallocate(rc, rt, ry, f2, aa, bb, cnt)
        allocate(cnt(self%nc), source=0)
        do e = 1, self%ip_p(nptcls+1) - 1
            cnt(self%col_p(e)) = cnt(self%col_p(e)) + 1
        end do
        allocate(self%ip_c(self%nc+1))
        self%ip_c(1) = 1
        do c = 1, self%nc
            self%ip_c(c+1) = self%ip_c(c) + cnt(c)
        end do
        self%hasdata = cnt > 0
        allocate(self%pid_c(self%ip_c(self%nc+1)-1), self%t_c(self%ip_c(self%nc+1)-1), self%r_c(self%ip_c(self%nc+1)-1))
        cnt = 0
        do i = 1, nptcls
            do e = self%ip_p(i), self%ip_p(i+1) - 1
                c = self%col_p(e)
                self%pid_c(self%ip_c(c) + cnt(c)) = i
                self%t_c(self%ip_c(c) + cnt(c))   = self%t_p(e)
                self%r_c(self%ip_c(c) + cnt(c))   = self%r_p(e)
                cnt(c) = cnt(c) + 1
            end do
        end do
        write(logfhandle,'(A,I0,A,F8.1,A,I0,A,3F8.2)') '>>> SPARSE_SITES operator: ', self%ip_p(nptcls+1) - 1, ' samples (', &
            &real(self%ip_p(nptcls+1) - 1) / real(nptcls), ' per particle), ', count(self%hasdata), &
            &' cells with data; amplitude 5/50/95 % ', quant(self%amp, 0.05d0), quant(self%amp, 0.5d0), quant(self%amp, 0.95d0)
        call flush(logfhandle)

      contains

        real function quant( v, q )
            real(dp), intent(in) :: v(:), q
            real, allocatable :: s(:)
            integer :: m
            s = real(v)
            call hpsort(s)
            m = max(1, min(size(s), nint(q * real(size(s) - 1)) + 1))
            quant = s(m)
        end function quant

    end subroutine ops_new

    pure integer function ops_cell_of( self, lg ) result( c )
        class(flex_sites_ops), intent(in) :: self
        integer,               intent(in) :: lg(3)
        c = lg(1) + 1 + self%nh * (modulo(lg(2), self%box) + self%box * modulo(lg(3), self%box))
    end function ops_cell_of

    subroutine ops_to_ft( self, x, w )
        class(flex_sites_ops), intent(inout) :: self
        real(dp),              intent(in)    :: x(:)
        complex(dp),           intent(out)   :: w(:)
        complex(kind=c_float_complex), pointer :: cp(:,:,:)
        integer :: ithr
        ithr = 1
        !$ ithr = omp_get_thread_num() + 1
        call self%fimg(ithr)%set_rmat(reshape(real(x), [self%box,self%box,self%box]), .false.)
        call self%fimg(ithr)%fft()
        call self%fimg(ithr)%get_cmat_ptr(cp)
        w = real(self%n3, dp) * cmplx(reshape(cp, [self%nc]), kind=dp)
    end subroutine ops_to_ft

    subroutine ops_to_real( self, w, x )
        class(flex_sites_ops), intent(inout) :: self
        complex(dp),           intent(in)    :: w(:)
        real(dp),              intent(out)   :: x(:)
        complex(kind=c_float_complex), pointer :: cp(:,:,:)
        integer :: ithr
        ithr = 1
        !$ ithr = omp_get_thread_num() + 1
        call self%fimg(ithr)%get_cmat_ptr(cp)
        cp = reshape(cmplx(w / real(self%n3, dp), kind=c_float_complex), shape(cp))
        call self%fimg(ithr)%set_ft(.true.)
        call self%fimg(ithr)%ifft()
        x = real(reshape(self%fimg(ithr)%get_rmat(), [self%n3]), dp)
    end subroutine ops_to_real

    subroutine ops_kill( self )
        class(flex_sites_ops), intent(inout) :: self
        integer :: i
        if( allocated(self%fimg) )then
            do i = 1, size(self%fimg)
                call self%fimg(i)%kill
            end do
            deallocate(self%fimg)
        endif
        if( allocated(self%cell_shell) ) deallocate(self%cell_shell)
        if( allocated(self%cell_lg) )    deallocate(self%cell_lg)
        if( allocated(self%pind) )       deallocate(self%pind)
        if( allocated(self%eo) )         deallocate(self%eo)
        if( allocated(self%col_p) )      deallocate(self%col_p)
        if( allocated(self%pid_c) )      deallocate(self%pid_c)
        if( allocated(self%ip_p) )       deallocate(self%ip_p)
        if( allocated(self%ip_c) )       deallocate(self%ip_c)
        if( allocated(self%cell_fw) )    deallocate(self%cell_fw)
        if( allocated(self%beta) )       deallocate(self%beta)
        if( allocated(self%sig2) )       deallocate(self%sig2)
        if( allocated(self%amp) )        deallocate(self%amp)
        if( allocated(self%mu_real) )    deallocate(self%mu_real)
        if( allocated(self%t_p) )        deallocate(self%t_p)
        if( allocated(self%t_c) )        deallocate(self%t_c)
        if( allocated(self%r_p) )        deallocate(self%r_p)
        if( allocated(self%r_c) )        deallocate(self%r_c)
        if( allocated(self%mu) )         deallocate(self%mu)
        if( allocated(self%hasdata) )    deallocate(self%hasdata)
        self%n = 0
    end subroutine ops_kill

    subroutine modes_ft( ops, x, w )
        type(flex_sites_ops), intent(inout) :: ops
        real(dp),             intent(in)    :: x(:,:)
        complex(dp),          intent(out)   :: w(:,:)
        integer :: j
        !$omp parallel do default(shared) private(j) schedule(dynamic) proc_bind(close)
        do j = 1, size(x,2)
            call ops%to_ft(x(:,j), w(:,j))
        end do
        !$omp end parallel do
    end subroutine modes_ft

    subroutine modes_real( ops, w, x )
        type(flex_sites_ops), intent(inout) :: ops
        complex(dp),          intent(in)    :: w(:,:)
        real(dp),             intent(out)   :: x(:,:)
        integer :: j
        !$omp parallel do default(shared) private(j) schedule(dynamic) proc_bind(close)
        do j = 1, size(w,2)
            call ops%to_real(w(:,j), x(:,j))
        end do
        !$omp end parallel do
    end subroutine modes_real

    subroutine split_modes( w, wr, wi )
        complex(dp),           intent(in)  :: w(:,:)
        real(dp), allocatable, intent(out) :: wr(:,:), wi(:,:)
        wr = transpose(real(w, dp))
        wi = transpose(aimag(w))
    end subroutine split_modes

    subroutine pstats( ops, i, q, wr, wi, b, g, buf )
        type(flex_sites_ops), intent(in)    :: ops
        integer,              intent(in)    :: i, q
        real(dp),             intent(in)    :: wr(:,:), wi(:,:)
        real(dp),             intent(out)   :: b(q), g(q,q)
        real(dp),             intent(inout) :: buf(:,:)
        integer(8) :: p
        integer    :: c, m, j
        real(dp)   :: s
        b = 0.d0
        m = 0
        do p = ops%ip_p(i), ops%ip_p(i+1) - 1
            c = ops%col_p(p)
            s = sqrt(real(ops%t_p(p), dp))
            b = b + real(ops%r_p(p), dp) * wr(:,c) + real(aimag(ops%r_p(p)), dp) * wi(:,c)
            buf(:,m+1) = s * wr(:,c)
            buf(:,m+2) = s * wi(:,c)
            m = m + 2
        end do
        call dsyrk('U', 'N', q, m, 1.d0, buf, q, 0.d0, g, q)
        do j = 1, q
            g(j+1:q,j) = g(j,j+1:q)
        end do
    end subroutine pstats

    subroutine chol_ll( q, h, b, m, ld, ok )
        integer,  intent(in)    :: q
        real(dp), intent(inout) :: h(q,q)
        real(dp), intent(in)    :: b(q)
        real(dp), intent(out)   :: m(q), ld
        logical,  intent(out)   :: ok
        integer :: info, j
        call dpotrf('U', q, h, q, info)
        ok = info == 0
        m  = 0.d0
        ld = 0.d0
        if( .not. ok ) return
        do j = 1, q
            ld = ld + 2.d0 * log(h(j,j))
        end do
        m = b
        call dpotrs('U', q, 1, h, q, m, q, info)
    end subroutine chol_ll

    subroutine prior_precision( st, qs, qn, pp, ldp )
        type(flex_sites_state), intent(in)  :: st
        integer,                intent(in)  :: qs, qn
        real(dp),               intent(out) :: pp(qs+qn,qs+qn), ldp
        real(dp) :: gi(qn,qn)
        integer  :: j, info
        pp  = 0.d0
        ldp = 0.d0
        do j = 1, qs
            pp(j,j) = 1.d0
        end do
        if( qn == 0 ) return
        gi = st%gn
        call dpotrf('U', qn, gi, qn, info)
        if( info /= 0 ) THROW_HARD('sparse sites: the nuisance prior is not positive definite')
        do j = 1, qn
            ldp = ldp - 2.d0 * log(gi(j,j))
        end do
        call dpotri('U', qn, gi, qn, info)
        do j = 1, qn
            gi(j+1:qn,j) = gi(j,j+1:qn)
        end do
        pp(qs+1:,qs+1:) = gi
    end subroutine prior_precision

    integer function pk( j, k )
        integer, intent(in) :: j, k
        pk = j + (k * (k - 1)) / 2
    end function pk

    subroutine sites_fit_nuisance( ops, wn, fm, gn )
        type(flex_sites_ops),  intent(inout) :: ops
        complex(dp),           intent(in)    :: wn(:,:)
        logical,               intent(in)    :: fm(:)
        real(dp), allocatable, intent(inout) :: gn(:,:)
        real(dp), allocatable :: wr(:,:), wi(:,:), bs(:,:), gs(:,:,:), buf(:,:), acc(:,:,:)
        real(dp) :: gi(size(wn,2),size(wn,2)), h(size(wn,2),size(wn,2)), m(size(wn,2)), gnew(size(wn,2),size(wn,2)), ld, ldp
        integer  :: qn, i, it, j, ithr, nthr, info, maxm
        logical  :: ok
        qn = size(wn,2)
        allocate(gn(qn,qn), source=0.d0)
        if( qn == 0 ) return
        call split_modes(wn, wr, wi)
        maxm = 2 * int(maxval(ops%ip_p(2:) - ops%ip_p(:ops%n)))
        allocate(bs(qn,ops%n), gs(qn,qn,ops%n))
        !$omp parallel default(shared) private(i,buf) proc_bind(close)
        allocate(buf(qn,maxm))
        !$omp do schedule(dynamic,64)
        do i = 1, ops%n
            if( fm(i) ) call pstats(ops, i, qn, wr, wi, bs(:,i), gs(:,:,i), buf)
        end do
        !$omp end do
        deallocate(buf)
        !$omp end parallel
        do j = 1, qn
            gn(j,j) = 1.d0 / max(median(pack(gs(j,j,:), fm)), 1.d-300)
        end do
        nthr = 1
        !$ nthr = omp_get_max_threads()
        allocate(acc(qn,qn,nthr))
        do it = 1, 200
            gi = gn
            call dpotrf('U', qn, gi, qn, info)
            if( info /= 0 ) exit
            ldp = 0.d0
            do j = 1, qn
                ldp = ldp - 2.d0 * log(gi(j,j))
            end do
            call dpotri('U', qn, gi, qn, info)
            do j = 1, qn
                gi(j+1:qn,j) = gi(j,j+1:qn)
            end do
            acc = 0.d0
            !$omp parallel do default(shared) private(i,ithr,h,m,ld,ok,j) schedule(dynamic,256) proc_bind(close)
            do i = 1, ops%n
                if( .not. fm(i) ) cycle
                ithr = 1
                !$ ithr = omp_get_thread_num() + 1
                h = gs(:,:,i) + gi
                call chol_ll(qn, h, bs(:,i), m, ld, ok)
                if( .not. ok ) cycle
                call dpotri('U', qn, h, qn, j)
                do j = 1, qn
                    h(j+1:qn,j) = h(j,j+1:qn)
                end do
                acc(:,:,ithr) = acc(:,:,ithr) + h + spread(m,2,qn) * spread(m,1,qn)
            end do
            !$omp end parallel do
            gnew = sum(acc, 3) / real(count(fm), dp)
            gnew = 0.5d0 * (gnew + transpose(gnew))
            ld = maxval(abs(log([(gnew(j,j) / gn(j,j), j=1,qn)])))
            gn = gnew
            if( ld < 1.d-4 ) exit
        end do
        write(logfhandle,'(A,I0,A,*(1X,ES10.3))') '>>> SPARSE_SITES nuisance prior fitted alone in ', it, ' its; prior sd', &
            &(sqrt(gn(j,j)), j=1,qn)
        call flush(logfhandle)
    end subroutine sites_fit_nuisance

    function sites_diag_precision( ops, w ) result( dp_med )
        type(flex_sites_ops), intent(in) :: ops
        complex(dp),          intent(in) :: w(:,:)
        real(dp), allocatable :: dp_med(:), v(:,:)
        integer(8) :: p
        integer    :: i, q
        q = size(w,2)
        allocate(v(q,ops%n), source=0.d0)
        !$omp parallel do default(shared) private(i,p) schedule(dynamic,256) proc_bind(close)
        do i = 1, ops%n
            do p = ops%ip_p(i), ops%ip_p(i+1) - 1
                v(:,i) = v(:,i) + real(ops%t_p(p), dp) * abs(w(ops%col_p(p),:))**2
            end do
        end do
        !$omp end parallel do
        allocate(dp_med(q))
        do i = 1, q
            dp_med(i) = median(v(i,:))
        end do
    end function sites_diag_precision

    function group_norm2( grd, x, blk, nb ) result( out )
        type(flex_sites_grid), intent(in) :: grd
        real(dp),              intent(in) :: x(:,:)
        integer,               intent(in) :: blk(:), nb
        real(dp), allocatable :: out(:,:)
        integer :: b, j, s
        allocate(out(nb,grd%ng), source=0.d0)
        !$omp parallel do default(shared) private(b,j,s) schedule(dynamic) proc_bind(close)
        do b = 1, nb
            do j = 1, size(x,2)
                if( blk(j) /= b ) cycle
                do s = 1, grd%nsupp
                    out(b,grd%gid(s)) = out(b,grd%gid(s)) + x(grd%sv(s),j)**2
                end do
            end do
        end do
        !$omp end parallel do
    end function group_norm2

    subroutine estep_fit( ops, st, fm, qs, qn, wr, wi, slot, ms, ess, ens, llsum, nfm )
        type(flex_sites_ops),   intent(in)    :: ops
        type(flex_sites_state), intent(inout) :: st
        logical,                intent(in)    :: fm(:)
        integer,                intent(in)    :: qs, qn
        real(dp),               intent(in)    :: wr(:,:), wi(:,:)
        integer,  allocatable,  intent(out)   :: slot(:)
        real(dp), allocatable,  intent(out)   :: ms(:,:), ess(:,:), ens(:,:)
        real(dp),               intent(out)   :: llsum
        integer,                intent(out)   :: nfm
        real(dp), allocatable :: buf(:,:), snn(:,:,:)
        real(dp) :: pp(qs+qn,qs+qn), ldp, b(qs+qn), g(qs+qn,qs+qn), h(qs+qn,qs+qn), m(qs+qn), ld, ll
        integer  :: i, j, k, q, nthr, ithr, maxm, info
        logical  :: ok
        q = qs + qn
        call prior_precision(st, qs, qn, pp, ldp)
        allocate(slot(ops%n), source=0)
        nfm = 0
        do i = 1, ops%n
            if( fm(i) )then
                nfm = nfm + 1
                slot(i) = nfm
            endif
        end do
        allocate(ms(qs,nfm), ess((qs*(qs+1))/2,nfm), ens(qs*qn,nfm), source=0.d0)
        nthr = 1
        !$ nthr = omp_get_max_threads()
        allocate(snn(qn,qn,nthr), source=0.d0)
        maxm  = 2 * int(maxval(ops%ip_p(2:) - ops%ip_p(:ops%n)))
        llsum = 0.d0
        !$omp parallel default(shared) private(i,j,k,ithr,b,g,h,m,ld,ll,ok,info,buf) proc_bind(close)
        allocate(buf(q,maxm))
        !$omp do schedule(dynamic,64) reduction(+:llsum)
        do i = 1, ops%n
            if( slot(i) == 0 ) cycle
            ithr = 1
            !$ ithr = omp_get_thread_num() + 1
            call pstats(ops, i, q, wr, wi, b, g, buf)
            h = g + pp
            call chol_ll(q, h, b, m, ld, ok)
            if( .not. ok ) cycle
            ll = 0.5d0 * dot_product(b, m) - 0.5d0 * ld + 0.5d0 * ldp
            llsum = llsum + ll
            call dpotri('U', q, h, q, info)
            do j = 1, q
                h(j+1:q,j) = h(j,j+1:q)
            end do
            ms(:,slot(i)) = m(:qs)
            do k = 1, qs
                do j = 1, k
                    ess(pk(j,k),slot(i)) = h(j,k) + m(j) * m(k)
                end do
            end do
            do k = 1, qn
                ens((k-1)*qs+1:k*qs,slot(i)) = h(qs+k,:qs) + m(qs+k) * m(:qs)
            end do
            if( qn > 0 ) snn(:,:,ithr) = snn(:,:,ithr) + h(qs+1:,qs+1:) + spread(m(qs+1:),2,qn) * spread(m(qs+1:),1,qn)
        end do
        !$omp end do
        deallocate(buf)
        !$omp end parallel
        if( qn > 0 )then
            st%gn = sum(snn, 3) / real(nfm, dp)
            st%gn = 0.5d0 * (st%gn + transpose(st%gn))
        endif
    end subroutine estep_fit

    subroutine cell_terms( ops, wn, slot, ms, ens, wgt, qs, qn, a, rhs, ess )
        type(flex_sites_ops),     intent(in)  :: ops
        complex(dp),              intent(in)  :: wn(:,:)
        integer,                  intent(in)  :: slot(:), qs, qn
        real(dp),                 intent(in)  :: ms(:,:), ens(:,:), wgt(:)
        real(dp), allocatable,    intent(out) :: a(:,:)
        complex(dp), allocatable, intent(out) :: rhs(:,:)
        real(dp), optional,       intent(in)  :: ess(:,:)
        real(dp)   :: br(qs), bi(qs), cm(qs*qn), t, wv
        integer(8) :: p
        integer    :: c, i, s, k
        allocate(rhs(ops%nc,qs))
        if( present(ess) )then
            allocate(a(size(ess,1),ops%nc), source=0.d0)
        else
            allocate(a(0,0))
        endif
        !$omp parallel do default(shared) private(c,p,i,s,k,t,wv,br,bi,cm) schedule(dynamic,128) proc_bind(close)
        do c = 1, ops%nc
            br = 0.d0; bi = 0.d0; cm = 0.d0
            do p = ops%ip_c(c), ops%ip_c(c+1) - 1
                i = ops%pid_c(p)
                s = slot(i)
                if( s == 0 ) cycle
                wv = wgt(i)
                if( wv == 0.d0 ) cycle
                t  = real(ops%t_c(p), dp) * wv
                if( present(ess) ) a(:,c) = a(:,c) + t * ess(:,s)
                if( qn > 0 ) cm = cm + t * ens(:,s)
                br = br + real(ops%r_c(p), dp) * wv * ms(:,s)
                bi = bi + real(aimag(ops%r_c(p)), dp) * wv * ms(:,s)
            end do
            rhs(c,:) = cmplx(br, bi, kind=dp)
            do k = 1, qn
                rhs(c,:) = rhs(c,:) - wn(c,k) * cm((k-1)*qs+1:k*qs)
            end do
        end do
        !$omp end parallel do
    end subroutine cell_terms

    subroutine fill_alpha( al )
        real(dp), intent(inout) :: al(0:)
        integer :: s0, s1, s
        s0 = -1; s1 = -1
        do s = 0, ubound(al,1)
            if( al(s) > 0.d0 )then
                if( s0 < 0 ) s0 = s
                s1 = s
            endif
        end do
        if( s0 < 0 ) THROW_HARD('sparse sites: no shell carries data')
        al(:s0-1) = al(s0)
        al(s1+1:) = 10.d0 * al(s1)
    end subroutine fill_alpha

    subroutine alpha_init( ops, a, qs, alpha )
        type(flex_sites_ops),  intent(in)  :: ops
        real(dp),              intent(in)  :: a(:,:)
        integer,               intent(in)  :: qs
        real(dp), allocatable, intent(out) :: alpha(:)
        real(dp), allocatable :: d(:)
        integer :: c, j, s
        allocate(d(ops%nc), alpha(0:ops%nyq), source=0.d0)
        do c = 1, ops%nc
            d(c) = sum([(a(pk(j,j),c), j=1,qs)]) / real(qs,dp) / ops%cell_fw(c)
        end do
        do s = 0, ops%nyq
            if( any(ops%cell_shell == s .and. ops%hasdata .and. d > 0.d0) ) &
                &alpha(s) = median(pack(d, ops%cell_shell == s .and. ops%hasdata .and. d > 0.d0))
        end do
        call fill_alpha(alpha)
    end subroutine alpha_init

    real(dp) function rho_start( ops, a, qs, prior_c )
        type(flex_sites_ops), intent(in) :: ops
        real(dp),             intent(in) :: a(:,:), prior_c(:)
        integer,              intent(in) :: qs
        real(dp), allocatable :: d(:)
        integer :: c, j
        allocate(d(ops%nc))
        do c = 1, ops%nc
            d(c) = (sum([(a(pk(j,j),c), j=1,qs)]) / real(qs,dp) + prior_c(c)) / ops%cell_fw(c)
        end do
        rho_start = median(pack(d, ops%hasdata .and. d > 0.d0)) * real(ops%n3, dp)
    end function rho_start

    subroutine cell_inverse( a, extra, qs, minv )
        real(dp),              intent(in)  :: a(:,:), extra(:)
        integer,               intent(in)  :: qs
        real(dp), allocatable, intent(inout) :: minv(:,:)
        real(dp) :: f(qs,qs)
        integer  :: c, j, k, info
        if( .not. allocated(minv) ) allocate(minv(size(a,1),size(a,2)))
        !$omp parallel do default(shared) private(c,j,k,f,info) schedule(static) proc_bind(close)
        do c = 1, size(a,2)
            do k = 1, qs
                do j = 1, k
                    f(j,k) = a(pk(j,k),c)
                end do
                f(k,k) = f(k,k) + extra(c)
            end do
            call dpotrf('U', qs, f, qs, info)
            call dpotri('U', qs, f, qs, info)
            do k = 1, qs
                do j = 1, k
                    minv(pk(j,k),c) = f(j,k)
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine cell_inverse

    subroutine admm( ops, grd, a, rhs, prior_c, z, u, rho, blk, nb, maxit, lam, bsupp, its )
        type(flex_sites_ops),  intent(inout) :: ops
        type(flex_sites_grid), intent(in)    :: grd
        real(dp),              intent(in)    :: a(:,:), prior_c(:)
        complex(dp),           intent(in)    :: rhs(:,:)
        real(dp),              intent(inout) :: z(:,:), u(:,:), rho
        integer,               intent(in)    :: blk(:), nb, maxit
        real(dp), optional,    intent(in)    :: lam(:,:)
        logical,  optional,    intent(in)    :: bsupp(:,:)
        integer,               intent(out)   :: its
        real(dp),    parameter   :: TOL = 2.d-3, RELAX = 1.6d0
        integer,     parameter   :: MINIT = 5
        real(dp),    allocatable :: minv(:,:), x(:,:), y(:,:), zold(:,:), nrm(:,:), shr(:,:), rc(:)
        complex(dp), allocatable :: w(:,:)
        real(dp) :: r, s, ep, ed, f, xre(size(z,2)), xim(size(z,2)), yre(size(z,2)), yim(size(z,2))
        integer  :: it, qs, c, j, sv
        logical  :: l_pen
        qs    = size(z,2)
        l_pen = present(lam)
        if( l_pen ) l_pen = any(lam > 0.d0)
        allocate(rc(ops%nc), w(ops%nc,qs), x(ops%n3,qs), y(ops%n3,qs), zold(ops%n3,qs))
        rc = rho * ops%cell_fw / real(ops%n3, dp)
        call cell_inverse(a, prior_c + rc, qs, minv)
        its = 0
        do it = 1, maxit
            its = it
            x = z - u
            call modes_ft(ops, x, w)
            !$omp parallel do default(shared) private(c,xre,xim,yre,yim) schedule(static) proc_bind(close)
            do c = 1, ops%nc
                xre = real(rhs(c,:), dp)  + rc(c) * real(w(c,:), dp)
                xim = aimag(rhs(c,:)) + rc(c) * aimag(w(c,:))
                call dspmv('U', qs, 1.d0, minv(:,c), xre, 1, 0.d0, yre, 1)
                call dspmv('U', qs, 1.d0, minv(:,c), xim, 1, 0.d0, yim, 1)
                w(c,:) = cmplx(yre, yim, kind=dp)
            end do
            !$omp end parallel do
            call modes_real(ops, w, x)
            y    = RELAX * x + (1.d0 - RELAX) * z + u
            zold = z
            z    = 0.d0
            if( l_pen )then
                nrm = sqrt(group_norm2(grd, y, blk, nb))
                shr = max(0.d0, 1.d0 - (lam / rho) / max(nrm, 1.d-300))
                do j = 1, qs
                    do sv = 1, grd%nsupp
                        z(grd%sv(sv),j) = y(grd%sv(sv),j) * shr(blk(j),grd%gid(sv))
                    end do
                end do
            else if( present(bsupp) )then
                do j = 1, qs
                    do sv = 1, grd%nsupp
                        if( bsupp(blk(j),grd%gid(sv)) ) z(grd%sv(sv),j) = y(grd%sv(sv),j)
                    end do
                end do
            else
                do j = 1, qs
                    z(grd%sv,j) = y(grd%sv,j)
                end do
            endif
            u  = y - z
            r  = sqrt(sum((x - z)**2))
            s  = rho * sqrt(sum((z - zold)**2))
            ep = TOL * max(sqrt(sum(x**2)), sqrt(sum(z**2)), 1.d-300)
            ed = TOL * max(rho * sqrt(sum(u**2)), 1.d-300)
            if( it >= MINIT .and. r <= ep .and. s <= ed ) exit
            if( mod(it,10) == 0 )then
                f = 1.d0
                if( r / ep > 10.d0 * s / ed )then
                    f = 2.d0
                else if( s / ed > 10.d0 * r / ep )then
                    f = 0.5d0
                endif
                if( f /= 1.d0 )then
                    rho = rho * f
                    u   = u / f
                    rc  = rho * ops%cell_fw / real(ops%n3, dp)
                    call cell_inverse(a, prior_c + rc, qs, minv)
                endif
            endif
        end do
    end subroutine admm

    subroutine noise_lambda( ops, grd, wn, slot, ms, ens, qs, qn, blk, nb, lam )
        type(flex_sites_ops),  intent(inout) :: ops
        type(flex_sites_grid), intent(in)    :: grd
        complex(dp),           intent(in)    :: wn(:,:)
        integer,               intent(in)    :: slot(:), qs, qn, blk(:), nb
        real(dp),              intent(in)    :: ms(:,:), ens(:,:)
        real(dp), allocatable, intent(out)   :: lam(:,:)
        real(dp),    allocatable :: sg(:), a(:,:), xr(:,:)
        complex(dp), allocatable :: rhs(:,:)
        integer :: f, i, j
        allocate(lam(nb,grd%ng), source=0.d0)
        allocate(sg(ops%n), xr(ops%n3,qs))
        do f = 1, 2
            do i = 1, ops%n
                sg(i) = merge(1.d0, -1.d0, ran3() < 0.5)
            end do
            call cell_terms(ops, wn, slot, ms, ens, sg, qs, qn, a, rhs)
            do j = 1, qs
                rhs(:,j) = real(ops%n3, dp) * rhs(:,j) / ops%cell_fw
            end do
            call modes_real(ops, rhs, xr)
            lam = lam + group_norm2(grd, xr, blk, nb)
        end do
        lam = sqrt(lam / 2.d0)
    end subroutine noise_lambda

    subroutine alpha_update( ops, wft, a, alpha, live, frac )
        type(flex_sites_ops), intent(in)    :: ops
        complex(dp),          intent(in)    :: wft(:,:)
        real(dp),             intent(in)    :: a(:,:), frac(:)
        real(dp),             intent(inout) :: alpha(0:)
        logical,              intent(in)    :: live(:)
        real(dp), allocatable :: e(:), num(:), den(:)
        logical,  allocatable :: hd(:)
        integer,  allocatable :: idx(:)
        real(dp) :: f(count(live),count(live)), ev(count(live)), work(3*count(live)), acell, fm, fs
        integer  :: c, j, k, nl, info, s
        nl = count(live)
        if( nl == 0 ) return
        idx = pack([(j, j=1,size(live))], live)
        fm  = sum(frac(idx)) / real(nl, dp)
        fs  = sum(frac(idx))
        allocate(e(ops%nc))
        !$omp parallel do default(shared) private(c,j,k,f,ev,work,info,acell) schedule(static) proc_bind(close)
        do c = 1, ops%nc
            do k = 1, nl
                do j = 1, k
                    f(j,k) = a(pk(min(idx(j),idx(k)),max(idx(j),idx(k))),c)
                end do
            end do
            call dsyev('N', 'U', nl, f, nl, ev, work, 3*nl, info)
            acell = alpha(ops%cell_shell(c)) * ops%cell_fw(c)
            e(c)  = ops%cell_fw(c) * (sum(abs(wft(c,idx))**2) + 2.d0 * fm * sum(1.d0 / (max(ev, 0.d0) + acell)))
        end do
        !$omp end parallel do
        allocate(num(0:ops%nyq), den(0:ops%nyq), source=0.d0)
        allocate(hd(0:ops%nyq), source=.false.)
        do c = 1, ops%nc
            s = ops%cell_shell(c)
            num(s) = num(s) + ops%cell_fw(c) * fs
            den(s) = den(s) + e(c)
            hd(s)  = hd(s) .or. ops%hasdata(c)
        end do
        do s = 0, ops%nyq
            alpha(s) = 0.d0
            if( hd(s) .and. den(s) > 0.d0 ) alpha(s) = num(s) / max(den(s), 1.d-300)
        end do
        call fill_alpha(alpha)
    end subroutine alpha_update

    subroutine rotate_blocks( grd, xs, lam, blk, nb, sweeps, rm, p0, p1, nsw )
        type(flex_sites_grid), intent(in)  :: grd
        real(dp),              intent(in)  :: xs(:,:), lam(:,:)
        integer,               intent(in)  :: blk(:), nb, sweeps
        real(dp),              intent(out) :: rm(:,:), p0, p1
        integer,               intent(out) :: nsw
        integer,  parameter   :: NTH = 145
        real(dp), allocatable :: q(:,:,:), bn(:,:), pen(:), nj(:), nk(:), qj(:,:), qk(:,:), rj(:), rk(:)
        real(dp) :: th(NTH), cs(NTH), sn(NTH), ct, st
        integer  :: qs, j, k, s, t, bj, bk, i0, nrot, sw, l
        qs = size(xs,2)
        allocate(q(grd%ng,qs,qs), source=0.d0)
        !$omp parallel do default(shared) private(j,k,s) schedule(dynamic) proc_bind(close)
        do k = 1, qs
            do j = 1, k
                do s = 1, grd%nsupp
                    q(grd%gid(s),j,k) = q(grd%gid(s),j,k) + xs(s,j) * xs(s,k)
                end do
                q(:,k,j) = q(:,j,k)
            end do
        end do
        !$omp end parallel do
        rm = 0.d0
        do j = 1, qs
            rm(j,j) = 1.d0
        end do
        do t = 1, NTH
            th(t) = -0.5d0 * PI_DP + PI_DP * real(t - 1, dp) / real(NTH - 1, dp)
        end do
        cs = cos(th); sn = sin(th)
        i0 = (NTH + 1) / 2
        allocate(bn(grd%ng,nb), pen(NTH))
        call block_norms
        p0  = penalty()
        nsw = 0
        do sw = 1, sweeps
            nsw  = sw
            nrot = 0
            do j = 1, qs
                do k = j + 1, qs
                    bj = blk(j); bk = blk(k)
                    if( bj == bk ) cycle
                    nj = bn(:,bj) - q(:,j,j)
                    nk = bn(:,bk) - q(:,k,k)
                    !$omp parallel do default(shared) private(t) schedule(static) proc_bind(close)
                    do t = 1, NTH
                        pen(t) = sum(lam(bj,:) * sqrt(max(nj + cs(t)**2 * q(:,j,j) - 2.d0 * cs(t) * sn(t) * q(:,j,k) + &
                            &sn(t)**2 * q(:,k,k), 0.d0))) + sum(lam(bk,:) * sqrt(max(nk + sn(t)**2 * q(:,j,j) + &
                            &2.d0 * cs(t) * sn(t) * q(:,j,k) + cs(t)**2 * q(:,k,k), 0.d0)))
                    end do
                    !$omp end parallel do
                    t = minloc(pen, 1)
                    if( t == i0 .or. pen(t) >= pen(i0) * (1.d0 - 1.d-7) ) cycle
                    ct = cs(t); st = sn(t)
                    qj = q(:,:,j); qk = q(:,:,k)
                    q(:,:,j) = ct * qj - st * qk
                    q(:,:,k) = st * qj + ct * qk
                    qj = q(:,j,:); qk = q(:,k,:)
                    q(:,j,:) = ct * qj - st * qk
                    q(:,k,:) = st * qj + ct * qk
                    rj = rm(:,j); rk = rm(:,k)
                    rm(:,j) = ct * rj - st * rk
                    rm(:,k) = st * rj + ct * rk
                    call block_norms
                    nrot = nrot + 1
                end do
            end do
            if( nrot == 0 ) exit
        end do
        p1 = penalty()

      contains

        subroutine block_norms
            bn = 0.d0
            do l = 1, qs
                bn(:,blk(l)) = bn(:,blk(l)) + q(:,l,l)
            end do
        end subroutine block_norms

        real(dp) function penalty()
            integer :: b
            penalty = 0.d0
            do b = 1, nb
                penalty = penalty + sum(lam(b,:) * sqrt(max(bn(:,b), 0.d0)))
            end do
        end function penalty

    end subroutine rotate_blocks

    subroutine sites_em_iter( ops, grd, st, fm, kappa, rot, tag, bsupp )
        type(flex_sites_ops),   intent(inout) :: ops
        type(flex_sites_grid),  intent(in)    :: grd
        type(flex_sites_state), intent(inout) :: st
        logical,                intent(in)    :: fm(:)
        real(dp),               intent(in)    :: kappa
        logical,                intent(in)    :: rot
        character(len=*),       intent(in)    :: tag
        logical, optional,      intent(in)    :: bsupp(:,:)
        complex(dp), allocatable :: wft(:,:), rhs(:,:)
        real(dp),    allocatable :: wr(:,:), wi(:,:), ms(:,:), ess(:,:), ens(:,:), a(:,:), prior_c(:), lam(:,:), lamr(:,:)
        real(dp),    allocatable :: xp(:,:), gnr(:,:), en(:), frac(:), wgt(:), rm(:,:), f(:,:), tmp(:,:)
        integer,     allocatable :: slot(:)
        logical,     allocatable :: act(:,:), live(:)
        real(dp) :: llsum, p0, p1, emax
        integer  :: qs, qn, nfm, c, j, k, b, nsw, its
        character(len=32) :: rtxt
        qs = st%qs
        qn = size(st%wn,2)
        allocate(wft(ops%nc,qs+qn))
        call modes_ft(ops, st%x, wft(:,:qs))
        if( qn > 0 ) wft(:,qs+1:) = st%wn
        call split_modes(wft, wr, wi)
        call estep_fit(ops, st, fm, qs, qn, wr, wi, slot, ms, ess, ens, llsum, nfm)
        deallocate(wr, wi)
        allocate(wgt(ops%n), source=1.d0)
        call cell_terms(ops, st%wn, slot, ms, ens, wgt, qs, qn, a, rhs, ess=ess)
        deallocate(ess)
        if( .not. allocated(st%alpha) ) call alpha_init(ops, a, qs, st%alpha)
        prior_c = st%alpha(ops%cell_shell) * ops%cell_fw
        if( st%rho  < 0.d0 ) st%rho  = rho_start(ops, a, qs, prior_c)
        if( st%rhor < 0.d0 ) st%rhor = st%rho
        rtxt = '-'
        if( present(bsupp) )then
            call admm(ops, grd, a, rhs, prior_c, st%zr, st%ur, st%rhor, st%blk, st%nb, grd%admm_its, bsupp=bsupp, its=its)
            st%z = st%zr
            st%u = 0.d0
        else
            call noise_lambda(ops, grd, st%wn, slot, ms, ens, qs, qn, st%blk, st%nb, lam)
            lam = kappa * lam
            call admm(ops, grd, a, rhs, prior_c, st%z, st%u, st%rho, st%blk, st%nb, grd%admm_its, lam=lam, its=its)
            xp = st%z
            if( rot .and. st%nb > 1 )then
                allocate(lamr(st%nb,grd%ng), rm(qs,qs))
                do b = 1, st%nb
                    lamr(b,:) = sum(lam, 1) / real(st%nb, dp)
                end do
                call rotate_blocks(grd, xp(grd%sv,:), lamr, st%blk, st%nb, grd%rot_sweeps, rm, p0, p1, nsw)
                xp   = matmul(xp, rm)
                st%z = matmul(st%z, rm); st%u = matmul(st%u, rm); st%zr = matmul(st%zr, rm); st%ur = matmul(st%ur, rm)
                rhs  = matmul(rhs, rm)
                allocate(f(qs,qs), tmp(qs,qs))
                !$omp parallel do default(shared) private(c,j,k,f,tmp) schedule(static) proc_bind(close)
                do c = 1, ops%nc
                    do k = 1, qs
                        do j = 1, k
                            f(j,k) = a(pk(j,k),c)
                            f(k,j) = f(j,k)
                        end do
                    end do
                    tmp = matmul(transpose(rm), matmul(f, rm))
                    do k = 1, qs
                        do j = 1, k
                            a(pk(j,k),c) = tmp(j,k)
                        end do
                    end do
                end do
                !$omp end parallel do
                write(rtxt,'(ES10.4,A,ES10.4)') p0, ' -> ', p1
            endif
            act = sqrt(group_norm2(grd, xp, st%blk, st%nb)) > 0.d0
            call admm(ops, grd, a, rhs, prior_c, st%zr, st%ur, st%rhor, st%blk, st%nb, grd%admm_its, bsupp=act, its=its)
        endif
        st%x = st%zr
        gnr  = group_norm2(grd, st%x, st%blk, st%nb)
        en   = sum(gnr, 2)
        emax = max(maxval(en), 1.d-300)
        allocate(live(qs), frac(qs))
        do j = 1, qs
            live(j) = en(st%blk(j)) > 1.d-12 * emax
            frac(j) = real(sum(merge(grd%gvox, 0, gnr(st%blk(j),:) > 0.d0)), dp) / real(ops%n3, dp)
        end do
        call modes_ft(ops, st%x, wft(:,:qs))
        call alpha_update(ops, wft(:,:qs), a, st%alpha, live, frac)
        write(logfhandle,'(A,A,F11.5,A,I0,A,*(1X,I0))') '>>> SPARSE_SITES ', tag//': LL/ptcl ', llsum / real(max(nfm,1), dp), &
            &' on ', nfm, ' particles | active groups/block', (count(gnr(b,:) > 0.d0), b=1,st%nb)
        write(logfhandle,'(A,*(1X,ES9.3))') '>>> SPARSE_SITES    energy/block', (en(b), b=1,st%nb)
        write(logfhandle,'(A,A,A,I0)') '>>> SPARSE_SITES    rotation penalty ', trim(rtxt), ', admm its ', its
        call flush(logfhandle)
    end subroutine sites_em_iter

    subroutine sites_gains( ops, st, wc, cb, ncand, held, gain, se, base_ll )
        type(flex_sites_ops),   intent(inout) :: ops
        type(flex_sites_state), intent(in)    :: st
        complex(dp),            intent(in)    :: wc(:,:)
        integer,                intent(in)    :: cb(:), ncand
        logical,                intent(in)    :: held(:)
        real(dp),               intent(out)   :: gain(ncand), se(ncand), base_ll
        complex(dp), allocatable :: wft(:,:)
        real(dp),    allocatable :: wr(:,:), wi(:,:), buf(:,:), sd(:,:), sd2(:,:)
        integer,     allocatable :: base(:), cj(:,:), nr(:)
        real(dp) :: pp0(st%qs+size(st%wn,2),st%qs+size(st%wn,2)), ldp, lls
        integer  :: qs, qn, nsm, q, nb0, i, k, ithr, nthr, maxm, maxr, nh
        qs  = st%qs
        qn  = size(st%wn,2)
        nsm = size(wc,2)
        q   = qs + nsm + qn
        nb0 = qs + qn
        allocate(wft(ops%nc,q))
        if( qs > 0 ) call modes_ft(ops, st%x, wft(:,:qs))
        wft(:,qs+1:qs+nsm) = wc
        if( qn > 0 ) wft(:,qs+nsm+1:) = st%wn
        call split_modes(wft, wr, wi)
        deallocate(wft)
        allocate(base(nb0))
        base(:qs) = [(i, i=1,qs)]
        base(qs+1:) = [(qs+nsm+i, i=1,qn)]
        allocate(nr(ncand), source=0)
        do i = 1, nsm
            nr(cb(i)) = nr(cb(i)) + 1
        end do
        maxr = maxval(nr)
        allocate(cj(maxr,ncand), source=0)
        nr = 0
        do i = 1, nsm
            nr(cb(i)) = nr(cb(i)) + 1
            cj(nr(cb(i)),cb(i)) = qs + i
        end do
        call prior_precision(st, qs, qn, pp0, ldp)
        nthr = 1
        !$ nthr = omp_get_max_threads()
        allocate(sd(ncand,nthr), sd2(ncand,nthr), source=0.d0)
        maxm = 2 * int(maxval(ops%ip_p(2:) - ops%ip_p(:ops%n)))
        lls  = 0.d0
        nh   = count(held)
        !$omp parallel default(shared) private(i,k,ithr,buf) proc_bind(close)
        allocate(buf(q,maxm))
        !$omp do schedule(dynamic,64) reduction(+:lls)
        do i = 1, ops%n
            if( .not. held(i) ) cycle
            ithr = 1
            !$ ithr = omp_get_thread_num() + 1
            call particle_gains(i, buf, lls, sd(:,ithr), sd2(:,ithr))
        end do
        !$omp end do
        deallocate(buf)
        !$omp end parallel
        do k = 1, ncand
            gain(k) = sum(sd(k,:))
            se(k)   = sqrt(max(sum(sd2(k,:)) - gain(k)**2 / real(max(nh,1), dp), 0.d0))
        end do
        base_ll = lls / real(max(nh,1), dp)

      contains

        subroutine particle_gains( ii, bf, llacc, dsum, d2sum )
            integer,  intent(in)    :: ii
            real(dp), intent(inout) :: bf(:,:), llacc, dsum(:), d2sum(:)
            real(dp) :: b(q), g(q,q), h(nb0,nb0), v0(nb0), ld, gx(nb0,maxr), s(maxr,maxr), u(maxr), us(maxr), d
            integer  :: kk, r, info, jj
            logical  :: ok
            call pstats(ops, ii, q, wr, wi, b, g, bf)
            h = g(base,base) + pp0
            call chol_ll(nb0, h, b(base), v0, ld, ok)
            if( .not. ok ) return
            llacc = llacc + 0.5d0 * dot_product(b(base), v0) - 0.5d0 * ld + 0.5d0 * ldp
            do kk = 1, ncand
                r = nr(kk)
                if( r == 0 ) cycle
                gx(:,:r) = g(base,cj(:r,kk))
                call dpotrs('U', nb0, r, h, nb0, gx, nb0, info)
                s(:r,:r) = g(cj(:r,kk),cj(:r,kk)) - matmul(transpose(g(base,cj(:r,kk))), gx(:,:r))
                do jj = 1, r
                    s(jj,jj) = s(jj,jj) + 1.d0
                end do
                u(:r) = b(cj(:r,kk)) - matmul(transpose(gx(:,:r)), b(base))
                call dpotrf('U', r, s, maxr, info)
                if( info /= 0 ) cycle
                us(:r) = u(:r)
                call dpotrs('U', r, 1, s, maxr, us, maxr, info)
                d = 0.5d0 * dot_product(u(:r), us(:r))
                do jj = 1, r
                    d = d - log(s(jj,jj))
                end do
                dsum(kk)  = dsum(kk)  + d
                d2sum(kk) = d2sum(kk) + d * d
            end do
        end subroutine particle_gains

    end subroutine sites_gains

    subroutine sites_posterior_all( ops, st, m, sdv, ll )
        type(flex_sites_ops),   intent(inout) :: ops
        type(flex_sites_state), intent(in)    :: st
        real(dp), allocatable,  intent(out)   :: m(:,:), sdv(:,:), ll(:)
        complex(dp), allocatable :: wft(:,:)
        real(dp),    allocatable :: wr(:,:), wi(:,:), buf(:,:)
        real(dp) :: pp(st%qs+size(st%wn,2),st%qs+size(st%wn,2)), ldp, b(st%qs+size(st%wn,2)), &
            &g(st%qs+size(st%wn,2),st%qs+size(st%wn,2)), h(st%qs+size(st%wn,2),st%qs+size(st%wn,2)), ld
        integer  :: qs, qn, q, i, j, maxm, info
        logical  :: ok
        qs = st%qs
        qn = size(st%wn,2)
        q  = qs + qn
        allocate(wft(ops%nc,q))
        call modes_ft(ops, st%x, wft(:,:qs))
        if( qn > 0 ) wft(:,qs+1:) = st%wn
        call split_modes(wft, wr, wi)
        deallocate(wft)
        call prior_precision(st, qs, qn, pp, ldp)
        allocate(m(q,ops%n), sdv(qs,ops%n), ll(ops%n), source=0.d0)
        maxm = 2 * int(maxval(ops%ip_p(2:) - ops%ip_p(:ops%n)))
        !$omp parallel default(shared) private(i,j,b,g,h,ld,ok,info,buf) proc_bind(close)
        allocate(buf(q,maxm))
        !$omp do schedule(dynamic,64)
        do i = 1, ops%n
            call pstats(ops, i, q, wr, wi, b, g, buf)
            h = g + pp
            call chol_ll(q, h, b, m(:,i), ld, ok)
            if( .not. ok ) cycle
            ll(i) = 0.5d0 * dot_product(b, m(:,i)) - 0.5d0 * ld + 0.5d0 * ldp
            call dpotri('U', q, h, q, info)
            do j = 1, qs
                sdv(j,i) = sqrt(max(h(j,j), 0.d0))
            end do
        end do
        !$omp end do
        deallocate(buf)
        !$omp end parallel
    end subroutine sites_posterior_all

end module simple_flex_sites_fit
