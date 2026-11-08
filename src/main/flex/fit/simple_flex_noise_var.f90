!@descr: flex_pca per-particle noise level: an empirical-Bayes noise variance per particle and band from the residual against the node mean
module simple_flex_noise_var
use, intrinsic :: iso_c_binding, only: c_float_complex
use simple_core_module_api,               only: dp, fplane_type, hpsort, logfhandle, maximgbatchsz, median, ori, simple_exception
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

public :: noise_var_on, estimate_noise_var, retire_noise_var, adopt_noise_var, accumulate_moments

integer,          parameter :: NBANDS = 6
integer,          parameter :: BAND_RING(NBANDS) = [2, 6, 8, 10, 12, 14]
integer,          parameter :: KB_REF = 15
character(len=*), parameter :: NOISE_VAR_FILE = 'flex_pca_noise_var.txt'

contains

    logical function noise_var_on()
        character(len=8) :: ev
        integer :: ios
        call get_environment_variable('SIMPLE_FLEX_NOISE_VAR', ev, status=ios)
        noise_var_on = .not. (ios == 0 .and. ev == 'no')
    end function noise_var_on

    subroutine estimate_noise_var( params, build, plane_store, pinds, nptcls )
        class(parameters),       intent(inout) :: params
        class(builder),          intent(inout) :: build
        class(flex_plane_store), intent(inout) :: plane_store
        integer,                 intent(in)    :: pinds(:), nptcls
        type(ptcl_read_session)        :: reads
        type(fplane_type), allocatable :: fpls(:)
        type(ori),         allocatable :: os(:)
        type(image)                    :: vol
        complex(kind=c_float_complex), pointer :: cp(:,:,:)
        complex(dp), allocatable :: mu(:)
        real(dp),    allocatable :: f2(:,:), aa(:,:), bb(:,:), cw(:), beta(:), sig2(:), amp(:), s(:,:), v(:,:)
        real,        allocatable :: band_var(:,:), rmats(:,:,:)
        integer,     allocatable :: hs(:), ks(:), rs(:), band_of(:)
        real(dp)    :: nu(NBANDS), alpha_b(NBANDS), beta_b(NBANDS), res
        integer     :: box, nh, nyq, kb, ns, h, k, r2, i, j, ib, bl(2), bsz, it, b, kft(2), x, z, lg(3), funit, edge(NBANDS)
        logical     :: l_pcache
        box = params%box_crop
        nh  = box / 2 + 1
        nyq = box / 2
        kft = projected_model_kfromto(params%box_crop, params%smpd_crop, params%lp)
        kb  = min(kft(2), nyq - 1)
        edge = noise_band_edges(kb)
        call vol%new([box,box,box], params%smpd_crop)
        call vol%read_and_crop(params%vols(1), params%smpd, box, params%smpd_crop)
        call vol%fft()
        do z = 1, box
            do h = 1, box
                do x = 1, nh
                    lg = vol%comp_addr_logi(x, h, z)
                    if( cell_of(lg, nh, box) /= x + nh * ((h - 1) + box * (z - 1)) ) &
                        &THROW_HARD('noise variance: cell layout differs from the Fourier layout of the image')
                end do
            end do
        end do
        call vol%get_cmat_ptr(cp)
        mu = cmplx(reshape(cp, [nh*box*box]), kind=dp)
        call vol%kill
        allocate(hs((2*kb+1)**2), ks((2*kb+1)**2), rs((2*kb+1)**2))
        ns = 0
        do k = -kb, 0
            do h = -kb, kb
                if( k == 0 .and. h <= 0 ) cycle
                r2 = h * h + k * k
                if( r2 < 4 .or. r2 > kb * (kb + 1) ) cycle
                ns = ns + 1
                hs(ns) = h; ks(ns) = k; rs(ns) = nint(sqrt(real(r2)))
            end do
        end do
        allocate(f2(0:kb,nptcls), aa(0:kb,nptcls), bb(0:kb,nptcls), source=0.d0)
        allocate(os(MAXIMGBATCHSZ), rmats(3,3,MAXIMGBATCHSZ))
        l_pcache = plane_store%cache_in_use()
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
                rmats(:,:,i) = os(i)%get_mat()
            end do
            call accumulate_moments(fpls(:bsz), rmats(:,:,:bsz), hs(:ns), ks(:ns), rs(:ns), mu, nh, box, &
                &f2(:,bl(1):bl(2)), aa(:,bl(1):bl(2)), bb(:,bl(1):bl(2)))
        end do
        call reads%kill
        call cleanup_rec_buffers(build, fpls)
        allocate(cw(0:kb), source=0.d0)
        do j = 1, ns
            cw(rs(j)) = cw(rs(j)) + 2.d0
        end do
        allocate(beta(0:kb), sig2(0:kb), source=1.d0)
        allocate(amp(nptcls), source=1.d0)
        do it = 1, 3
            do h = 0, kb
                beta(h) = max(sum(amp * aa(h,:)) / max(sum(amp**2 * bb(h,:)), 1.d-30), 1.d-12)
                sig2(h) = max(sum(f2(h,:) - 2.d0 * amp * beta(h) * aa(h,:) + (amp * beta(h))**2 * bb(h,:)) / &
                    &max(real(nptcls,dp) * cw(h), 1.d0), 1.d-30)
            end do
            do i = 1, nptcls
                amp(i) = sum(beta(2:) * aa(2:,i) / sig2(2:)) / max(sum(beta(2:)**2 * bb(2:,i) / sig2(2:)), 1.d-30)
            end do
            amp = min(max(amp, 0.05d0), 20.d0)
            amp = amp / median(amp)
        end do
        allocate(band_of(0:2*box))
        do h = 0, 2 * box
            band_of(h) = max(1, count(edge <= h))
        end do
        allocate(s(NBANDS,nptcls), source=0.d0)
        nu = 0.d0
        do h = 2, kb
            b = band_of(h)
            nu(b) = nu(b) + cw(h)
            do i = 1, nptcls
                res = f2(h,i) - 2.d0 * amp(i) * beta(h) * aa(h,i) + (amp(i) * beta(h))**2 * bb(h,i)
                s(b,i) = s(b,i) + res / sig2(h)
            end do
        end do
        allocate(v(NBANDS,nptcls), source=1.d0)
        do b = 1, NBANDS
            if( nu(b) <= 0.d0 ) cycle
            call invgamma_ml(s(b,:), nu(b), alpha_b(b), beta_b(b))
            v(b,:) = (beta_b(b) + 0.5d0 * s(b,:)) / (alpha_b(b) + 0.5d0 * nu(b) - 1.d0)
            v(b,:) = v(b,:) / (sum(v(b,:)) / real(nptcls,dp))
            write(logfhandle,'(A,I0,A,I0,A,I0,A,ES11.4,A,F8.4,A,F7.4,A,3F8.3)') '>>> FLEX_PCA NOISE VARIANCE band ', b, ' (rings ', &
                &edge(b), '-', merge(kb, edge(min(b+1,NBANDS))-1, b == NBANDS), '): inverse-gamma alpha ', alpha_b(b), &
                &' prior mean ', beta_b(b) / max(alpha_b(b) - 1.d0, 1.d-30), ' prior CV ', 1.d0 / sqrt(max(alpha_b(b) - 2.d0, 1.d-30)), &
                &'; multiplier 5/50/95 % ', quant(v(b,:), 0.05d0), quant(v(b,:), 0.5d0), quant(v(b,:), 0.95d0)
        end do
        call flush(logfhandle)
        allocate(band_var(NBANDS, build%spproj_field%get_noris()), source=1.)
        do i = 1, nptcls
            band_var(:,pinds(i)) = real(v(:,i))
        end do
        call plane_store%set_noise_bands(band_var, band_of)
        open(newunit=funit, file=NOISE_VAR_FILE, status='replace', action='write')
        write(funit,'(A,I0,A,I0)') '# bands ', NBANDS, ' kb ', kb
        do b = 1, NBANDS
            write(funit,'(A,I0,A,I0,2ES16.7,F12.4)') '# band ', b, ' ring ', edge(b), alpha_b(b), beta_b(b), nu(b)
        end do
        do i = 1, nptcls
            write(funit,'(I10,I3,6ES14.6)') pinds(i), build%spproj_field%get_eo(pinds(i)), v(:,i)
        end do
        close(funit)

      contains

        real function quant( vv, q )
            real(dp), intent(in) :: vv(:), q
            real, allocatable :: w(:)
            integer :: m
            w = real(vv)
            call hpsort(w)
            m = max(1, min(size(w), nint(q * real(size(w) - 1)) + 1))
            quant = w(m)
        end function quant

    end subroutine estimate_noise_var

    subroutine accumulate_moments( fpls, rmats, hs, ks, rs, mu, nh, box, f2, aa, bb )
        type(fplane_type), intent(in)    :: fpls(:)
        real,              intent(in)    :: rmats(:,:,:)
        integer,           intent(in)    :: hs(:), ks(:), rs(:), nh, box
        complex(dp),       intent(in)    :: mu(:)
        real(dp),          intent(inout) :: f2(0:,:), aa(0:,:), bb(0:,:)
        complex(dp) :: mus, y
        real        :: loc(3), t
        integer     :: i, j, h, k, p(3), nyq
        nyq = box / 2
        !$omp parallel do default(shared) private(i,j,h,k,loc,p,mus,y,t) schedule(static) proc_bind(close)
        do i = 1, size(fpls)
            do j = 1, size(hs)
                h = hs(j)
                k = ks(j)
                y = cmplx(fpls(i)%cmplx_plane(h,k), kind=dp)
                f2(rs(j),i) = f2(rs(j),i) + 2.d0 * abs(y)**2
                loc = real(h) * rmats(1,:,i) + real(k) * rmats(2,:,i)
                p   = nint(loc)
                if( any(abs(p) >= nyq) ) cycle
                if( p(1) < 0 )then
                    mus = conjg(mu(cell_of(-p, nh, box)))
                else
                    mus = mu(cell_of(p, nh, box))
                endif
                t = real(fpls(i)%transfer_plane(h,k))
                aa(rs(j),i) = aa(rs(j),i) + 2.d0 * real(conjg(mus) * real(t,dp) * y, dp)
                bb(rs(j),i) = bb(rs(j),i) + 2.d0 * abs(mus)**2 * real(t,dp)**2
            end do
        end do
        !$omp end parallel do
    end subroutine accumulate_moments

    pure integer function cell_of( lgv, nh, box )
        integer, intent(in) :: lgv(3), nh, box
        cell_of = lgv(1) + 1 + nh * (modulo(lgv(2), box) + box * modulo(lgv(3), box))
    end function cell_of

    pure function noise_band_edges( kb ) result( edge )
        integer, intent(in) :: kb
        integer :: edge(NBANDS), b
        edge(1) = BAND_RING(1)
        do b = 2, NBANDS
            edge(b) = max(edge(b-1) + 1, nint(real(BAND_RING(b) * kb) / real(KB_REF)))
        end do
    end function noise_band_edges

    subroutine invgamma_ml( s, nu, alpha, beta )
        real(dp), intent(in)  :: s(:), nu
        real(dp), intent(out) :: alpha, beta
        real(dp) :: lo, hi, a1, a2, p1, p2, gr
        integer  :: it
        gr = 0.5d0 * (sqrt(5.d0) - 1.d0)
        lo = log(1.0001d0)
        hi = log(1.d9)
        a1 = hi - gr * (hi - lo); p1 = profile(exp(a1))
        a2 = lo + gr * (hi - lo); p2 = profile(exp(a2))
        do it = 1, 200
            if( hi - lo < 1.d-8 ) exit
            if( p1 > p2 )then
                hi = a2; a2 = a1; p2 = p1
                a1 = hi - gr * (hi - lo); p1 = profile(exp(a1))
            else
                lo = a1; a1 = a2; p1 = p2
                a2 = lo + gr * (hi - lo); p2 = profile(exp(a2))
            endif
        end do
        alpha = exp(0.5d0 * (lo + hi))
        beta  = beta_of(alpha)

      contains

        real(dp) function beta_of( a )
            real(dp), intent(in) :: a
            real(dp) :: bl, bh, bm, target, m
            integer  :: jt
            m      = max(sum(s) / real(size(s), dp), 1.d-30)
            target = real(size(s), dp) * a / (a + 0.5d0 * nu)
            bl = log(1.d-12 * m)
            bh = log(1.d12 * m)
            do jt = 1, 200
                bm = 0.5d0 * (bl + bh)
                if( sum(exp(bm) / (exp(bm) + 0.5d0 * s)) < target )then
                    bl = bm
                else
                    bh = bm
                endif
                if( bh - bl < 1.d-12 ) exit
            end do
            beta_of = exp(0.5d0 * (bl + bh))
        end function beta_of

        real(dp) function profile( a )
            real(dp), intent(in) :: a
            real(dp) :: b
            b = beta_of(a)
            profile = sum(a * log(b) + log_gamma(a + 0.5d0 * nu) - log_gamma(a) - (a + 0.5d0 * nu) * log(b + 0.5d0 * s))
        end function profile

    end subroutine invgamma_ml

    subroutine retire_noise_var()
        integer :: funit
        logical :: l_exists
        inquire(file=NOISE_VAR_FILE, exist=l_exists)
        if( .not. l_exists ) return
        open(newunit=funit, file=NOISE_VAR_FILE, status='replace', action='write')
        write(funit,'(A)') '# bands 0'
        close(funit)
    end subroutine retire_noise_var

    subroutine adopt_noise_var( params, build, plane_store )
        class(parameters),       intent(in)    :: params
        class(builder),          intent(inout) :: build
        class(flex_plane_store), intent(inout) :: plane_store
        character(len=512) :: line
        character(len=8)   :: tag
        real,    allocatable :: band_var(:,:)
        integer, allocatable :: band_of(:)
        real    :: vb(NBANDS)
        integer :: funit, ios, pind, eo, h, nb, b, ring, edge(NBANDS)
        logical :: l_exists
        inquire(file=NOISE_VAR_FILE, exist=l_exists)
        if( .not. l_exists ) return
        allocate(band_var(NBANDS, build%spproj_field%get_noris()), source=1.)
        edge = 0
        nb   = -1
        open(newunit=funit, file=NOISE_VAR_FILE, status='old', action='read')
        do
            read(funit,'(A)', iostat=ios) line
            if( ios /= 0 ) exit
            if( line(1:8) == '# bands ' )then
                read(line(9:), *) nb
                if( nb == 0 ) exit
                if( nb /= NBANDS ) THROW_HARD('noise variance file '//NOISE_VAR_FILE//' holds another band count')
                cycle
            endif
            if( line(1:7) == '# band ' )then
                read(line(8:), *, iostat=ios) b, tag, ring
                if( ios /= 0 .or. b < 1 .or. b > NBANDS ) THROW_HARD('noise variance file '//NOISE_VAR_FILE//' has an unreadable band')
                edge(b) = ring
                cycle
            endif
            if( line(1:1) == '#' ) cycle
            read(line, *, iostat=ios) pind, eo, vb
            if( ios /= 0 ) THROW_HARD('noise variance file '//NOISE_VAR_FILE//' has an unreadable row')
            band_var(:,pind) = vb
        end do
        close(funit)
        if( nb == 0 ) return
        if( nb /= NBANDS .or. any(edge < 2) ) THROW_HARD('noise variance file '//NOISE_VAR_FILE//' lacks its band header')
        allocate(band_of(0:2*params%box_crop))
        do h = 0, ubound(band_of,1)
            band_of(h) = max(1, count(edge <= h))
        end do
        call plane_store%set_noise_bands(band_var, band_of)
        write(logfhandle,'(A,A)') '>>> FLEX_PCA NOISE VARIANCE adopted from ', NOISE_VAR_FILE
        call flush(logfhandle)
    end subroutine adopt_noise_var

end module simple_flex_noise_var
