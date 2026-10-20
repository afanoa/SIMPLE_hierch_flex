!@descr: kernel-latent EM (klem) core: band construction, band-limited particle loading at the particle's own pose,
!! Fourier-slice projection of coarse-grid volumes, the cross-half E-step, the nearest-voxel M-step, the lattice kernel
!! and the Sinkhorn balance. Numerically the twin of scripts/cls_expansion_proto/toy55/klem.py (validated E-step):
!!   half-plane band to shell kmax (ky outer, kx inner; kx > 0 or kx == 0, ky >= 0; wq 2/1), F = SIMPLE's centred FT
!!   x box^2 (undoing the 1/N) x contrast sign x exp(+2 pi i (h ss fs x + k ss fs y)/box), C = sin(pi wl k2 (df -
!!   0.5 wl^2 k2 cs) + atan(A/sqrt(1-A^2))) in pixel units (|C| when ctfabs), R = Rz(-e3) Ry(e2) Rz(-e1), loc = (kx,ky,0) M
!!   (M = R or R^T by conv), trilinear interpolation with clamping on the centred coarse grid, ll = -sum wq |F - C SC P|^2
!!   / sig2 against the OTHER half's volumes (eo 0 -> half 2), M-step = nearest-voxel accumulation with Friedel mates.
module simple_klem_core
use simple_core_module_api
use simple_sp_project, only: sp_project
use simple_image,      only: image
use simple_string,     only: string
use simple_stack_io,   only: stack_io
use simple_srch_sort_loc, only: hpsort
use simple_kbinterpol,    only: kbinterpol
!$ use omp_lib
implicit none
public :: klem_band, klem_data, klem_band_new, klem_load_particles, klem_project_cell, klem_estep_ll, klem_mstep_accum
public :: klem_lattice_kernel, klem_sinkhorn, klem_rotmat, klem_write_volume, klem_pose_step, klem_rotmat2euler
public :: klem_set_gridding, klem_gridding_on, klem_deapod_envelope
private
#include "simple_local_flags.inc"

type klem_band
    integer              :: kmax = 0, nband = 0, n3 = 0
    integer, allocatable :: h(:), k(:), shell(:)
    real,    allocatable :: kx(:), ky(:), wq(:)
end type klem_band

type klem_data
    integer :: n = 0, nband = 0
    complex, allocatable :: F(:,:)     ! (nband, n) band coefficients with shift and sign folded in
    real,    allocatable :: C(:,:)     ! (nband, n) CTF (or |CTF|)
    real,    allocatable :: R(:,:,:)   ! (3,3,n) rotation matrices
    integer, allocatable :: eo(:)      ! 0 even, 1 odd
end type klem_data

! Kaiser-Bessel gridding (M-step insertion and E-step slice interpolation) instead of nearest voxel / trilinear; off by default
type(kbinterpol), save :: kbw
logical,          save :: l_grid = .false.
integer,          save :: grid_mode = 0   ! 1 = Kaiser-Bessel 3x3x3, 2 = trilinear 2x2x2

contains

    subroutine klem_set_gridding( mode )
        integer, intent(in) :: mode
        grid_mode = mode; l_grid = mode /= 0
        if( l_grid ) kbw = kbinterpol(KBWINSZ, KBALPHA)
    end subroutine klem_set_gridding

    !> stencil weight at integer offset d from a fractional position t in [0,1): KB (mode 1) or linear (mode 2)
    pure real function gw( x )
        real, intent(in) :: x
        if( grid_mode == 2 )then
            gw = max(0., 1. - abs(x))
        else
            gw = kbw%apod(x)
        endif
    end function gw

    pure logical function klem_gridding_on()
        klem_gridding_on = l_grid
    end function klem_gridding_on

    !> real-space deapodisation envelope of the discrete KB stencil (half-width KBWINSZ, offsets -1..1, period box):
    !! e(x) = sum_d w(d) cos(2 pi d x / box) per axis, normalised, squared (M-step insertion and E-step interpolation),
    !! inverted and capped at gmax
    subroutine klem_deapod_envelope( box, gmax, dm )
        integer, intent(in)  :: box
        real,    intent(in)  :: gmax
        real,    intent(out) :: dm(box,box,box)
        real    :: w(-1:1), e1d(box), x, e
        integer :: d, i, j, k, c
        if( grid_mode == 2 )then
            w = [0.25, 0.5, 0.25]
        else
            do d = -1,1
                w(d) = kbw%apod(real(d))
            enddo
        endif
        w = w / sum(w); c = box / 2 + 1
        do i = 1,box
            x = real(i - c) / real(box); e = 0.
            do d = -1,1
                e = e + w(d) * cos(2. * PI * real(d) * x)
            enddo
            e1d(i) = e
        enddo
        do k = 1,box
            do j = 1,box
                do i = 1,box
                    e = (e1d(i) * e1d(j) * e1d(k))**2
                    dm(i,j,k) = min(gmax, 1. / max(e, 1.e-6))
                enddo
            enddo
        enddo
    end subroutine klem_deapod_envelope

    !> Kaiser-Bessel interpolation on the centred grid (3x3x3 stencil, indices clamped), normalised by the stencil sum
    pure function klem_kbinterp( Vc, kmax, loc ) result( val )
        integer, intent(in) :: kmax
        complex, intent(in) :: Vc(-kmax:kmax,-kmax:kmax,-kmax:kmax)
        real,    intent(in) :: loc(3)
        complex :: val
        real    :: x(3), w, wx, wy, wsum
        integer :: ix, iy, iz, cx, cy, cz, d
        do d = 1,3
            x(d) = min(max(loc(d), -real(kmax)), real(kmax))
        enddo
        cx = nint(x(1)); cy = nint(x(2)); cz = nint(x(3)); val = cmplx(0.,0.); wsum = 0.
        do ix = max(cx-1,-kmax), min(cx+1,kmax)
            wx = gw(x(1) - real(ix))
            do iy = max(cy-1,-kmax), min(cy+1,kmax)
                wy = wx * gw(x(2) - real(iy))
                do iz = max(cz-1,-kmax), min(cz+1,kmax)
                    w = wy * gw(x(3) - real(iz)); val = val + w * Vc(ix,iy,iz); wsum = wsum + w
                enddo
            enddo
        enddo
        if( wsum > 0. ) val = val / wsum
    end function klem_kbinterp

    !> the half-plane band of a box to shell kmax
    subroutine klem_band_new( bd, kmax )
        type(klem_band), intent(inout) :: bd
        integer,         intent(in)    :: kmax
        integer :: h, k, q
        real    :: rr
        bd%kmax = kmax; bd%n3 = 2 * kmax + 1; q = 0
        do k = -kmax,kmax
            do h = -kmax,kmax
                rr = sqrt(real(h*h + k*k))
                if( rr <= real(kmax) .and. (h > 0 .or. (h == 0 .and. k >= 0)) ) q = q + 1
            enddo
        enddo
        bd%nband = q
        if( allocated(bd%h) ) deallocate(bd%h, bd%k, bd%shell, bd%kx, bd%ky, bd%wq)
        allocate(bd%h(q), bd%k(q), bd%shell(q), bd%kx(q), bd%ky(q), bd%wq(q))
        q = 0
        do k = -kmax,kmax
            do h = -kmax,kmax
                rr = sqrt(real(h*h + k*k))
                if( rr <= real(kmax) .and. (h > 0 .or. (h == 0 .and. k >= 0)) )then
                    q = q + 1; bd%h(q) = h; bd%k(q) = k; bd%kx(q) = real(h); bd%ky(q) = real(k)
                    bd%wq(q) = merge(2., 1., h > 0); bd%shell(q) = nint(rr)
                endif
            enddo
        enddo
    end subroutine klem_band_new

    !> R = Rz(-e3) Ry(e2) Rz(-e1) (the prototype's euler2m([-e1, e2, -e3]))
    pure function klem_rotmat( e ) result( R )
        real, intent(in) :: e(3)
        real :: R(3,3), A(3,3), B(3,3), Cm(3,3), ca, sa
        ca = cos(-e(3) * real(PI) / 180.); sa = sin(-e(3) * real(PI) / 180.)
        A = 0.; A(1,1) = ca; A(1,2) = -sa; A(2,1) = sa; A(2,2) = ca; A(3,3) = 1.
        ca = cos(e(2) * real(PI) / 180.);  sa = sin(e(2) * real(PI) / 180.)
        B = 0.; B(1,1) = ca; B(1,3) = -sa; B(3,1) = sa; B(3,3) = ca; B(2,2) = 1.
        ca = cos(-e(1) * real(PI) / 180.); sa = sin(-e(1) * real(PI) / 180.)
        Cm = 0.; Cm(1,1) = ca; Cm(1,2) = -sa; Cm(2,1) = sa; Cm(2,2) = ca; Cm(3,3) = 1.
        R = matmul(matmul(A, B), Cm)
    end function klem_rotmat

    !> band coefficients, CTF and rotation of the listed project rows (serial reads in batches, parallel transforms)
    subroutine klem_load_particles( spproj, bd, pinds, eos, box, smpd, sshift, fscale, csign, ctfabs, dat )
        class(sp_project), intent(inout) :: spproj
        type(klem_band),   intent(in)    :: bd
        integer,           intent(in)    :: pinds(:), eos(:), box
        real,              intent(in)    :: smpd, sshift, fscale, csign
        logical,           intent(in)    :: ctfabs
        type(klem_data),   intent(inout) :: dat
        type(image), allocatable :: imgs(:)
        type(string)    :: stkname
        type(ctfparams) :: ctfp
        type(stack_io)  :: stkio
        real(dp) :: wl, cs_px, dfxp, dfyp, df, k2, ph, ang, phase, afac
        real     :: e(3), sh(2)
        integer  :: n, i, j, q, hh, kk, nb, ib, i0, i1, t, cur_stk
        integer, allocatable :: stk_of(:), ind_of(:), ord(:), key(:)
        complex  :: F_q
        real(dp) :: t_read, t_comp, t_a
        n = size(pinds); dat%n = n; dat%nband = bd%nband
        if( allocated(dat%F) ) deallocate(dat%F, dat%C, dat%R, dat%eo)
        allocate(dat%F(bd%nband,n), source=cmplx(0.,0.)); allocate(dat%C(bd%nband,n), source=0.); allocate(dat%R(3,3,n), source=0.); allocate(dat%eo(n), source=0)
        dat%eo = eos
        nb = 512; allocate(imgs(nb)); t_read = 0.d0; t_comp = 0.d0; t_a = wtime_ld()
        do i = 1,nb
            call imgs(i)%new([box,box,1], smpd)
        enddo
        ! particles sorted into stack order (stack index, index in stack) so that every stack file is opened once and
        ! read sequentially through the buffered stack reader (project-row order is a permutation of stack order:
        ! one open/close + a random seek per particle was the loading cost)
        allocate(stk_of(n), ind_of(n), ord(n), key(n))
        do i = 1,n
            call spproj%map_ptcl_ind2stk_ind('ptcl3D', pinds(i), stk_of(i), ind_of(i))
            key(i) = stk_of(i) * 10000000 + ind_of(i); ord(i) = i
        enddo
        call hpsort(key, ord)
        cur_stk = -1
        do ib = 1,n,nb
            i0 = ib; i1 = min(ib + nb - 1, n)
            do j = i0,i1
                i = ord(j)
                if( stk_of(i) /= cur_stk )then
                    if( cur_stk > 0 ) call stkio%close
                    stkname = spproj%os_stk%get_str(stk_of(i), 'stk')
                    call stkio%open(stkname, smpd, 'read', bufsz=nb)
                    cur_stk = stk_of(i)
                endif
                call stkio%read(ind_of(i), imgs(j-i0+1))
                t_read = t_read + (wtime_ld() - t_a); t_a = wtime_ld()
            enddo
            !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
            !$omp private(i,j,t,q,hh,kk,ctfp,e,sh,F_q,wl,cs_px,dfxp,dfyp,df,k2,ph,ang,phase,afac)
            do j = i0,i1
                i = ord(j); t = j - i0 + 1
                call imgs(t)%fft()
                ctfp = spproj%get_ctfparams('ptcl3D', pinds(i))
                e    = spproj%os_ptcl3D%get_euler(pinds(i))
                sh   = spproj%os_ptcl3D%get_2Dshift(pinds(i))
                dat%R(:,:,i) = klem_rotmat(e)
                wl    = 12.26d0 / sqrt(real(ctfp%kv,dp) * 1.d3 + 0.9784d0 * (real(ctfp%kv,dp) * 1.d3)**2 * 1.d-6) / real(smpd,dp)
                cs_px = real(ctfp%cs,dp) * 1.d7 / real(smpd,dp)
                dfxp  = real(ctfp%dfx,dp) * 1.d4 / real(smpd,dp)
                dfyp  = real(ctfp%dfy,dp) * 1.d4 / real(smpd,dp)
                afac  = atan(real(ctfp%fraca,dp) / sqrt(1.d0 - real(ctfp%fraca,dp)**2))
                do q = 1,bd%nband
                    hh = bd%h(q); kk = bd%k(q)
                    F_q   = imgs(t)%get_fcomp2D(hh,kk) * real(box)**2
                    phase = 2.d0 * PI * (real(hh,dp) * real(sshift*fscale*sh(1),dp) + real(kk,dp) * real(sshift*fscale*sh(2),dp)) / real(box,dp)
                    dat%F(q,i) = F_q * csign * cmplx(real(cos(phase)), real(sin(phase)))
                    ang = atan2(real(kk,dp), real(hh,dp))
                    df  = 0.5d0 * (dfxp + dfyp + cos(2.d0 * (ang - real(ctfp%angast,dp) * PI / 180.d0)) * (dfxp - dfyp))
                    k2  = real(hh*hh + kk*kk, dp) / real(box,dp)**2
                    ph  = PI * wl * k2 * (df - 0.5d0 * wl * wl * k2 * cs_px) + afac
                    dat%C(q,i) = real(sin(ph))
                    if( ctfabs ) dat%C(q,i) = abs(dat%C(q,i))
                enddo
            enddo
            !$omp end parallel do
            t_comp = t_comp + (wtime_ld() - t_a); t_a = wtime_ld()
        enddo
        if( cur_stk > 0 ) call stkio%close
        do i = 1,nb
            call imgs(i)%kill
        enddo
        write(logfhandle,'(A,I0,A,F8.1,A,F8.1,A)') '>>> klem_load_particles: ', n, ' particles; stack reads ', t_read, ' s, FFT + band extraction ', t_comp, ' s'; flush(logfhandle)
        deallocate(imgs, stk_of, ind_of, ord, key)
    end subroutine klem_load_particles

    function wtime_ld() result( t )
        real(dp) :: t
        integer(8) :: cnt_, rate_
        call system_clock(cnt_, rate_); t = real(cnt_,dp) / real(max(rate_,1_8),dp)
    end function wtime_ld

    !> trilinear interpolation on the centred grid with the coordinates clamped (scipy mode='nearest')
    pure function klem_trilinear( Vc, kmax, loc ) result( val )
        integer, intent(in) :: kmax
        complex, intent(in) :: Vc(-kmax:kmax,-kmax:kmax,-kmax:kmax)
        real,    intent(in) :: loc(3)
        complex :: val
        real    :: x(3), f(3)
        integer :: a(3), b(3), d
        do d = 1,3
            x(d) = min(max(loc(d), -real(kmax)), real(kmax))
            a(d) = floor(x(d)); b(d) = min(a(d) + 1, kmax); f(d) = x(d) - real(a(d))
        enddo
        val = (1.-f(1))*(1.-f(2))*(1.-f(3)) * Vc(a(1),a(2),a(3)) + f(1)*(1.-f(2))*(1.-f(3)) * Vc(b(1),a(2),a(3)) &
            + (1.-f(1))*f(2)*(1.-f(3)) * Vc(a(1),b(2),a(3)) + (1.-f(1))*(1.-f(2))*f(3) * Vc(a(1),a(2),b(3)) &
            + f(1)*f(2)*(1.-f(3)) * Vc(b(1),b(2),a(3)) + f(1)*(1.-f(2))*f(3) * Vc(b(1),a(2),b(3)) &
            + (1.-f(1))*f(2)*f(3) * Vc(a(1),b(2),b(3)) + f(1)*f(2)*f(3) * Vc(b(1),b(2),b(3))
    end function klem_trilinear

    !> the band coefficients of a cell volume projected at a rotation (loc = (kx,ky,0) M)
    subroutine klem_project_cell( bd, kc, Vc, R, conv, P )
        type(klem_band), intent(in)  :: bd
        integer,         intent(in)  :: kc
        complex,         intent(in)  :: Vc(-kc:kc,-kc:kc,-kc:kc)
        real,            intent(in)  :: R(3,3)
        integer,         intent(in)  :: conv
        complex,         intent(out) :: P(bd%nband)
        real    :: M(3,3), loc(3)
        integer :: q
        if( conv == 0 )then
            M = R
        else
            M = transpose(R)
        endif
        P = cmplx(0.,0.)
        do q = 1,bd%nband
            if( bd%shell(q) > kc ) cycle
            loc(1) = bd%kx(q) * M(1,1) + bd%ky(q) * M(2,1)
            loc(2) = bd%kx(q) * M(1,2) + bd%ky(q) * M(2,2)
            loc(3) = bd%kx(q) * M(1,3) + bd%ky(q) * M(2,3)
            if( l_grid )then
                P(q) = klem_kbinterp(Vc, kc, loc)
            else
                P(q) = klem_trilinear(Vc, kc, loc)
            endif
        enddo
    end subroutine klem_project_cell

    !> cross-half E-step: ll(k,i) = -sum_q wq |F - C SC P_k^{other half}|^2 / sig2
    !> bd = the band the data was loaded at; kc = the current band limit (coefficients beyond it are skipped); V on the kc grid
    subroutine klem_estep_ll( bd, kc, dat, V, ncl, conv, sig2, scl, ll )
        type(klem_band), intent(in)  :: bd
        integer,         intent(in)  :: kc
        type(klem_data), intent(in)  :: dat
        integer,         intent(in)  :: ncl, conv
        complex,         intent(in)  :: V(-kc:kc,-kc:kc,-kc:kc,ncl,2)
        real,            intent(in)  :: sig2(bd%nband), scl(bd%nband)
        real,            intent(out) :: ll(ncl,dat%n)
        complex, allocatable :: P(:)
        real(dp) :: acc
        integer  :: i, ic, q, ihalf
        !$omp parallel do default(shared) schedule(dynamic,16) proc_bind(close) private(i,ic,q,ihalf,P,acc)
        do i = 1,dat%n
            allocate(P(bd%nband))
            ihalf = merge(2, 1, dat%eo(i) == 0)
            do ic = 1,ncl
                call klem_project_cell(bd, kc, V(:,:,:,ic,ihalf), dat%R(:,:,i), conv, P)
                acc = 0.d0
                do q = 1,bd%nband
                    if( bd%shell(q) > kc ) cycle
                    acc = acc + real(bd%wq(q),dp) * real(abs(dat%F(q,i) - dat%C(q,i) * scl(q) * P(q)),dp)**2 / real(sig2(q),dp)
                enddo
                ll(ic,i) = -real(acc)
            enddo
            deallocate(P)
        enddo
        !$omp end parallel do
    end subroutine klem_estep_ll

    !> nearest-voxel M-step accumulation: num/den per component and half (the prototype's klem_mstep kernel)
    subroutine klem_mstep_accum( bd, kc, dat, W, ncl, conv, num, den )
        type(klem_band), intent(in)    :: bd
        integer,         intent(in)    :: kc
        type(klem_data), intent(in)    :: dat
        integer,         intent(in)    :: ncl, conv
        real,            intent(in)    :: W(ncl,dat%n)
        complex(dp),     intent(inout) :: num(-kc:kc,-kc:kc,-kc:kc,ncl,2)
        real(dp),        intent(inout) :: den(-kc:kc,-kc:kc,-kc:kc,ncl,2)
        real     :: M(3,3), loc(3), wt, gx, gy, g
        real(dp) :: c2
        complex  :: cf
        integer  :: ic, i, q, ih, ix, iy, iz, jx, jy, jz, kmax
        kmax = kc
        !$omp parallel do default(shared) schedule(dynamic,1) proc_bind(close) private(ic,i,q,ih,ix,iy,iz,jx,jy,jz,M,loc,wt,c2,cf,gx,gy,g)
        do ic = 1,ncl
            do i = 1,dat%n
                wt = W(ic,i)
                if( wt < 1.e-4 ) cycle
                ih = merge(1, 2, dat%eo(i) == 0)
                if( conv == 0 )then
                    M = dat%R(:,:,i)
                else
                    M = transpose(dat%R(:,:,i))
                endif
                do q = 1,bd%nband
                    if( bd%shell(q) > kc ) cycle
                    loc(1) = bd%kx(q) * M(1,1) + bd%ky(q) * M(2,1)
                    loc(2) = bd%kx(q) * M(1,2) + bd%ky(q) * M(2,2)
                    loc(3) = bd%kx(q) * M(1,3) + bd%ky(q) * M(2,3)
                    ix = nint(loc(1)); iy = nint(loc(2)); iz = nint(loc(3))
                    cf = wt * dat%C(q,i) * dat%F(q,i); c2 = real(wt * dat%C(q,i) * dat%C(q,i), dp)
                    if( l_grid )then
                        do jx = ix-1,ix+1
                            if( abs(jx) > kmax ) cycle
                            gx = gw(loc(1) - real(jx)); if( gx < 1.e-6 ) cycle
                            do jy = iy-1,iy+1
                                if( abs(jy) > kmax ) cycle
                                gy = gx * gw(loc(2) - real(jy)); if( gy < 1.e-6 ) cycle
                                do jz = iz-1,iz+1
                                    if( abs(jz) > kmax ) cycle
                                    g = gy * gw(loc(3) - real(jz))
                                    if( g < 1.e-6 ) cycle
                                    num(jx,jy,jz,ic,ih) = num(jx,jy,jz,ic,ih) + g * cf
                                    den(jx,jy,jz,ic,ih) = den(jx,jy,jz,ic,ih) + real(g,dp) * c2
                                    num(-jx,-jy,-jz,ic,ih) = num(-jx,-jy,-jz,ic,ih) + g * conjg(cf)
                                    den(-jx,-jy,-jz,ic,ih) = den(-jx,-jy,-jz,ic,ih) + real(g,dp) * c2
                                enddo
                            enddo
                        enddo
                    else if( abs(ix) <= kmax .and. abs(iy) <= kmax .and. abs(iz) <= kmax )then
                        num(ix,iy,iz,ic,ih) = num(ix,iy,iz,ic,ih) + cf
                        den(ix,iy,iz,ic,ih) = den(ix,iy,iz,ic,ih) + c2
                        num(-ix,-iy,-iz,ic,ih) = num(-ix,-iy,-iz,ic,ih) + conjg(cf)
                        den(-ix,-iy,-iz,ic,ih) = den(-ix,-iy,-iz,ic,ih) + c2
                    endif
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine klem_mstep_accum

    !> column-normalised Gaussian kernel on the n1 x n2 lattice (column m = the kernel of target cell m)
    subroutine klem_lattice_kernel( n1, n2, sigma, Ksm, i1, i2 )
        integer, intent(in)  :: n1, n2
        real,    intent(in)  :: sigma
        real,    intent(out) :: Ksm(n1*n2,n1*n2)
        integer, intent(out) :: i1(n1*n2), i2(n1*n2)
        integer :: m, l, ngrid
        real    :: d2, s
        ngrid = n1 * n2
        do m = 1,ngrid
            i1(m) = (m - 1) / n2; i2(m) = mod(m - 1, n2)
        enddo
        do m = 1,ngrid
            s = 0.
            do l = 1,ngrid
                d2 = real((i1(l) - i1(m))**2 + (i2(l) - i2(m))**2)
                Ksm(l,m) = exp(-0.5 * d2 / max(sigma, 1.e-6)**2); s = s + Ksm(l,m)
            enddo
            Ksm(:,m) = Ksm(:,m) / s
        enddo
    end subroutine klem_lattice_kernel

    !> OT-balanced responsibilities (Sinkhorn on the log-likelihoods): grid columns -> equal shares of the grid mass,
    !! rows -> 1, free components unconstrained, dual potential tempered by lambda (the prototype's --balance)
    subroutine klem_sinkhorn( L, ncl, ngrid, n, lambda, r )
        integer,  intent(in)    :: ncl, ngrid, n
        real,     intent(in)    :: L(ncl,n)
        real,     intent(in)    :: lambda
        real,     intent(inout) :: r(ncl,n)
        real(dp), allocatable :: u(:), v(:), col(:)
        real(dp) :: gmass, tgt, mx, s
        integer  :: it, i, ic
        allocate(u(n), v(ncl), col(ncl)); u = 0.d0; v = 0.d0
        gmass = sum(real(r(1:ngrid,:),dp)); tgt = log(max(gmass, 1.d-30) / real(ngrid,dp))
        do it = 1,30
            do ic = 1,ngrid
                mx = -huge(1.d0)
                do i = 1,n
                    mx = max(mx, real(L(ic,i),dp) + u(i))
                enddo
                s = 0.d0
                do i = 1,n
                    s = s + exp(real(L(ic,i),dp) + u(i) - mx)
                enddo
                v(ic) = tgt - (mx + log(s))
            enddo
            !$omp parallel do default(shared) private(i,ic,mx,s,col)
            do i = 1,n
                do ic = 1,ncl
                    col(ic) = real(L(ic,i),dp) + v(ic)
                enddo
                mx = maxval(col); s = sum(exp(col - mx)); u(i) = -(mx + log(s))
            enddo
            !$omp end parallel do
        enddo
        !$omp parallel do default(shared) private(i,ic,mx,s,col)
        do i = 1,n
            do ic = 1,ncl
                col(ic) = real(L(ic,i),dp) + u(i) + real(lambda,dp) * v(ic)
            enddo
            mx = maxval(col); s = sum(exp(col - mx))
            do ic = 1,ncl
                r(ic,i) = real(exp(col(ic) - mx) / s)
            enddo
        enddo
        !$omp end parallel do
        deallocate(u, v, col)
    end subroutine klem_sinkhorn

    !> write a centred coarse-grid Fourier volume as a real-space MRC at the full (even) box, zero-padded in Fourier space
    subroutine klem_write_volume( Vc, kmax, box, smpd, fname )
        integer,          intent(in) :: kmax, box
        complex,          intent(in) :: Vc(-kmax:kmax,-kmax:kmax,-kmax:kmax)
        real,             intent(in) :: smpd
        character(len=*), intent(in) :: fname
        type(image)  :: vol
        type(string) :: fn
        integer :: h, k, l, phys(3)
        call vol%new([box,box,box], smpd)
        call vol%zero_and_flag_ft()
        do l = -kmax,kmax
            do k = -kmax,kmax
                do h = 0,kmax   ! the stored half: h >= 0 (the Friedel mates are implied)
                    phys = vol%comp_addr_phys(h,k,l)
                    call vol%set_fcomp([h,k,l], phys, Vc(h,k,l) / real(box)**3)   ! SIMPLE's inverse FFT has no 1/N
                enddo
            enddo
        enddo
        call vol%ifft()
        fn = trim(fname)
        call vol%write(fn)
        call vol%kill
    end subroutine klem_write_volume

    !> SIMPLE Euler angles from a rotation matrix (inverse of klem_rotmat): the prototype's m2euler_simple with the sign flips
    pure function klem_rotmat2euler( R ) result( e )
        real, intent(in) :: R(3,3)
        real :: e(3), e1, e2, e3
        e2 = acos(min(max(R(3,3), -1.), 1.)) * 180. / real(PI)
        if( abs(R(3,3)) < 1. - 1.e-6 )then
            e1 = atan2(-R(3,2), R(3,1)) * 180. / real(PI); e3 = atan2(-R(2,3), -R(1,3)) * 180. / real(PI)
        else
            e1 = atan2(R(2,1), R(1,1)) * 180. / real(PI) * merge(1., -1., R(3,3) > 0.); e3 = 0.
        endif
        e(1) = modulo(-e1, 360.); e(2) = e2; e(3) = modulo(-e3, 360.)
    end function klem_rotmat2euler

    !> rotation matrix of a rotation vector (Rodrigues)
    pure function klem_rotvec( v ) result( Rm )
        real, intent(in) :: v(3)
        real :: Rm(3,3), K(3,3), a, kx, ky, kz
        a = sqrt(sum(v**2)); Rm = 0.; Rm(1,1) = 1.; Rm(2,2) = 1.; Rm(3,3) = 1.
        if( a < 1.e-12 ) return
        kx = v(1) / a; ky = v(2) / a; kz = v(3) / a
        K = 0.; K(1,2) = -kz; K(1,3) = ky; K(2,1) = kz; K(2,3) = -kx; K(3,1) = -ky; K(3,2) = kx
        Rm = Rm + sin(a) * K + (1. - cos(a)) * matmul(K, K)
    end function klem_rotvec

    !> the whole-reference local pose step: every particle's rotation moved towards the posterior-mean rotation over a
    !! 3x3x3 rotation-vector grid (step, then step/2) under a Gaussian prior of width prior (deg), scored against its
    !! argmax component of the OTHER half (rotation only, no shifts); the prototype's _pose_chunk
    subroutine klem_pose_step( bd, kc, dat, V, ncl, conv, sig2, scl, am, step, prior, nmoved, gain )
        type(klem_band), intent(in)    :: bd
        integer,         intent(in)    :: kc, ncl, conv
        type(klem_data), intent(inout) :: dat
        complex,         intent(in)    :: V(-kc:kc,-kc:kc,-kc:kc,ncl,2)
        real,            intent(in)    :: sig2(bd%nband), scl(bd%nband)
        integer,         intent(in)    :: am(dat%n)
        real,            intent(in)    :: step, prior
        integer,         intent(out)   :: nmoved
        real,            intent(out)   :: gain
        complex, allocatable :: P(:)
        real     :: g(3,27), Rcur(3,3), Rc(3,3), L(27), Lp(27), wd(27), dv(3), stp, ll0, ll1, mx, sw
        real(dp) :: acc, gsum
        integer  :: i, q, ic, ih, ig, i1, i2, i3, pass, nm
        ig = 0
        do i1 = -1,1
            do i2 = -1,1
                do i3 = -1,1
                    ig = ig + 1; g(:,ig) = real([i1, i2, i3])
                enddo
            enddo
        enddo
        nm = 0; gsum = 0.d0
        !$omp parallel do default(shared) schedule(dynamic,8) proc_bind(close) private(i,q,ic,ih,ig,pass,P,Rcur,Rc,L,Lp,wd,dv,stp,ll0,ll1,mx,sw,acc) reduction(+:nm,gsum)
        do i = 1,dat%n
            allocate(P(bd%nband)); ic = am(i); ih = merge(2, 1, dat%eo(i) == 0); Rcur = dat%R(:,:,i); ll0 = 0.
            do pass = 1,2
                stp = merge(step, 0.5 * step, pass == 1) * real(PI) / 180.
                do ig = 1,27
                    Rc = matmul(Rcur, klem_rotvec(g(:,ig) * stp))
                    call klem_project_cell(bd, kc, V(:,:,:,ic,ih), Rc, conv, P)
                    acc = 0.d0
                    do q = 1,bd%nband
                        if( bd%shell(q) > kc ) cycle
                        acc = acc + real(bd%wq(q),dp) * real(abs(dat%F(q,i) - dat%C(q,i) * scl(q) * P(q)),dp)**2 / real(sig2(q),dp)
                    enddo
                    L(ig) = -real(acc)
                enddo
                if( pass == 1 ) ll0 = L(14)   ! the centre of the grid = the current pose
                do ig = 1,27
                    Lp(ig) = L(ig) - 0.5 * sum((g(:,ig) * stp)**2) / (prior * real(PI) / 180.)**2
                enddo
                mx = maxval(Lp); sw = 0.; dv = 0.
                do ig = 1,27
                    wd(ig) = exp(Lp(ig) - mx); sw = sw + wd(ig)
                enddo
                do ig = 1,27
                    dv = dv + wd(ig) / sw * g(:,ig) * stp
                enddo
                Rcur = matmul(Rcur, klem_rotvec(dv))
            enddo
            call klem_project_cell(bd, kc, V(:,:,:,ic,ih), Rcur, conv, P)
            acc = 0.d0
            do q = 1,bd%nband
                if( bd%shell(q) > kc ) cycle
                acc = acc + real(bd%wq(q),dp) * real(abs(dat%F(q,i) - dat%C(q,i) * scl(q) * P(q)),dp)**2 / real(sig2(q),dp)
            enddo
            ll1 = -real(acc)
            if( sum(abs(Rcur - dat%R(:,:,i))) > 1.e-6 ) nm = nm + 1
            gsum = gsum + real(ll1 - ll0, dp); dat%R(:,:,i) = Rcur
            deallocate(P)
        enddo
        !$omp end parallel do
        nmoved = nm; gain = real(gsum / real(max(dat%n,1),dp))
    end subroutine klem_pose_step

end module simple_klem_core
