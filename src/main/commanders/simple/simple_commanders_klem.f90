!@descr: kernel-latent EM (klem) E-step: per-particle log-likelihoods against ncl cell volumes of the OTHER data half
!! (cross-half), band-limited to lp, at the particle's own pose; the Fortran twin of klem.py's _estep_chunk so the two
!! can be compared number by number. The prototype's calibrated conventions come in as environment knobs:
!!   KLEM_SS (shift sign, -1), KLEM_FS (shift scale, 1), KLEM_CONV (0: loc = hk R, 1: loc = hk R^T), KLEM_SIGN (contrast
!!   sign, 1), KLEM_PHORIG (0: SIMPLE's fft already centres the phase origin; 1 adds numpy's checkerboard on top, for tests only), KLEM_PINDS (text file,
!!   one "row eo" pair per line: the particles to score and their half, 0 = even), KLEM_NSAMP (200: particles whose
!!   band coefficients are written for the comparison).
!! Inputs: projfile, dir (cell_NN_even.mrc / cell_NN_odd.mrc on the coarse grid 2*kmax+1), nstates (ncl), lp (A),
!! infile (per-shell noise and reference scale: "shell sig2 scale", shells 0..kmax), outfile (stream file: header,
!! ll(ncl) per particle, then the samples).
module simple_commanders_klem
use simple_core_module_api
use simple_commander_base, only: commander_base
use simple_cmdline,        only: cmdline
use simple_parameters,     only: parameters
use simple_sp_project,     only: sp_project
use simple_image,          only: image
use simple_string,         only: string
!$ use omp_lib
implicit none
public :: commander_klem_estep
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_klem_estep
  contains
    procedure :: execute => exec_klem_estep
end type commander_klem_estep

contains

    subroutine exec_klem_estep( self, cline )
        class(commander_klem_estep), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(parameters)         :: params
        type(sp_project)         :: spproj
        type(image)              :: vol
        type(image), allocatable :: imgs(:)
        type(string)             :: stkname, fn
        type(ctfparams)          :: ctfp
        complex,     allocatable :: V(:,:,:,:,:), Fs(:,:), Ps(:,:)
        real,        allocatable :: ll(:,:), Cs(:,:), kx(:), ky(:), wq(:), sig2(:), scl(:), sig2_sh(:), scl_sh(:)
        integer,     allocatable :: pinds(:), eos(:), shell(:), hs(:), ks(:)
        character(len=256)       :: envval, fname
        real(dp) :: acc, wl, cs_px, dfxp, dfyp, df, k2, ph, ang, phase, afac, two_pi
        real     :: sshift, fscale, csign, e(3), sh(2), R(3,3), M(3,3), loc(3), dfx, dfy, angast, smpd, smpd_c, rr
        integer  :: conv, phorig, nsamp, box, kmax, n3, c3, c2, nband, ncl, n, i, j, q, hh, kk, lz, ic, half, ind, u, io
        integer  :: nb, ib, i0, i1, t, nthr, ishell, ihalf, phys(3), envlen, envstat, lims(3,2), ldim(3), nimgs, nc
        complex  :: comp, F_q, P_q
        logical  :: l_exists
        two_pi = 2.d0 * PI
        call params%new(cline)
        call spproj%read(params%projfile)
        box  = params%box
        smpd = params%smpd
        kmax = floor(real(box) * smpd / params%lp)
        n3   = 2 * kmax + 1
        c3   = n3 / 2
        c2   = box / 2
        ncl    = params%nstates
        ! conventions
        sshift = -1.; fscale = 1.; conv = 0; csign = 1.; phorig = 0; nsamp = 200
        call get_environment_variable('KLEM_SS', envval, envlen, envstat);     if( envstat == 0 .and. envlen > 0 ) read(envval(1:envlen),*) sshift
        call get_environment_variable('KLEM_FS', envval, envlen, envstat);     if( envstat == 0 .and. envlen > 0 ) read(envval(1:envlen),*) fscale
        call get_environment_variable('KLEM_CONV', envval, envlen, envstat);   if( envstat == 0 .and. envlen > 0 ) read(envval(1:envlen),*) conv
        call get_environment_variable('KLEM_SIGN', envval, envlen, envstat);   if( envstat == 0 .and. envlen > 0 ) read(envval(1:envlen),*) csign
        call get_environment_variable('KLEM_PHORIG', envval, envlen, envstat); if( envstat == 0 .and. envlen > 0 ) read(envval(1:envlen),*) phorig
        call get_environment_variable('KLEM_NSAMP', envval, envlen, envstat);  if( envstat == 0 .and. envlen > 0 ) read(envval(1:envlen),*) nsamp
        ! the particles to score and their halves
        call get_environment_variable('KLEM_PINDS', envval, envlen, envstat)
        if( envstat /= 0 .or. envlen < 1 ) THROW_HARD('KLEM_PINDS (file of "row eo" lines) is required')
        open(newunit=u, file=envval(1:envlen), status='old', action='read')
        n = 0
        do
            read(u,*,iostat=io) i, j
            if( io /= 0 ) exit
            n = n + 1
        enddo
        rewind(u)
        allocate(pinds(n), eos(n))
        do i = 1,n
            read(u,*) pinds(i), eos(i)
        enddo
        close(u)
        ! the band: half plane to shell kmax, numpy meshgrid('xy') order (ky outer, kx inner), kx > 0 or (kx == 0, ky >= 0)
        nband = 0
        do kk = -kmax,kmax
            do hh = -kmax,kmax
                rr = sqrt(real(hh*hh + kk*kk))
                if( rr <= real(kmax) .and. (hh > 0 .or. (hh == 0 .and. kk >= 0)) ) nband = nband + 1
            enddo
        enddo
        allocate(hs(nband), ks(nband), kx(nband), ky(nband), wq(nband), shell(nband), sig2(nband), scl(nband))
        q = 0
        do kk = -kmax,kmax
            do hh = -kmax,kmax
                rr = sqrt(real(hh*hh + kk*kk))
                if( rr <= real(kmax) .and. (hh > 0 .or. (hh == 0 .and. kk >= 0)) )then
                    q = q + 1; hs(q) = hh; ks(q) = kk; kx(q) = real(hh); ky(q) = real(kk)
                    wq(q) = merge(2., 1., hh > 0); shell(q) = nint(rr)
                endif
            enddo
        enddo
        ! per-shell noise and reference scale
        allocate(sig2_sh(0:kmax), scl_sh(0:kmax), source=1.)
        open(newunit=u, file=params%infile%to_char(), status='old', action='read')
        do
            read(u,*,iostat=io) ishell, rr, ang
            if( io /= 0 ) exit
            if( ishell >= 0 .and. ishell <= kmax )then
                sig2_sh(ishell) = rr; scl_sh(ishell) = real(ang)
            endif
        enddo
        close(u)
        do q = 1,nband
            sig2(q) = sig2_sh(shell(q)); scl(q) = scl_sh(shell(q))
        enddo
        ! the cell volumes of both halves as centred complex grids (hh,kk,lz in -kmax..kmax)
        ! cell maps: any even box spanning the particle box's extent (the coarse grid zero-padded, or a box_crop
        ! reconstruction); the Fourier lattice is then the particle box's and the band fits inside it
        fname = trim(params%dir%to_char())//'/cell_01_even.mrc'; fn = trim(fname); call find_ldim_nptcls(fn, ldim, nimgs)
        nc = ldim(1)
        if( mod(nc,2) /= 0 ) THROW_HARD('cell maps must have an even box (zero-pad the coarse grid to the full box)')
        if( nc/2 < kmax ) THROW_HARD('cell map box too small for the band')
        smpd_c = real(box) * smpd / real(nc)
        allocate(V(-kmax:kmax,-kmax:kmax,-kmax:kmax,ncl,2), source=cmplx(0.,0.))
        call vol%new([nc,nc,nc], smpd_c)
        do half = 1,2
            do ic = 1,ncl
                fname = trim(params%dir%to_char())//'/cell_'//int2str_pad(ic,2)//merge('_even.mrc', '_odd.mrc ', half == 1)
                inquire(file=trim(fname), exist=l_exists)
                if( .not. l_exists ) THROW_HARD('missing cell volume: '//trim(fname))
                fn = trim(fname)
                call vol%read(fn)
                call vol%fft()
                do lz = -kmax,kmax
                    do kk = -kmax,kmax
                        do hh = -kmax,kmax
                            phys = vol%comp_addr_phys(hh,kk,lz)
                            comp = vol%get_fcomp([hh,kk,lz], phys)
                            if( phorig == 1 .and. mod(abs(hh + kk + lz), 2) == 1 ) comp = -comp   ! numpy-style centred origin on top of SIMPLE's
                            V(hh,kk,lz,ic,half) = comp * real(nc)**3   ! SIMPLE's forward FFT divides by N; the prototype's does not
                        enddo
                    enddo
                enddo
            enddo
        enddo
        call vol%kill
        write(logfhandle,'(A,I0,A,I0,A,I0,A,I0,A)') '>>> KLEM E-STEP: ', n, ' particles x ', ncl, ' cells, band kmax ', kmax, ' (', nband, ' coefficients)'
        write(logfhandle,'(A,F5.1,A,F5.2,A,I0,A,F5.1,A,I0)') '>>> KLEM conventions: shift sign ', sshift, ' scale ', fscale, ' conv ', conv, ' contrast sign ', csign, ' phase origin ', phorig
        ! outputs
        allocate(ll(ncl,n), source=0.)
        nsamp = min(nsamp, n)
        allocate(Fs(nband,nsamp), Ps(nband,nsamp), source=cmplx(0.,0.)); allocate(Cs(nband,nsamp), source=0.)
        nthr = 1
        !$ nthr = omp_get_max_threads()
        nb = 512
        allocate(imgs(nb))
        do i = 1,nb
            call imgs(i)%new([box,box,1], smpd)
        enddo
        do ib = 1,n,nb
            i0 = ib; i1 = min(ib + nb - 1, n)
            ! serial read of the batch
            do i = i0,i1
                call spproj%get_stkname_and_ind('ptcl3D', pinds(i), stkname, ind)
                call imgs(i-i0+1)%read(stkname, ind)
            enddo
            !$omp parallel do default(shared) schedule(dynamic) proc_bind(close) &
            !$omp private(i,t,q,ic,hh,kk,ihalf,ctfp,e,sh,R,M,loc,F_q,P_q,acc,wl,cs_px,dfxp,dfyp,df,k2,ph,ang,phase,afac)
            do i = i0,i1
                t = i - i0 + 1
                call imgs(t)%fft()
                ctfp = spproj%get_ctfparams('ptcl3D', pinds(i))
                e    = spproj%os_ptcl3D%get_euler(pinds(i))
                sh   = spproj%os_ptcl3D%get_2Dshift(pinds(i))
                ! R = rotmat(-e3, z) rotmat(e2, y) rotmat(-e1, z), the prototype's euler2m([-e1, e2, -e3])
                R = matmul(matmul(rotz(-e(3)), roty(e(2))), rotz(-e(1)))
                if( conv == 0 )then
                    M = R
                else
                    M = transpose(R)
                endif
                ! CTF constants in pixel units (the prototype's ctf_eval)
                wl    = 12.26d0 / sqrt(real(ctfp%kv,dp) * 1.d3 + 0.9784d0 * (real(ctfp%kv,dp) * 1.d3)**2 * 1.d-6) / real(smpd,dp)
                cs_px = real(ctfp%cs,dp) * 1.d7 / real(smpd,dp)
                dfxp  = real(ctfp%dfx,dp) * 1.d4 / real(smpd,dp)
                dfyp  = real(ctfp%dfy,dp) * 1.d4 / real(smpd,dp)
                afac  = atan(real(ctfp%fraca,dp) / sqrt(1.d0 - real(ctfp%fraca,dp)**2))
                ihalf = merge(2, 1, eos(i) == 0)   ! even particles (eo 0) are scored against the ODD cells and vice versa
                do ic = 1,ncl
                    acc = 0.d0
                    do q = 1,nband
                        hh = hs(q); kk = ks(q)
                        ! particle coefficient: fftshift-centred phase origin, contrast sign, the shift folded in
                        F_q = imgs(t)%get_fcomp2D(hh,kk) * real(box)**2   ! undo SIMPLE's 1/N forward-FFT normalisation
                        if( phorig == 1 .and. mod(abs(hh + kk), 2) == 1 ) F_q = -F_q
                        phase = two_pi * (real(hh,dp) * real(sshift*fscale*sh(1),dp) + real(kk,dp) * real(sshift*fscale*sh(2),dp)) / real(box,dp)
                        F_q   = F_q * csign * cmplx(real(cos(phase)), real(sin(phase)))
                        ! CTF
                        ang = atan2(real(kk,dp), real(hh,dp))
                        df  = 0.5d0 * (dfxp + dfyp + cos(2.d0 * (ang - real(ctfp%angast,dp) * PI / 180.d0)) * (dfxp - dfyp))
                        k2  = real(hh*hh + kk*kk, dp) / real(box,dp)**2
                        ph  = PI * wl * k2 * (df - 0.5d0 * wl * wl * k2 * cs_px) + afac
                        ! projection of the cell at the particle pose: loc = hk M, trilinear on the centred grid
                        loc(1) = kx(q) * M(1,1) + ky(q) * M(2,1)
                        loc(2) = kx(q) * M(1,2) + ky(q) * M(2,2)
                        loc(3) = kx(q) * M(1,3) + ky(q) * M(2,3)
                        P_q = trilinear(V(:,:,:,ic,ihalf), kmax, loc)
                        acc = acc + real(wq(q),dp) * (real(abs(F_q - real(sin(ph)) * scl(q) * P_q),dp))**2 / real(sig2(q),dp)
                        if( i <= nsamp .and. ic == 1 )then
                            Fs(q,i) = F_q; Cs(q,i) = real(sin(ph)); Ps(q,i) = P_q
                        endif
                    enddo
                    ll(ic,i) = -real(acc)
                enddo
            enddo
            !$omp end parallel do
        enddo
        ! stream output: header, ll per particle, samples, the centred grid of cell 1 (even) in numpy [lz,kk,hh] order
        open(newunit=u, file=params%outfile%to_char(), access='stream', form='unformatted', status='replace', action='write')
        write(u) n, ncl, nband, nsamp, n3
        do i = 1,n
            write(u) ll(:,i)
        enddo
        do i = 1,nsamp
            write(u) Fs(:,i)
        enddo
        do i = 1,nsamp
            write(u) Cs(:,i)
        enddo
        do i = 1,nsamp
            write(u) Ps(:,i)
        enddo
        do lz = -kmax,kmax
            do kk = -kmax,kmax
                do hh = -kmax,kmax
                    write(u) V(hh,kk,lz,1,1)
                enddo
            enddo
        enddo
        close(u)
        write(logfhandle,'(A)') '>>> KLEM E-STEP written to '//params%outfile%to_char()
        do i = 1,nb
            call imgs(i)%kill
        enddo
        deallocate(imgs, V, ll, Fs, Ps, Cs)
        call spproj%kill
        call simple_end('**** SIMPLE_KLEM_ESTEP NORMAL STOP ****')

      contains

        pure function rotz( ang ) result( Rm )
            real, intent(in) :: ang
            real :: Rm(3,3), ca, sa
            ca = cos(ang * real(PI) / 180.); sa = sin(ang * real(PI) / 180.)
            Rm = 0.; Rm(1,1) = ca; Rm(1,2) = -sa; Rm(2,1) = sa; Rm(2,2) = ca; Rm(3,3) = 1.
        end function rotz

        pure function roty( ang ) result( Rm )
            real, intent(in) :: ang
            real :: Rm(3,3), ca, sa
            ca = cos(ang * real(PI) / 180.); sa = sin(ang * real(PI) / 180.)
            Rm = 0.; Rm(1,1) = ca; Rm(1,3) = -sa; Rm(3,1) = sa; Rm(3,3) = ca; Rm(2,2) = 1.
        end function roty

        !> trilinear interpolation on the centred grid with the coordinates clamped to the grid (scipy mode='nearest')
        pure function trilinear( Vc, kmax, loc ) result( val )
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
        end function trilinear

    end subroutine exec_klem_estep

end module simple_commanders_klem
