#!/usr/bin/env python3
"""Adds the 'specmatch 0|1' board option to simple_commanders_klem_board.f90: after the M-step every grid cell's radial
amplitude spectrum is matched to its half's consensus, so sharpness cannot be a latent axis. usage: patch_specmatch.py FILE"""
import sys
fn = sys.argv[1]; s = open(fn).read()
def rep(old, new, count=1):
    global s
    assert s.count(old) == count, (old[:60], s.count(old))
    s = s.replace(old, new)
if 'specmatch' in s: print('already patched'); sys.exit(0)
rep("!!   devw 0 (0.05 = deviation-FSC shell weights) | pnoise 0|1 (per-particle noise factor) | init_file none|FILE (pind g1 g2)",
    "!!   devw 0 (0.05 = deviation-FSC shell weights) | pnoise 0|1 (per-particle noise factor) | init_file none|FILE (pind g1 g2)\n!!   specmatch 0|1 (every grid cell's radial amplitude spectrum matched to its half's consensus after the M-step)")
rep("        integer  :: pnoise, pscale, nchg\n", "        integer  :: pnoise, pscale, nchg, specmatch\n")
rep("pscale = 0; stop_frac = 0.03;", "pscale = 0; specmatch = 0; stop_frac = 0.03;")
rep("                case('pscale');        read(tok(2),*) pscale\n", "                case('pscale');        read(tok(2),*) pscale\n                case('specmatch');     read(tok(2),*) specmatch\n")
rep("""                    V(:,:,:,ic,ih) = cmplx((num(:,:,:,ic,ih) + real(eps,dp) * Gc(:,:,:,ih)) / (den(:,:,:,ic,ih) + real(eps,dp)))
                enddo
            enddo
""", """                    V(:,:,:,ic,ih) = cmplx((num(:,:,:,ic,ih) + real(eps,dp) * Gc(:,:,:,ih)) / (den(:,:,:,ic,ih) + real(eps,dp)))
                enddo
            enddo
            if( specmatch /= 0 ) call match_cell_spectra()
""")
rep("        subroutine apply_focus()   ! the distributed board's parallel version, one image per thread\n", """        subroutine match_cell_spectra()   ! grid cells cannot differ from the consensus by sharpness: per-shell amplitude match (0.5..2 clip), free components untouched
            real(dp) :: pv(0:kmax), pg(0:kmax)
            real     :: w(0:kmax)
            integer  :: ic, ih, h, k, l, ish
            do ih = 1,2
                pg = 0.d0
                do l = -kmax,kmax
                    do k = -kmax,kmax
                        do h = -kmax,kmax
                            ish = nint(sqrt(real(h*h + k*k + l*l)))
                            if( ish <= kmax ) pg(ish) = pg(ish) + real(abs(Gc(h,k,l,ih)),dp)**2
                        enddo
                    enddo
                enddo
                !$omp parallel do default(shared) private(ic,h,k,l,ish,pv,w)
                do ic = 1,ngrid
                    pv = 0.d0
                    do l = -kmax,kmax
                        do k = -kmax,kmax
                            do h = -kmax,kmax
                                ish = nint(sqrt(real(h*h + k*k + l*l)))
                                if( ish <= kmax ) pv(ish) = pv(ish) + real(abs(V(h,k,l,ic,ih)),dp)**2
                            enddo
                        enddo
                    enddo
                    w = 1.
                    do ish = 2,kmax
                        if( pv(ish) > 0.d0 .and. pg(ish) > 0.d0 ) w(ish) = real(min(max(sqrt(pg(ish) / pv(ish)), 0.5d0), 2.d0))
                    enddo
                    do l = -kmax,kmax
                        do k = -kmax,kmax
                            do h = -kmax,kmax
                                ish = nint(sqrt(real(h*h + k*k + l*l)))
                                if( ish <= kmax ) V(h,k,l,ic,ih) = V(h,k,l,ic,ih) * w(ish)
                            enddo
                        enddo
                    enddo
                enddo
                !$omp end parallel do
            enddo
        end subroutine match_cell_spectra

        subroutine apply_focus()   ! the distributed board's parallel version, one image per thread
""")
open(fn, 'w').write(s); print('patched', fn)
