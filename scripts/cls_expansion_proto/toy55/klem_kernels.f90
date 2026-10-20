! OpenMP kernels for the kernel-latent EM prototype (klem.py), called through ctypes.
! Layouts follow numpy C-order arrays as passed: F(nc, np) <- numpy (np, nc); R(3,3,np) <- numpy (np,3,3) so that
! R(a, b, p) = numpy R[p, b-1, a-1]; refs(n3,n3,n3,K,2) <- numpy (2, K, n3, n3, n3) [h][kk][iz][iy][ix];
! ll(K, np) <- numpy (np, K); W(K, np) <- numpy (np, K); num/den(n3,n3,n3,K,2) <- numpy (2, K, n3, n3, n3).
! Build: gfortran -O3 -fopenmp -shared -fPIC -o libklem.so klem_kernels.f90
module klem_kernels
    use iso_c_binding
    use omp_lib
    implicit none
contains

    subroutine klem_set_threads(n) bind(C, name='klem_set_threads')
        integer(c_int), value :: n
        call omp_set_num_threads(n)
    end subroutine klem_set_threads

    !> trilinear interpolation of a 3D array at 0-based fractional coordinates, clamped (scipy order=1, mode='nearest')
    pure subroutine interp3(V, n3, x, y, z, val)
        integer(c_int), intent(in) :: n3
        real(c_float),  intent(in) :: V(n3, n3, n3)
        real(c_double), intent(in) :: x, y, z
        real(c_double), intent(out) :: val
        integer :: ix, iy, iz
        real(c_double) :: xc, yc, zc, fx, fy, fz
        xc = min(max(x, 0.0d0), dble(n3 - 1)); yc = min(max(y, 0.0d0), dble(n3 - 1)); zc = min(max(z, 0.0d0), dble(n3 - 1))
        ix = min(int(xc), n3 - 2); iy = min(int(yc), n3 - 2); iz = min(int(zc), n3 - 2)
        fx = xc - ix; fy = yc - iy; fz = zc - iz
        val = (1 - fz) * ((1 - fy) * ((1 - fx) * V(ix + 1, iy + 1, iz + 1) + fx * V(ix + 2, iy + 1, iz + 1)) &
                        + fy * ((1 - fx) * V(ix + 1, iy + 2, iz + 1) + fx * V(ix + 2, iy + 2, iz + 1))) &
            + fz * ((1 - fy) * ((1 - fx) * V(ix + 1, iy + 1, iz + 2) + fx * V(ix + 2, iy + 1, iz + 2)) &
                    + fy * ((1 - fx) * V(ix + 1, iy + 2, iz + 2) + fx * V(ix + 2, iy + 2, iz + 2)))
    end subroutine interp3

    !> 3D frequency of the 2D coefficient (kx, ky) under the particle rotation: loc = hk . M, M = R (conv 0) or R^T (conv 1)
    pure subroutine rotate_coef(R, conv, kx, ky, x, y, z)
        real(c_double), intent(in) :: R(3, 3)      ! R(a, b) = numpy R[b-1, a-1]
        integer(c_int), intent(in) :: conv
        real(c_double), intent(in) :: kx, ky
        real(c_double), intent(out) :: x, y, z
        if (conv == 0) then   ! loc_i = sum_j hk_j M[j, i] with M = R: M[j, i] = numpy R[j, i] = R(i+1, j+1)
            x = kx * R(1, 1) + ky * R(1, 2); y = kx * R(2, 1) + ky * R(2, 2); z = kx * R(3, 1) + ky * R(3, 2)
        else                  ! M = R^T: M[j, i] = numpy R[i, j] = R(j+1, i+1)
            x = kx * R(1, 1) + ky * R(2, 1); y = kx * R(1, 2) + ky * R(2, 2); z = kx * R(1, 3) + ky * R(2, 3)
        end if
    end subroutine rotate_coef

    !> E-step: ll(kk, p) = -sum_q wq |F - C . slice_k(R_p) . SC|^2 / sig2 against the references of the OTHER half
    subroutine klem_estep(np, K, nc, n3, fre, fim, C, R, eo, refre, refim, kx, ky, wq, sig2, SC, conv, ll) bind(C, name='klem_estep')
        integer(c_int), value :: np, K, nc, n3, conv
        real(c_double), intent(in) :: fre(nc, np), fim(nc, np), C(nc, np), R(3, 3, np)
        integer(c_int), intent(in) :: eo(np)
        real(c_float),  intent(in) :: refre(n3, n3, n3, K, 2), refim(n3, n3, n3, K, 2)
        real(c_double), intent(in) :: kx(nc), ky(nc), wq(nc), sig2(nc), SC(nc)
        real(c_double), intent(out) :: ll(K, np)
        integer :: p, kk, q, hr
        real(c_double) :: x, y, z, pre, pim, s, cen, pr, pi
        cen = dble(n3 / 2)
        !$omp parallel do private(p, kk, q, hr, x, y, z, pre, pim, s, pr, pi) schedule(dynamic, 8)
        do p = 1, np
            hr = 2 - eo(p)   ! numpy half index 1 - eo -> Fortran 2 - eo
            do kk = 1, K
                s = 0.0d0
                do q = 1, nc
                    call rotate_coef(R(:, :, p), conv, kx(q), ky(q), x, y, z)
                    call interp3(refre(:, :, :, kk, hr), n3, x + cen, y + cen, z + cen, pre)
                    call interp3(refim(:, :, :, kk, hr), n3, x + cen, y + cen, z + cen, pim)
                    pr = C(q, p) * pre * SC(q); pi = C(q, p) * pim * SC(q)
                    s = s + wq(q) * ((fre(q, p) - pr) ** 2 + (fim(q, p) - pi) ** 2) / sig2(q)
                end do
                ll(kk, p) = -s
            end do
        end do
        !$omp end parallel do
    end subroutine klem_estep

    !> M-step accumulation (in place): nearest-voxel backprojection with Friedel mates, weights W(kk, p), into the
    !! particle's OWN half; parallel over components (no reduction)
    subroutine klem_mstep(np, K, nc, n3, kmax, fre, fim, C, R, eo, W, kx, ky, conv, numre, numim, den) bind(C, name='klem_mstep')
        integer(c_int), value :: np, K, nc, n3, kmax, conv
        real(c_double), intent(in) :: fre(nc, np), fim(nc, np), C(nc, np), R(3, 3, np), W(K, np), kx(nc), ky(nc)
        integer(c_int), intent(in) :: eo(np)
        real(c_double), intent(inout) :: numre(n3, n3, n3, K, 2), numim(n3, n3, n3, K, 2), den(n3, n3, n3, K, 2)
        integer :: p, kk, q, h, ix, iy, iz, jx, jy, jz
        real(c_double) :: x, y, z, wt, cf_re, cf_im, c2
        !$omp parallel do private(p, kk, q, h, ix, iy, iz, jx, jy, jz, x, y, z, wt, cf_re, cf_im, c2) schedule(dynamic, 1)
        do kk = 1, K
            do p = 1, np
                wt = W(kk, p)
                if (wt < 1.0d-4) cycle
                h = eo(p) + 1
                do q = 1, nc
                    call rotate_coef(R(:, :, p), conv, kx(q), ky(q), x, y, z)
                    ix = nint(x) + kmax + 1; iy = nint(y) + kmax + 1; iz = nint(z) + kmax + 1
                    cf_re = wt * C(q, p) * fre(q, p); cf_im = wt * C(q, p) * fim(q, p); c2 = wt * C(q, p) * C(q, p)
                    if (ix >= 1 .and. ix <= n3 .and. iy >= 1 .and. iy <= n3 .and. iz >= 1 .and. iz <= n3) then
                        numre(ix, iy, iz, kk, h) = numre(ix, iy, iz, kk, h) + cf_re
                        numim(ix, iy, iz, kk, h) = numim(ix, iy, iz, kk, h) + cf_im
                        den(ix, iy, iz, kk, h) = den(ix, iy, iz, kk, h) + c2
                    end if
                    jx = -nint(x) + kmax + 1; jy = -nint(y) + kmax + 1; jz = -nint(z) + kmax + 1   ! Friedel mate: conjugate
                    if (jx >= 1 .and. jx <= n3 .and. jy >= 1 .and. jy <= n3 .and. jz >= 1 .and. jz <= n3) then
                        numre(jx, jy, jz, kk, h) = numre(jx, jy, jz, kk, h) + cf_re
                        numim(jx, jy, jz, kk, h) = numim(jx, jy, jz, kk, h) - cf_im
                        den(jx, jy, jz, kk, h) = den(jx, jy, jz, kk, h) + c2
                    end if
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine klem_mstep

end module klem_kernels
