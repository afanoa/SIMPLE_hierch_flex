#!/usr/bin/env python3
"""Convert the Fortran klem_board outputs (klem_resp.bin, klem_halves.bin) into the npz files the Python tooling reads
(klem_resp.npz: pinds, r, i1, i2; klem_halves.npz: D0, D1, w, i1, i2, n3, band at <= 12 A as klem.py writes them).
usage: klem_fortran2npz.py BOARD_DIR"""
import sys, numpy as np
d = sys.argv[1]
with open(d + '/klem_resp.bin', 'rb') as f:
    n, ncl, ngrid, n1, n2 = np.fromfile(f, np.int32, 5); pinds = np.fromfile(f, np.int32, n); i1 = np.fromfile(f, np.int32, ngrid); i2 = np.fromfile(f, np.int32, ngrid)
    r = np.fromfile(f, np.float32, n * ncl).reshape(n, ncl)
np.savez_compressed(d + '/klem_resp.npz', pinds=pinds, r=r, i1=i1, i2=i2)
with open(d + '/klem_halves.bin', 'rb') as f:
    ncl2, ngrid2, kmax, n3 = np.fromfile(f, np.int32, 4); w = np.fromfile(f, np.float32, ngrid2)
    D = np.fromfile(f, np.complex64, 2 * ngrid2 * n3 ** 3).reshape(2, ngrid2, n3, n3, n3)   # fortran order h fastest -> numpy [l, k, h] = klem.py's (z, y, x) layout
kk = np.arange(n3) - n3 // 2; rr = np.sqrt(sum(q ** 2 for q in np.meshgrid(kk, kk, kk, indexing='ij'))); bm = (rr <= n3 // 2) & (rr >= 1); band = np.where(bm.ravel())[0]
D0 = np.array([D[0, k].ravel()[band] for k in range(ngrid2)]).astype(np.complex64); D1 = np.array([D[1, k].ravel()[band] for k in range(ngrid2)]).astype(np.complex64)
np.savez_compressed(d + '/klem_halves.npz', D0=D0, D1=D1, w=w, i1=i1, i2=i2, n3=n3, band=band)
print('%s: klem_resp.npz (%d particles x %d components, grid %d x %d) and klem_halves.npz (%d cells, %d band coefficients on the %d^3 grid)' % (d, n, ncl, n1, n2, ngrid2, len(band), n3))
