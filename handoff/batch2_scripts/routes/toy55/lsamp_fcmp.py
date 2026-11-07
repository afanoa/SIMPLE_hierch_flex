#!/usr/bin/env python3
"""lsamp_fcmp.py: the Fortran landscape_sample outputs of one node against the Python prototype's pole-pair coordinates on the same
particles: |corr| of every Fortran pair coordinate (full and per fold) with every Python pair coordinate (lsamp_pq_coords.py), the
Fortran sample counts, and its held-out maps zero-padded to box NB (Fourier) for arm_readout.py.
usage: lsamp_fcmp.py FORTRAN_DIR PY_PQ_DIR RUN NB SMPD"""
import sys, os, glob, numpy as np
fd, pq, run, NB, smpd = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), float(sys.argv[5])
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
def write_mrc(fn, v, apix):
    v = np.ascontiguousarray(v, np.float32); nz, ny, nx = v.shape; h = np.zeros(256, np.int32); h[0:3] = [nx, ny, nz]; h[3] = 2; h[7:10] = [nx, ny, nz]
    hf = h.view(np.float32); hf[10:13] = [nx * apix, ny * apix, nz * apix]; hf[13:16] = 90.; h[16:19] = [1, 2, 3]; hf[19], hf[20], hf[21] = v.min(), v.max(), v.mean()
    h[52] = np.frombuffer(b'MAP ', np.int32)[0]; h[53] = np.frombuffer(bytes([0x44, 0x44, 0, 0]), np.int32)[0]
    with open(fn, 'wb') as f: f.write(h.tobytes()); f.write(v.tobytes())
def pad(v, nb):
    n = v.shape[0]; F = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v))); G = np.zeros((nb, nb, nb), complex); o = (nb - n) // 2; G[o:o + n, o:o + n, o:o + n] = F
    return np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(G)))) * (nb / n) ** 3
PY = {}
for t in ('c_oA', 'c_oB', 'oA_oB'):
    f = '%s/%s/%s/coords.npz' % (pq, t, run)
    if os.path.exists(f): PY[t] = np.load(f)
for line in open(fd + '/landscape_pairs.txt'):
    w = line.split(); a, b = w[0], w[1]
    print('== %s pair %s-%s: %s samples (folds %s / %s)%s' % (run, a, b, w[2], w[3], w[4], '; steps CV %s / %s' % (w[5], w[6]) if len(w) > 6 else ''))
    C = np.loadtxt('%s/landscape_%s_%s_coords.txt' % (fd, a, b)); r = C[:, 0].astype(int) - 1
    for t, P in PY.items():
        pr = P['rows']; com, i1, i2 = np.intersect1d(r, pr, return_indices=True)
        cs = [np.corrcoef(C[i1, k], P[key][i2, 0])[0, 1] for k, key in ((1, 'tfull'), (2, 'tA'), (3, 'tB'))]
        print('   vs Python %-6s (%d shared): corr full %+.3f, fold A %+.3f, fold B %+.3f' % (t, len(com), *cs))
    for sub in ('heldout_A', 'heldout_B', 'final'):
        D = '%s/landscape_%s_%s/%s' % (fd, a, b, sub); O = D + '_pad'; os.makedirs(O, exist_ok=True)
        for fn in sorted(glob.glob(D + '/state_*.mrc')): write_mrc(O + '/' + os.path.basename(fn), pad(rd(fn), NB), smpd)
