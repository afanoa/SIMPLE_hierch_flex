#!/usr/bin/env python3
"""Before any tree on real data: build the bench's FFT cache for both halves (splitbench3 format: images x soft 230 A
mask, fft2, cropped to CROP) and reconstruct each half's consensus from NSUB random particles with the bench's insertion
(calibrated conv: R k, ramp +, signed or |CTF|); the cross-half FSC per shell must be high when poses, keys and CTF
handling are right. usage: geometry_check.py SIM_PREFIX LAT_PREFIX CACHE_PREFIX [NSUB]"""
import sys, os, numpy as np
simp, latp, cachep = sys.argv[1:4]; nsub = int(sys.argv[4]) if len(sys.argv) > 4 else 40000; c = 48; smpd = 1.464; kv_ = 300.; cs = 2.7; fraca = 0.1
lam = 12.2643247 / np.sqrt(kv_ * 1e3 * (1. + kv_ * 1e3 * 0.978466e-6)); rng = np.random.default_rng(0)
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
def rotm(e1, e2, e3):
    a1, a2, a3 = np.radians(e1), np.radians(e2), np.radians(e3); c1, si1, c2, si2, c3, si3 = np.cos(a1), np.sin(a1), np.cos(a2), np.sin(a2), np.cos(a3), np.sin(a3)
    Z1 = np.zeros((len(e1), 3, 3)); Z1[:, 0, 0] = c1; Z1[:, 0, 1] = -si1; Z1[:, 1, 0] = si1; Z1[:, 1, 1] = c1; Z1[:, 2, 2] = 1
    Y2 = np.zeros((len(e1), 3, 3)); Y2[:, 0, 0] = c2; Y2[:, 0, 2] = si2; Y2[:, 2, 0] = -si2; Y2[:, 2, 2] = c2; Y2[:, 1, 1] = 1
    Z3 = np.zeros((len(e1), 3, 3)); Z3[:, 0, 0] = c3; Z3[:, 0, 1] = -si3; Z3[:, 1, 0] = si3; Z3[:, 1, 1] = c3; Z3[:, 2, 2] = 1
    return Z1 @ Y2 @ Z3
res = {}
for half in 'AB':
    simdir = '%s_%s' % (simp, half); h = np.frombuffer(open(simdir + '/all.mrcs', 'rb').read(1024), np.int32); n, nz, nsym = int(h[0]), int(h[2]), int(h[23])
    x = np.arange(n) - n // 2; rr2 = np.hypot(x[:, None], x[None, :]); R = 230. / 2 / smpd; msk2 = np.clip((R + 4 - rr2) / 8., 0, 1).astype(np.float32); sl = np.r_[0:c // 2, n - c // 2:n]
    cache = '%s_%s_c%d.npy' % (cachep, half, c)
    if not os.path.isfile(cache):
        I = np.memmap(simdir + '/all.mrcs', np.float32, 'r', offset=1024 + nsym, shape=(nz, n, n)); Fall = np.zeros((nz, c, c), np.complex64)
        for b in range(0, nz, 4000):
            blk = np.array(I[b:b + 4000], np.float32) * msk2; Fall[b:b + len(blk)] = np.fft.fft2(np.fft.ifftshift(blk, axes=(1, 2)))[:, sl][:, :, sl]
        np.save(cache, Fall); del Fall; print('half %s: cache written (%d particles)' % (half, nz), flush=True)
    P = np.array([[float(kv(l).get(k, 0)) for k in ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast')] for l in open('%s_%s/ptcl3D_r0.txt' % (latp, half)) if l.strip()])
    sel = np.sort(rng.choice(len(P), min(nsub, len(P)), replace=False)); F = np.load(cache, mmap_mode='r')[sel]; pp = P[sel]
    fc = np.fft.fftfreq(c) * c / n; fy2, fx2 = np.meshgrid(fc, fc, indexing='ij'); kk0 = np.stack([fx2, fy2, np.zeros_like(fx2)], -1).reshape(-1, 3) * n
    kf = np.hypot(fx2, fy2) / smpd; ang = np.arctan2(fy2, fx2)
    df = 0.5 * ((pp[:, 5] + pp[:, 6])[:, None, None] + (pp[:, 5] - pp[:, 6])[:, None, None] * np.cos(2. * (ang[None] - np.radians(pp[:, 7])[:, None, None]))) * 1e4
    chi = np.pi * lam * df * kf[None] ** 2 - 0.5 * np.pi * cs * 1e7 * lam ** 3 * kf[None] ** 4; C = (np.sqrt(1. - fraca ** 2) * np.sin(chi) + fraca * np.cos(chi)).astype(np.float32)
    ramp = np.exp(2j * np.pi * (fx2[None] * pp[:, 3, None, None] + fy2[None] * pp[:, 4, None, None])).astype(np.complex64)
    q = np.round(np.einsum('qj,nij->nqi', kk0, rotm(pp[:, 0], pp[:, 1], pp[:, 2]))).astype(np.int32); ok = (np.abs(q) <= c // 2 - 1).all(2); q = q % c
    idx = (q[..., 2] * c * c + q[..., 1] * c + q[..., 0]).astype(np.int64); idx[~ok] = c ** 3; ii = idx.ravel(); m3 = c ** 3 + 1
    for mode, CC in (('signed', C), ('abs', np.abs(C))):
        Wv = (CC * np.conj(ramp) * F).reshape(len(sel), -1); Wd = (CC ** 2).reshape(len(sel), -1)
        num = np.bincount(ii, Wv.real.ravel(), m3)[:-1] + 1j * np.bincount(ii, Wv.imag.ravel(), m3)[:-1]; den = np.bincount(ii, Wd.ravel(), m3)[:-1]
        res[(half, mode)] = num / (den + 0.05 * np.median(den[den > 0]))
f3 = np.fft.fftfreq(c) * c / n; kx, ky, kz = np.meshgrid(f3, f3, f3, indexing='ij'); sh = np.minimum(np.round(np.sqrt(kx**2 + ky**2 + kz**2) * n).astype(int), c // 2).ravel(); ns = c // 2 + 1
for mode in ('signed', 'abs'):
    X, Y = res[('A', mode)], res[('B', mode)]
    fsc = np.bincount(sh, np.real(X * np.conj(Y)), ns) / np.maximum(np.sqrt(np.bincount(sh, np.abs(X) ** 2, ns) * np.bincount(sh, np.abs(Y) ** 2, ns)), 1e-30)
    print('cross-half FSC of the %s-CTF consensus (%d particles per half): %s' % (mode, nsub, '  '.join('%.1fA %.2f' % (n * smpd / k, fsc[k]) for k in (2, 4, 6, 8, 12, 16, 20, 23))), flush=True)
