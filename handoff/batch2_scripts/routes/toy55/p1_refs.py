#!/usr/bin/env python3
"""References for held-out scoring of partition maps (problem 1). For every map k of the --rec dirs (in order; unfiltered half-maps
E_k, O_k, soft sphere): FSC_k per shell and the posterior-mean references of ONE half, w = max(FSC, 0): refsE/vol_JJ = w O_k scores the
even particles, refsO/vol_JJ = w E_k the odd ones (no particle meets a map that contains it). --self N: also refsS/vol_JJ for the
first N maps = the self-inclusive control, w2 (E + O) / 2 with w2 = 2 FSC / (1 + FSC) (the average map, as refine3D with lp set uses).
refine3D with lp set matches both halves against volN itself (read_mask_filter_refvols), so these files are passed as volN.
usage: p1_refs.py OUT --rec DIR [--rec DIR ...] --nb 8 [--self 8 --mskdiam 220]"""
import os, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--rec', action='append', required=True); ap.add_argument('--nb', type=int, default=8)
ap.add_argument('--self', type=int, default=0); ap.add_argument('--mskdiam', type=float, default=220.); a = ap.parse_args()
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
def wr(p, x):
    h = np.zeros(256, np.int32); h[0:3] = x.shape; h[3] = 2; h[7:10] = x.shape; hf = h.view(np.float32); hf[10:13] = np.array(x.shape) * smpd; hf[13:16] = 90.
    h[16:19] = [1, 2, 3]; hf[19], hf[20], hf[21] = x.min(), x.max(), x.mean(); h[52] = int.from_bytes(b'MAP ', 'little'); h[53] = 0x4444
    with open(p, 'wb') as f: f.write(h.tobytes()); f.write(x.astype(np.float32).tobytes())
for d in ('refsE', 'refsO', 'refsS'): os.makedirs('%s/%s' % (a.out, d), exist_ok=True)
f0 = '%s/recvol_state01_even_unfil.mrc' % a.rec[0]; n = int(np.fromfile(f0, np.int32, 1)[0]); smpd = float(np.fromfile(f0, np.float32, 13)[10]) / n
g = np.arange(n) - n // 2; R = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2); r0 = a.mskdiam / 2 / smpd
M = np.where(R <= r0, 1., np.where(R <= r0 + 6, 0.5 + 0.5 * np.cos(np.pi * (R - r0) / 6), 0.))
f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n), indexing='ij'); sh = np.minimum(np.round(np.sqrt(kx ** 2 + ky ** 2 + kz ** 2) * n).astype(int), n // 2)
wgt = np.where((kz == 0) | (np.abs(kz) == 0.5), 1., 2.); res = n * smpd / np.maximum(np.arange(n // 2 + 1), 1)
j = 0; rows = []
for d in a.rec:
    for k in range(1, a.nb + 1):
        j += 1; EF = np.fft.rfftn(M * rd('%s/recvol_state%02d_even_unfil.mrc' % (d, k))); OF = np.fft.rfftn(M * rd('%s/recvol_state%02d_odd_unfil.mrc' % (d, k)))
        num = np.bincount(sh.ravel(), (wgt * np.real(EF * OF.conj())).ravel(), n // 2 + 1); den = np.sqrt(np.bincount(sh.ravel(), (wgt * np.abs(EF) ** 2).ravel(), n // 2 + 1) * np.bincount(sh.ravel(), (wgt * np.abs(OF) ** 2).ravel(), n // 2 + 1))
        fsc = num / np.maximum(den, 1e-300); fsc[0] = 1.; w = np.clip(fsc, 0, 1)[sh]
        wr('%s/refsE/vol_%02d.mrc' % (a.out, j), np.fft.irfftn(w * OF, s=(n, n, n))); wr('%s/refsO/vol_%02d.mrc' % (a.out, j), np.fft.irfftn(w * EF, s=(n, n, n)))
        if j <= a.self:
            w2 = np.clip(2 * fsc / (1 + fsc), 0, 1)[sh]; wr('%s/refsS/vol_%02d.mrc' % (a.out, j), np.fft.irfftn(w2 * (EF + OF) / 2, s=(n, n, n)))
        b = np.where(fsc[1:] < 0.5)[0]; b2 = np.where(fsc[1:] < 0.143)[0]
        rows.append('%d:%.1f/%.1f' % (j, res[b[0]] if len(b) else res[-1], res[b2[0]] if len(b2) else res[-1]))
print('%d references at box %d (%.3f A/px); FSC 0.5 / 0.143 per map: %s' % (j, n, smpd, ' '.join(rows)))
