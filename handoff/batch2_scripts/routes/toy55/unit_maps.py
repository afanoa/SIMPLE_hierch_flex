#!/usr/bin/env python3
"""Viewing maps of merged units: the particle-weighted average of each unit's member leaves' half-map means (unfiltered), padded in
Fourier space to BOX and low-passed to LP (raised-cosine edge), unsharpened: OUTDIR/unit_NN.mrc (full maps only).
usage: unit_maps.py OUTDIR LEAVES_INT (id even odd count) UNIT_LEAVES (unit leaf,leaf,..) [--box 160 --lp 10 --junk JUNKFILE]"""
import os, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('leaves'); ap.add_argument('units'); ap.add_argument('--box', type=int, default=160)
ap.add_argument('--lp', type=float, default=10.); ap.add_argument('--junk', default=''); a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64), float(np.fromfile(p, np.float32, 13)[10]) / n
def wr(p, x, smpd):
    h = np.zeros(256, np.int32); h[0:3] = x.shape; h[3] = 2; h[7:10] = x.shape; hf = h.view(np.float32); hf[10:13] = np.array(x.shape) * smpd; hf[13:16] = 90.
    h[16:19] = [1, 2, 3]; hf[19], hf[20], hf[21] = x.min(), x.max(), x.mean(); h[52] = int.from_bytes(b'MAP ', 'little'); h[53] = 0x4444
    with open(p, 'wb') as f: f.write(h.tobytes()); f.write(x.astype(np.float32).tobytes())
L = {int(l.split()[0]): (l.split()[1], l.split()[2], float(l.split()[3])) for l in open(a.leaves) if l.strip() and not l.startswith('#')}
U = [(int(l.split()[0]), [int(x) for x in l.split()[1].split(',')]) for l in open(a.units) if l.strip() and not l.startswith('#')]
junk = {int(x) for x in open(a.junk).read().split()} if a.junk and os.path.isfile(a.junk) else set()
for u, mem in U:
    acc, wsum, smpd = 0., 0., None
    for lf in mem:
        e, o, c = L[lf]; E, smpd = rd(e); O, _ = rd(o); acc = acc + c * (E + O) / 2; wsum += c
    V = acc / wsum; n = V.shape[0]
    if a.box > n:
        F = np.fft.fftshift(np.fft.fftn(V)); P = np.zeros((a.box,) * 3, complex); lo = a.box // 2 - n // 2; P[lo:lo + n, lo:lo + n, lo:lo + n] = F
        V = np.real(np.fft.ifftn(np.fft.ifftshift(P))) * (a.box / n) ** 3; smpd = smpd * n / a.box; n = a.box
    f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n), indexing='ij'); k = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2) / smpd
    k1, k2 = 1. / (a.lp + 1.), 1. / (a.lp - 0.5); filt = np.where(k <= k1, 1., np.where(k >= k2, 0., 0.5 + 0.5 * np.cos(np.pi * (k - k1) / (k2 - k1))))
    wr('%s/%sunit_%02d.mrc' % (a.out, 'junk_' if u in junk else '', u), np.fft.irfftn(np.fft.rfftn(V) * filt, s=V.shape), smpd)
print('%d unit maps in %s (%d junk)' % (len(U), a.out, len(junk)))
