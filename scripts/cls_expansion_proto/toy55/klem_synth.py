#!/usr/bin/env python3
"""Synthetic 2-D family for klem: blob B moves along y with z1 (5 levels, +-6 px = large motion), blob C is present or
absent with z2 (composition), 10 % junk (a random blob elsewhere). usage: klem_synth.py OUTDIR [--n 3000 --box 48 --snr 0.3]"""
import sys, os, argparse, numpy as np
os.environ.setdefault('TOY_SMPD', '1.464'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model_class import write_mrc, fftc
from latent_windows import slice_batch, shift_phase
from view_tangent_toy import euler2m, Band, ctf_eval
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--n', type=int, default=3000); ap.add_argument('--box', type=int, default=48); ap.add_argument('--snr', type=float, default=0.3)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); rng = np.random.default_rng(0); B = a.box; zz, yy, xx = np.mgrid[:B, :B, :B]
def blob(c, s=3.5): return np.exp(-0.5 * ((zz - c[0]) ** 2 + (yy - c[1]) ** 2 + (xx - c[2]) ** 2) / s ** 2)
def vol(z1, z2): return blob((24, 24, 24), 4.5) + 0.8 * blob((24, 30 + 3 * z1, 15)) + 0.8 * z2 * blob((18, 16, 30), 3.)
levels = [-2, -1, 0, 1, 2]; n = a.n; nj = int(0.1 * n); z1 = rng.choice(levels, n); z2 = rng.integers(0, 2, n); junk = np.zeros(n, bool); junk[rng.choice(n, nj, replace=False)] = True
F = {}; 
for l in levels:
    for c in (0, 1): F[(l, c)] = fftc(vol(l, c))
kmax = B // 2 - 2; Bd = Band(B, kmax); cc = B // 2; ix = (Bd.kx + cc).astype(int); iy = (Bd.ky + cc).astype(int)
e = np.c_[rng.uniform(0, 360, n), rng.uniform(0, 180, n), rng.uniform(0, 360, n)]; sx = rng.uniform(-2, 2, n); sy = rng.uniform(-2, 2, n); dfx = rng.uniform(1.0, 2.0, n)
imgs = np.zeros((n, B, B), np.float32); sig = None
for i in range(n):
    V = fftc(0.9 * blob((rng.uniform(14, 34), rng.uniform(14, 34), rng.uniform(14, 34)), 5.)) if junk[i] else F[(z1[i], z2[i])]
    R = euler2m([-e[i, 0], e[i, 1], -e[i, 2]])
    P = slice_batch(np.ascontiguousarray(V.real, np.float32), np.ascontiguousarray(V.imag, np.float32), R[None], Bd.kx, Bd.ky, 0)[0] * ctf_eval(Bd.kx, Bd.ky, dfx[i], dfx[i], 0., B) * shift_phase(Bd, sx[i:i + 1], sy[i:i + 1], B, 1., 1.)[0]
    Fi = np.zeros((B, B), complex); Fi[iy, ix] = P; Fi[(-Bd.ky + cc).astype(int), (-Bd.kx + cc).astype(int)] = np.conj(P)
    im = np.real(np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(Fi))))
    if sig is None: sig = im.std() / np.sqrt(a.snr)
    imgs[i] = im + sig * rng.standard_normal(im.shape)
h = np.zeros(256, np.int32); h[0:3] = (B, B, n); h[3] = 2; open(a.out + '/stk.mrcs', 'wb').write(h.tobytes() + imgs.tobytes())
with open(a.out + '/tab.txt', 'w') as f:
    for i in range(n): f.write('e1=%.3f e2=%.3f e3=%.3f x=%.3f y=%.3f dfx=%.4f dfy=%.4f angast=0 eo=%d state=1\n' % (e[i, 0], e[i, 1], e[i, 2], sx[i], sy[i], dfx[i], dfx[i], i % 2))
conf = np.where(junk, 0, 1 + (z1 + 2) + 5 * z2)
with open(a.out + '/gt.txt', 'w') as f:
    f.write('# particle conformer closed e1 e2 (closed = blob C present)\n'); [f.write('%d %d %d 0 0\n' % (i + 1, conf[i], int(z2[i] == 1 and not junk[i]))) for i in range(n)]
write_mrc(a.out + '/cons.mrc', sum(vol(l, c) for l in levels for c in (0, 1)) / 10., 1.464)
print('states 1..10 = (z1 level -2..2) + 5*z2; axis = z1; written', a.out, '; conformer counts', np.bincount(conf)[1:])
