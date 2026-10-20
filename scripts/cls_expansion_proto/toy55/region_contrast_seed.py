#!/usr/bin/env python3
"""2D-route seeding of a REGION board: per direction bin, the principal within-bin contrast of the sub-class averages
inside the region's projected mask (pooled signal, one view); the sign of that contrast synchronised between neighbouring
bins (Z2 spectral sync on in-plane-calibrated comparisons); the signed contrasts back-projected along the known bin
directions into a 3D deviation mode of the region; every particle gets the coordinate of its sub-class along its bin's
contrast (signed), which seeds the board's lattice (init file "pind g1 g2") and is read out against the truth directly.
usage: region_contrast_seed.py STITCH_DIR CLSEXP_DIR BINS_NPZ CAVGS_BINS REGION_MASK CONSENSUS OUT [--nn 6 --maxang 25 --crop 80 --lp 12 --hp 60 --n1 7]"""
import sys, os, argparse, numpy as np
from scipy.ndimage import rotate as ndrotate, affine_transform
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc, write_mrc
ap = argparse.ArgumentParser(); [ap.add_argument(k) for k in ('stitch', 'clsexp', 'bins', 'cavgs', 'mask', 'cons', 'out')]
ap.add_argument('--nn', type=int, default=6); ap.add_argument('--maxang', type=float, default=25.); ap.add_argument('--crop', type=int, default=80); ap.add_argument('--lp', type=float, default=12.); ap.add_argument('--hp', type=float, default=60.)
ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--n1', type=int, default=7); ap.add_argument('--n2', type=int, default=1); ap.add_argument('--npart', type=int, default=0); a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
def read_stack(fn):
    h = np.frombuffer(open(fn, 'rb').read(1024), np.int32); nx, ny, nz, nsym = [int(v) for v in (h[0], h[1], h[2], h[23])]
    return np.array(np.memmap(fn, np.float32, 'r', offset=1024 + nsym, shape=(nz, ny, nx)), np.float32)
S = read_stack(a.stitch + '/avgs.mrcs'); M = np.load(a.stitch + '/avgs_meta.npz'); N = int(M['N']); full, even, odd = S[:N], S[N:2 * N], S[2 * N:3 * N]; pop = M['pop'].astype(float); parent = M['parent'].astype(int)
B = np.load(a.bins); C = B['centres']; nb = len(C); binavg = read_stack(a.cavgs); n = full.shape[1]; c = a.crop; sl = np.r_[0:c // 2, n - c // 2:n]
E2 = np.degrees(np.arccos(np.clip(C[:, 2], -1, 1))); E1 = np.degrees(np.arctan2(C[:, 1], C[:, 0])) % 360.
def rot(e1d, e2d):
    a1, a2 = np.radians([e1d, e2d]); c1, s1, c2, s2 = np.cos(a1), np.sin(a1), np.cos(a2), np.sin(a2)
    return np.array([[c1, -s1, 0], [s1, c1, 0], [0, 0, 1]]) @ np.array([[c2, 0, s2], [0, 1, 0], [-s2, 0, c2]])
def to_box(vol, nn):   # Fourier-crop a cubic volume to box nn (the averages' box) if needed
    m = vol.shape[0]
    if m == nn: return vol
    F = np.fft.fftn(np.fft.ifftshift(vol)); s_ = np.r_[0:nn // 2, m - nn // 2:m]; return np.fft.fftshift(np.real(np.fft.ifftn(F[np.ix_(s_, s_, s_)]))).astype(np.float32) * (nn / m) ** 3
V = np.clip(to_box(read_mrc(a.mask).astype(np.float32), n), 0, None); ctr = np.array([n // 2] * 3)
def project(vol, Rm):   # validated convention (diag_geom: corr 0.93 with the bin averages): rotate so the direction becomes axis 0, sum
    Mz = Rm[::-1, ::-1]; return affine_transform(vol, Mz, offset=ctr - Mz @ ctr, order=1).sum(0)
f = np.fft.fftfreq(n); fy, fx = np.meshgrid(f, f, indexing='ij'); k = np.hypot(fx, fy) / a.smpd
bp = (np.exp(-0.5 * (k * a.lp / 1.2) ** 2 * 2.) * (1 - np.exp(-0.5 * (k * a.hp / 1.2) ** 2 * 2.))).astype(np.float32)
def prep(img): F = np.fft.fft2(np.fft.ifftshift(img)) * bp; F = F[np.ix_(sl, sl)]; return np.fft.fftshift(np.real(np.fft.ifft2(F))).astype(np.float32)
def crop_real(img): F = np.fft.fft2(np.fft.ifftshift(img))[np.ix_(sl, sl)]; return np.fft.fftshift(np.real(np.fft.ifft2(F))).astype(np.float32) * (c / n) ** 2
PE, PO, PF = [np.array([prep(im) for im in X]) for X in (even, odd, full)]; BA = np.array([prep(im) for im in binavg])
members = {b: np.where(parent == b + 1)[0] for b in range(nb)}
# projected region masks per bin (soft), bin-centred contrasts, principal contrast per bin from the FULL averages, coefficients from the halves
PM = np.array([np.clip(crop_real(project(V, rot(E1[b], E2[b]))), 0, None) for b in range(nb)]); PM = np.array([m / (m.max() + 1e-12) for m in PM]); PM = np.clip(PM * 2., 0, 1)
v_b, coef, ok = [None] * nb, np.zeros(N), np.zeros(nb, bool); ev = np.zeros(nb); v2_b, coef2 = [None] * nb, np.zeros(N)
for b in range(nb):
    mem = members[b]
    if len(mem) < 3: continue
    w = pop[mem] / pop[mem].sum(); D = PF[mem] - (w[:, None, None] * PF[mem]).sum(0); X = (D * PM[b]).reshape(len(mem), -1) * np.sqrt(w)[:, None]
    U, s, Vt = np.linalg.svd(X, full_matrices=False); v_b[b] = Vt[0].reshape(c, c); ev[b] = s[0] ** 2 / max((s ** 2).sum(), 1e-30); ok[b] = True
    DE = (PE[mem] - (w[:, None, None] * PE[mem]).sum(0)).reshape(len(mem), -1); DO = (PO[mem] - (w[:, None, None] * PO[mem]).sum(0)).reshape(len(mem), -1)
    coef[mem] = 0.5 * (DE @ Vt[0] + DO @ Vt[0])   # cross-half-averaged coefficient of each sub-class along its bin's contrast
    if len(mem) >= 4 and a.n2 > 1: v2_b[b] = Vt[1].reshape(c, c); coef2[mem] = 0.5 * (DE @ Vt[1] + DO @ Vt[1])
print('%d bins with a principal contrast; explained within-bin variance in the region: median %.2f' % (ok.sum(), np.median(ev[ok])), flush=True)
# Z2 sign synchronisation between neighbouring bins
ang = np.degrees(np.arccos(np.clip(C @ C.T, -1, 1))); angs = np.arange(-180, 180, 4.); Wm = np.zeros((nb, nb))
def corr(u, v): u, v = u - u.mean(), v - v.mean(); return float((u * v).sum() / (np.sqrt((u * u).sum() * (v * v).sum()) + 1e-30))
for i in range(nb):
    if not ok[i]: continue
    for j in np.argsort(ang[i])[1:a.nn + 1]:
        if not ok[j] or ang[i, j] > a.maxang or Wm[i, j] != 0: continue
        psi = max(((corr(BA[i], ndrotate(BA[j], p, reshape=False, order=1)), p) for p in angs))[1]
        s_ij = corr(v_b[i] * PM[i], ndrotate(v_b[j] * PM[j], psi, reshape=False, order=1)); Wm[i, j] = Wm[j, i] = s_ij
vals, vecs = np.linalg.eigh(Wm); sign = np.sign(vecs[:, -1]); sign[sign == 0] = 1
agree = np.mean([np.sign(Wm[i, j]) == sign[i] * sign[j] for i in range(nb) for j in range(nb) if Wm[i, j] != 0]) if (Wm != 0).any() else float('nan')
print('sign synchronisation: top eigenvalue %.1f (next %.1f), %d edges, %.0f%% of edges consistent with the synchronised signs' % (vals[-1], vals[-2], int((Wm != 0).sum() / 2), 100 * agree), flush=True)
for b in range(nb):
    if ok[b]: coef[members[b]] *= sign[b]
if a.n2 > 1:   # second contrast: its own sign sync (same neighbour pairs, same in-plane calibration)
    W2 = np.zeros((nb, nb))
    for i in range(nb):
        for j in range(nb):
            if Wm[i, j] != 0 and v2_b[i] is not None and v2_b[j] is not None:
                psi = max(((corr(BA[i], ndrotate(BA[j], p, reshape=False, order=1)), p) for p in angs))[1]; W2[i, j] = corr(v2_b[i] * PM[i], ndrotate(v2_b[j] * PM[j], psi, reshape=False, order=1))
    vals2, vecs2 = np.linalg.eigh(W2); sign2 = np.sign(vecs2[:, -1]); sign2[sign2 == 0] = 1
    for b in range(nb):
        if ok[b]: coef2[members[b]] *= sign2[b]
    print('second contrast: sign sync top eigenvalue %.1f (next %.1f)' % (vals2[-1], vals2[-2]), flush=True)
# 3D deviation mode: back-projection of the signed contrasts (Fourier slices, nearest voxel, validated convention R^T, e1 +)
fc = np.fft.fftfreq(c) * c / n; fy2, fx2 = np.meshgrid(fc, fc, indexing='ij'); kk = np.stack([fx2, fy2, np.zeros_like(fx2)], -1).reshape(-1, 3) * n
num = np.zeros(c ** 3, np.complex64); den = np.zeros(c ** 3, np.float32)
for b in range(nb):
    if not ok[b]: continue
    q = np.round(kk @ rot(E1[b], E2[b]).T).astype(int); good = (np.abs(q) <= c // 2 - 1).all(1); q = q % c; idx = q[:, 2] * c * c + q[:, 1] * c + q[:, 0]
    F = np.fft.fft2(np.fft.ifftshift(sign[b] * v_b[b] * PM[b])).ravel() * good; np.add.at(num, idx, F); np.add.at(den, idx, good.astype(np.float32))
mode = np.fft.fftshift(np.real(np.fft.ifftn((num / np.maximum(den, 1)).reshape(c, c, c)))).astype(np.float32); write_mrc(a.out + '/mode3d.mrc', mode, a.smpd * n / c)
# particle coordinates from the sub-class coefficients -> lattice init (n1 bins along the coordinate, n2 = 1)
L = np.loadtxt(a.clsexp + '/cls_expansion_labels.txt', comments='#').astype(int); npart = max(int(L[:, 0].max()), a.npart); row_of = {int(g): j for j, g in enumerate(M['gid'])}; rows = np.array([row_of.get(int(g), 0) for g in L[:, 2]]); z = np.zeros(npart); z[L[:, 0] - 1] = coef[rows]
zs = z[L[:, 0] - 1]; r = np.argsort(np.argsort(zs)); g1 = (r * a.n1 // len(zs)) + 1
z2 = np.zeros(npart); z2[L[:, 0] - 1] = coef2[rows]; r2 = np.argsort(np.argsort(z2[L[:, 0] - 1])); g2 = (r2 * a.n2 // len(zs)) + 1 if a.n2 > 1 else np.ones(len(zs), int)
with open(a.out + '/init_positions.txt', 'w') as f:
    for p, gg, hh in zip(L[:, 0], g1, g2): f.write('%d %d %d\n' % (p, gg, hh))
np.savetxt(a.out + '/particle_coord.txt', z, fmt='%.5f'); np.savetxt(a.out + '/particle_bin.txt', np.r_[np.zeros(0, int), (lambda full_: full_)(np.zeros(npart, int))], fmt='%d')
lab = np.zeros(npart, int); lab[L[:, 0] - 1] = g1; np.savetxt(a.out + '/particle_bin.txt', lab, fmt='%d')
print('seed written: %s/init_positions.txt (%d particles, %d x %d lattice along the synchronised contrast coordinates)' % (a.out, len(zs), a.n1, a.n2))
