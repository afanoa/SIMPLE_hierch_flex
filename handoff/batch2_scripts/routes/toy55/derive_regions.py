#!/usr/bin/env python3
"""Data-derived focal regions of a klem board: the cross-half reproducible variability of the surviving cells
(sum_k dev0_k x dev1_k in real space on the coarse grid, klem's own --focus map), thresholded and split into connected
components; each component dilated, softened and written as a full-box mask (region_k.mrc), largest variability mass
first. The deposited P12/P10.1 mask may be given as --ref for a reference overlap readout only.
usage: derive_regions.py BOARD OUT [--thr 0.5 --minvox 40 --dilate 2 --maxreg 6 --box 160 --ref MASK]"""
import sys, os, argparse, numpy as np
from scipy.ndimage import gaussian_filter, label, binary_dilation, binary_opening, zoom
ap = argparse.ArgumentParser(); ap.add_argument('board'); ap.add_argument('out'); ap.add_argument('--thr', type=float, default=0.5); ap.add_argument('--minvox', type=int, default=40)
ap.add_argument('--dilate', type=int, default=2); ap.add_argument('--maxreg', type=int, default=6); ap.add_argument('--box', type=int, default=160); ap.add_argument('--ref', default=''); ap.add_argument('--minrep', type=float, default=0.5); a = ap.parse_args()
SMPD = 1.464; os.makedirs(a.out, exist_ok=True)
def read_mrc(p):
    h = np.fromfile(p, dtype=np.int32, count=256); nx, ny, nz = h[0], h[1], h[2]; return np.fromfile(p, dtype=np.float32, offset=1024 + h[23]).reshape(nz, ny, nx).astype(float)
def write_mrc(fn, v, smpd=SMPD):
    v = np.ascontiguousarray(v, np.float32); nz, ny, nx = v.shape; h = np.zeros(256, np.int32); hf = h.view(np.float32)
    h[0:3] = (nx, ny, nz); h[3] = 2; h[7:10] = (nx, ny, nz); hf[10:13] = (nx * smpd, ny * smpd, nz * smpd); hf[13:16] = 90.; h[16:19] = (1, 2, 3)
    hf[19:22] = (v.min(), v.max(), v.mean()); h[52] = int.from_bytes(b'MAP ', 'little'); h[53] = int.from_bytes(bytes([68, 65, 0, 0]), 'little'); hf[54] = v.std(); open(fn, 'wb').write(h.tobytes() + v.tobytes())
H = np.load(a.board + '/klem_halves.npz'); n3 = int(H['n3']); band = H['band']; D0, D1 = H['D0'], H['D1']; st = np.loadtxt(a.board + '/klem_states.txt', comments='#'); crep = st[:len(D0), 2]; keep = crep >= a.minrep
ifft = lambda X: np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(X))))
def full(d):
    X = np.zeros(n3 ** 3, complex); X[band] = d; return ifft(X.reshape(n3, n3, n3))
v = np.zeros((n3, n3, n3))
for k in np.where(keep)[0]: v += full(D0[k]) * full(D1[k])
v = gaussian_filter(np.maximum(v, 0), 1.5); fm = np.clip(v / max(np.percentile(v, 99), 1e-30), 0, 1) ** 0.5
smpd_c = a.box * SMPD / n3; zz, yy, xx = np.indices((n3, n3, n3)); rc = np.hypot(np.hypot(xx - n3 / 2, yy - n3 / 2), zz - n3 / 2) * smpd_c
m = binary_opening(fm > a.thr, iterations=1); lab, nl = label(m); sizes = np.array([(lab == i).sum() for i in range(1, nl + 1)]); mass = np.array([v[lab == i].sum() for i in range(1, nl + 1)])
order = [i for i in np.argsort(-mass) if sizes[i] >= a.minvox][:a.maxreg]
print('board %s: %d of %d cells survive; coarse grid %d^3 (%.1f A/px); variability above %.2f: %d voxels in %d components, %d kept (>= %d voxels)' % (a.board.rstrip('/').split('/')[-1], keep.sum(), len(keep), n3, smpd_c, a.thr, int(m.sum()), nl, len(order), a.minvox))
ref = read_mrc(a.ref) > 0.5 if a.ref else None
union = np.zeros((a.box,) * 3)
for j, i in enumerate(order):
    comp = binary_dilation(lab == i + 1, iterations=a.dilate); soft = gaussian_filter(comp.astype(float), 1.0); soft = np.clip(soft / max(soft.max(), 1e-30), 0, 1)
    fullm = np.clip(zoom(soft, a.box / n3, order=1), 0, 1); fullm = fullm[:a.box, :a.box, :a.box] if fullm.shape[0] >= a.box else np.pad(fullm, [(0, a.box - fullm.shape[0])] * 3)
    write_mrc('%s/region_%d.mrc' % (a.out, j + 1), fullm); union = np.maximum(union, fullm)
    cz, cy, cx = [float(np.mean(q[lab == i + 1])) for q in (zz, yy, xx)]; rad = np.hypot(np.hypot(cx - n3 / 2, cy - n3 / 2), cz - n3 / 2) * smpd_c
    extra = ''
    if ref is not None: rm = fullm > 0.5; extra = ' | overlap with the reference mask: %.0f %% of it covered, %.0f %% of the region inside it' % (100. * (rm & ref).sum() / max(ref.sum(), 1), 100. * (rm & ref).sum() / max(rm.sum(), 1))
    print('  region %d: %d coarse voxels (%.1f %% of the variability mass), centroid at %.0f A radius, max radius %.0f A%s' % (j + 1, sizes[i], 100. * mass[i] / mass.sum(), rad, rc[lab == i + 1].max(), extra))
write_mrc(a.out + '/all_moving.mrc', union); write_mrc(a.out + '/variability_coarse.mrc', fm, smpd_c); print('masks -> %s/region_*.mrc (%d), all_moving.mrc, variability_coarse.mrc' % (a.out, len(order)))
