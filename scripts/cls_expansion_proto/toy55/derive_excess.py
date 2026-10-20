#!/usr/bin/env python3
"""Candidate ADDED-DENSITY sites of a klem board (truth-free): the cross-half reproducible POSITIVE deviation of the
surviving cells, sum_k max(dev0_k, 0) x max(dev1_k, 0) in real space on the coarse grid, optionally down-weighted where
the consensus already has density (added density sits where the mean map is weak) and with the main moving region
excluded (--exclude MASK, e.g. the node's region_1); thresholded, connected components ranked by mass, each written as a
full-box soft mask (excess_k.mrc). --ref MASK gives an overlap readout only (the truth site on a toy).
usage: derive_excess.py BOARD OUT [--thr 0.5 --minvox 20 --dilate 2 --maxreg 4 --box 160 --exclude MASK --ref MASK --envelope 1]"""
import sys, os, argparse, numpy as np
from scipy.ndimage import gaussian_filter, label, binary_dilation, binary_opening, zoom
ap = argparse.ArgumentParser(); ap.add_argument('board'); ap.add_argument('out'); ap.add_argument('--thr', type=float, default=0.5); ap.add_argument('--minvox', type=int, default=20)
ap.add_argument('--dilate', type=int, default=2); ap.add_argument('--maxreg', type=int, default=4); ap.add_argument('--box', type=int, default=160); ap.add_argument('--ref', default=''); ap.add_argument('--minrep', type=float, default=0.5)
ap.add_argument('--exclude', default='', help='full-box mask of the main moving region: voxels inside it are not candidates'); ap.add_argument('--envelope', type=float, default=1., help='down-weight by (1 - envelope x normalised consensus density): added density sits where the mean map is weak'); a = ap.parse_args()
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
def crop_mask(fn):   # full-box mask -> coarse grid
    mk = read_mrc(fn); nb = mk.shape[0]; Fm = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(mk))); c = nb // 2; h = n3 // 2
    return np.clip(ifft(Fm[c - h:c - h + n3, c - h:c - h + n3, c - h:c - h + n3]) * (n3 / nb) ** 3, 0, 1)
v = np.zeros((n3, n3, n3)); cons = np.zeros((n3, n3, n3)); nk = 0
for k in np.where(keep)[0]: d0, d1 = full(D0[k]), full(D1[k]); v += np.maximum(d0, 0) * np.maximum(d1, 0); nk += 1
v = gaussian_filter(v, 1.0)
if a.envelope > 0 and os.path.exists(a.board + '/consensus_coarse.npy'):
    cons = np.load(a.board + '/consensus_coarse.npy'); cons = np.clip(cons / max(np.percentile(cons, 99.5), 1e-30), 0, 1); v = v * (1 - a.envelope * cons)
if a.exclude: v = v * (1 - crop_mask(a.exclude))
fm = np.clip(v / max(np.percentile(v, 99.5), 1e-30), 0, 1)
smpd_c = a.box * SMPD / n3; zz, yy, xx = np.indices((n3, n3, n3)); rc = np.hypot(np.hypot(xx - n3 / 2, yy - n3 / 2), zz - n3 / 2) * smpd_c
m = binary_opening(fm > a.thr, iterations=1); lab, nl = label(m); sizes = np.array([(lab == i).sum() for i in range(1, nl + 1)]); mass = np.array([v[lab == i].sum() for i in range(1, nl + 1)])
order = [i for i in np.argsort(-mass) if sizes[i] >= a.minvox][:a.maxreg]
print('board %s: %d of %d cells survive; coarse grid %d^3 (%.1f A/px); reproducible positive excess above %.2f: %d voxels in %d components, %d kept (>= %d voxels)%s' % (a.board.rstrip('/').split('/')[-1], keep.sum(), len(keep), n3, smpd_c, a.thr, int(m.sum()), nl, len(order), a.minvox, ', main moving region excluded' if a.exclude else ''))
ref = read_mrc(a.ref) > 0.5 if a.ref else None
for j, i in enumerate(order):
    comp = binary_dilation(lab == i + 1, iterations=a.dilate); soft = gaussian_filter(comp.astype(float), 1.0); soft = np.clip(soft / max(soft.max(), 1e-30), 0, 1)
    fullm = np.clip(zoom(soft, a.box / n3, order=1), 0, 1); fullm = fullm[:a.box, :a.box, :a.box] if fullm.shape[0] >= a.box else np.pad(fullm, [(0, a.box - fullm.shape[0])] * 3)
    write_mrc('%s/excess_%d.mrc' % (a.out, j + 1), fullm)
    cz, cy, cx = [float(np.mean(q[lab == i + 1])) for q in (zz, yy, xx)]; rad = np.hypot(np.hypot(cx - n3 / 2, cy - n3 / 2), cz - n3 / 2) * smpd_c; extra = ''
    if ref is not None: rm = fullm > 0.5; extra = ' | overlap with the reference site: %.0f %% of it covered, %.0f %% of the candidate inside it' % (100. * (rm & ref).sum() / max(ref.sum(), 1), 100. * (rm & ref).sum() / max(rm.sum(), 1))
    print('  candidate %d: %d coarse voxels (%.1f %% of the excess mass), centroid at %.0f A radius%s' % (j + 1, sizes[i], 100. * mass[i] / max(mass.sum(), 1e-30), rad, extra))
write_mrc(a.out + '/excess_coarse.mrc', fm, smpd_c); print('masks -> %s/excess_*.mrc (%d), excess_coarse.mrc' % (a.out, len(order)))
