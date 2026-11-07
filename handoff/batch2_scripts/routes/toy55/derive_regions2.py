#!/usr/bin/env python3
"""Data-derived focal regions by MOTION PATTERN rather than by spatial connectivity: every moving voxel of a klem board
(cross-half reproducible variability above --thr) carries a profile over the surviving cells (its half-averaged
deviation per cell); voxels are k-means clustered on the normalised profiles (things that move together), clusters
with correlated centroids merged, each cluster split into its spatially connected parts, parts below --minvox dropped,
the rest dilated, softened and written as full-box masks (largest variability mass first). --ref MASK: overlap readout only.
usage: derive_regions2.py BOARD OUT [--thr 0.5 --k 5 --merge 0.8 --minvox 40 --dilate 2 --maxreg 6 --ref MASK]"""
import sys, os, argparse, numpy as np
from scipy.ndimage import gaussian_filter, label, binary_dilation, zoom
from scipy.cluster.vq import kmeans2
ap = argparse.ArgumentParser(); ap.add_argument('board'); ap.add_argument('out'); ap.add_argument('--thr', type=float, default=0.5); ap.add_argument('--k', type=int, default=5); ap.add_argument('--merge', type=float, default=0.8)
ap.add_argument('--minvox', type=int, default=40); ap.add_argument('--dilate', type=int, default=2); ap.add_argument('--maxreg', type=int, default=6); ap.add_argument('--box', type=int, default=160); ap.add_argument('--ref', default=''); ap.add_argument('--minrep', type=float, default=0.5); ap.add_argument('--within', default='', help='full-box mask: only moving voxels inside it are segmented (the parent node region, inherited)'); a = ap.parse_args()
SMPD = 1.464; os.makedirs(a.out, exist_ok=True); rng = np.random.default_rng(0)
def read_mrc(p):
    h = np.fromfile(p, dtype=np.int32, count=256); nx, ny, nz = h[0], h[1], h[2]; return np.fromfile(p, dtype=np.float32, offset=1024 + h[23]).reshape(nz, ny, nx).astype(float)
def write_mrc(fn, v, smpd=SMPD):
    v = np.ascontiguousarray(v, np.float32); nz, ny, nx = v.shape; h = np.zeros(256, np.int32); hf = h.view(np.float32)
    h[0:3] = (nx, ny, nz); h[3] = 2; h[7:10] = (nx, ny, nz); hf[10:13] = (nx * smpd, ny * smpd, nz * smpd); hf[13:16] = 90.; h[16:19] = (1, 2, 3)
    hf[19:22] = (v.min(), v.max(), v.mean()); h[52] = int.from_bytes(b'MAP ', 'little'); h[53] = int.from_bytes(bytes([68, 65, 0, 0]), 'little'); hf[54] = v.std(); open(fn, 'wb').write(h.tobytes() + v.tobytes())
H = np.load(a.board + '/klem_halves.npz'); n3 = int(H['n3']); band = H['band']; D0, D1 = H['D0'], H['D1']; st = np.loadtxt(a.board + '/klem_states.txt', comments='#'); keep = np.where(st[:len(D0), 2] >= a.minrep)[0]
ifft = lambda X: np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(X))))
def full(d):
    X = np.zeros(n3 ** 3, complex); X[band] = d; return ifft(X.reshape(n3, n3, n3))
dev0 = np.array([full(D0[k]) for k in keep]); dev1 = np.array([full(D1[k]) for k in keep]); v = gaussian_filter(np.maximum(np.sum(dev0 * dev1, 0), 0), 1.5); fm = np.clip(v / max(np.percentile(v, 99), 1e-30), 0, 1) ** 0.5
mov = fm > a.thr
if a.within:
    wm = read_mrc(a.within); nb = wm.shape[0]; Fw = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(wm))); c = nb // 2; h = n3 // 2
    wmc = np.clip(np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(Fw[c - h:c - h + n3, c - h:c - h + n3, c - h:c - h + n3])))) * (n3 / nb) ** 3, 0, 1)
    mov = mov & (wmc > 0.5); print('restricted to the inherited region: %d moving voxels inside it' % mov.sum())
idx = np.where(mov.ravel())[0]; prof = (0.5 * (dev0 + dev1)).reshape(len(keep), -1)[:, idx].T; prof = prof - prof.mean(1, keepdims=True); prof /= np.maximum(np.linalg.norm(prof, axis=1, keepdims=True), 1e-30)
best = None
for seed in range(12):
    c, l = kmeans2(prof, a.k, minit='++', seed=seed); inertia = np.sum((prof - c[l]) ** 2)
    if best is None or inertia < best[0]: best = (inertia, c, l)
_, cent, lab = best; C = cent / np.maximum(np.linalg.norm(cent, axis=1, keepdims=True), 1e-30); S = C @ C.T
groups = list(range(a.k)); parent = list(range(a.k))
for i in range(a.k):
    for j in range(i + 1, a.k):
        if S[i, j] > a.merge: parent[j] = parent[i]
lab = np.array([parent[x] for x in lab]); smpd_c = a.box * SMPD / n3; zz, yy, xx = np.indices((n3, n3, n3)); rc = np.hypot(np.hypot(xx - n3 / 2, yy - n3 / 2), zz - n3 / 2) * smpd_c
print('board %s: %d cells; %d moving coarse voxels (fm > %.2f); k-means k=%d on motion profiles -> %d groups after merging centroids above %.2f (centroid corr matrix off-diagonal max %.2f)' % (a.board.rstrip('/').split('/')[-1], len(keep), len(idx), a.thr, a.k, len(set(lab)), a.merge, (S - np.eye(a.k)).max()))
parts = []
for g in sorted(set(lab)):
    m = np.zeros(n3 ** 3, bool); m[idx[lab == g]] = True; m = m.reshape(n3, n3, n3); L, nl = label(m)
    for i in range(1, nl + 1):
        comp = L == i
        if comp.sum() >= a.minvox: parts.append((v[comp].sum(), g, comp))
parts.sort(key=lambda t: -t[0]); parts = parts[:a.maxreg]; ref = read_mrc(a.ref) > 0.5 if a.ref else None; union = np.zeros((a.box,) * 3)
for j, (mass, g, comp) in enumerate(parts):
    soft = gaussian_filter(binary_dilation(comp, iterations=a.dilate).astype(float), 1.0); soft = np.clip(soft / max(soft.max(), 1e-30), 0, 1); fullm = np.clip(zoom(soft, a.box / n3, order=1), 0, 1)
    fullm = fullm[:a.box, :a.box, :a.box] if fullm.shape[0] >= a.box else np.pad(fullm, [(0, a.box - fullm.shape[0])] * 3); write_mrc('%s/region_%d.mrc' % (a.out, j + 1), fullm); union = np.maximum(union, fullm)
    extra = ''
    if ref is not None: rm = fullm > 0.5; extra = ' | reference mask: %.0f %% of it covered, %.0f %% of the region inside it' % (100. * (rm & ref).sum() / max(ref.sum(), 1), 100. * (rm & ref).sum() / max(rm.sum(), 1))
    print('  region %d (motion group %d): %d coarse voxels, %.1f %% of the variability mass, centroid radius %.0f A, max radius %.0f A%s' % (j + 1, g, comp.sum(), 100. * mass / v[mov].sum(), np.hypot(np.hypot(xx[comp].mean() - n3 / 2, yy[comp].mean() - n3 / 2), zz[comp].mean() - n3 / 2) * smpd_c, rc[comp].max(), extra))
write_mrc(a.out + '/all_moving.mrc', union); print('masks -> %s/region_*.mrc (%d)' % (a.out, len(parts)))
