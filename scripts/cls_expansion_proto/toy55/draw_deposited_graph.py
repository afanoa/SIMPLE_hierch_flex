#!/usr/bin/env python3
"""(1) The deposited conformers as a state graph: deviation from their image-weighted mean, MDS layout of the deviation
distances in a band, 2-NN edges, labels = conformer, images (paper), closed flag. (2) Joint layout of OUR states (cells
of a klem run on the coarse grid) with the deposited conformers (Fourier-cropped to the same grid): where our path
lies relative to theirs. usage: draw_deposited_graph.py DEPOSITED_DIR CONSENSUS CSV CELLS_DIR LATENT OUT_PREFIX"""
import sys, glob, re, csv, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
dep, consf, csvf, cells, latent, out = sys.argv[1:7]; SMPD = 1.464; CLOSED = {18, 19, 30, 50, 51, 52, 53}
def read_mrc(p):
    h = np.fromfile(p, dtype=np.int32, count=256); nx, ny, nz = h[0], h[1], h[2]; return np.fromfile(p, dtype=np.float32, offset=1024 + h[23]).reshape(nz, ny, nx).astype(float)
rows = list(csv.DictReader(open(csvf))); ids = [int(r['conformer']) for r in rows]; nimg = {int(r['conformer']): int(r['images']) for r in rows}
C = {i: read_mrc(glob.glob('%s/conf%02d_*.mrc' % (dep, i))[0]) for i in ids}; cons = read_mrc(consf); n = cons.shape[0]
zz, yy, xx = np.indices((n, n, n)); r = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2); msk = 1. / (1. + np.exp((r - 100 / SMPD) / 2.))
F = lambda v: np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v * msk)))
def corrmat(X):
    Xn = X / np.maximum(np.sqrt(np.sum(np.abs(X) ** 2, 1))[:, None], 1e-30); return np.real(Xn @ np.conj(Xn).T)
def mds(S):
    K = len(S); D2 = np.maximum(2 * (1 - S), 0); np.fill_diagonal(D2, 0); J = np.eye(K) - 1. / K; B = -0.5 * J @ D2 @ J; w, v = np.linalg.eigh(B); pos = v[:, -2:][:, ::-1] * np.sqrt(np.maximum(w[-2:][::-1], 0)); return (pos - pos.mean(0)) / pos.std(0)
def draw(ax, pos, thumbs, frames, labels, S, knn, zoom):
    K = len(pos); W = np.zeros((K, K))
    for i in range(K):
        for j in np.argsort(-(S[i] - 9 * np.eye(K)[i]))[:knn]: W[i, j] = W[j, i] = max(S[i, j], 0)
    for i in range(K):
        for j in range(i + 1, K):
            if W[i, j] > 0: ax.plot([pos[i, 0], pos[j, 0]], [pos[i, 1], pos[j, 1]], color='k', alpha=0.3 + 0.6 * max(W[i, j], 0), lw=0.6 + 2.5 * max(W[i, j], 0), zorder=1)
    for k in range(K):
        img = thumbs[k]; vk = np.abs(img).max(); ob = OffsetImage(img, cmap='RdBu_r', zoom=zoom, norm=matplotlib.colors.Normalize(-vk, vk)); ab = AnnotationBbox(ob, pos[k], frameon=True, bboxprops=dict(edgecolor=frames[k], lw=4), zorder=3); ax.add_artist(ab)
        ax.text(pos[k, 0], pos[k, 1] - 0.17, labels[k], ha='center', va='top', fontsize=9, zorder=4, bbox=dict(facecolor='white', alpha=0.75, lw=0, pad=1))
    ax.set_xlim(pos[:, 0].min() - 0.5, pos[:, 0].max() + 0.5); ax.set_ylim(pos[:, 1].min() - 0.5, pos[:, 1].max() + 0.5); ax.set_xticks([]); ax.set_yticks([])
# ---- 1. deposited alone, full box, band 60-12 A
band = (r >= n * SMPD / 60.) & (r <= n * SMPD / 12.); Dd = {i: C[i] - cons for i in ids}; Xd = np.array([F(Dd[i])[band] for i in ids]); Sd = corrmat(Xd); pos = mds(Sd)
print('deposited conformers, deviation correlations at 60-12 A:'); print('      ' + ' '.join('%5d' % i for i in ids))
for a_, i in enumerate(ids): print('%5d ' % i + ' '.join('%+.2f' % Sd[a_, b_] for b_ in range(len(ids))))
fig, ax = plt.subplots(figsize=(12, 10)); draw(ax, pos, [Dd[i].sum(1) for i in ids], ['tab:red' if i in CLOSED else 'goldenrod' for i in ids], ['conf %d%s\n%dk images' % (i, ' (closed)' if i in CLOSED else '', nimg[i] // 1000) for i in ids], Sd, 2, 0.55)
ax.set_title('the ten deposited conformers as a state graph: deviation from their image-weighted mean, layout by deviation distance (60-12 A), 2-NN edges; red frame = closed in the paper', fontsize=10); plt.tight_layout(); plt.savefig(out + '_deposited.png', dpi=85); print('figure ->', out + '_deposited.png')
# ---- 2. joint with our states on the coarse grid, band 60-20 A
fs = sorted(glob.glob(cells + '/cell_*.mrc')) + sorted(glob.glob(cells + '/state_*.mrc')); V = np.array([read_mrc(f) for f in fs]); K = len(V); nc = V.shape[1]; Dv = V - V.mean(0)
def crop(v):   # Fourier-crop a full-box map to the coarse grid
    Fv = F(v); c = n // 2; h = nc // 2; sub = Fv[c - h:c - h + nc, c - h:c - h + nc, c - h:c - h + nc]; return np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(sub))))
Dc = {i: crop(C[i] - cons) for i in ids}; smpd_c = n * SMPD / nc; zc, yc, xc = np.indices((nc, nc, nc)); rc = np.hypot(np.hypot(xc - nc / 2, yc - nc / 2), zc - nc / 2); bandc = (rc >= nc * smpd_c / 60.) & (rc <= nc * smpd_c / 20.)
Fc = lambda v: np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v)))
X = np.array([Fc(d)[bandc] for d in Dv] + [Fc(Dc[i])[bandc] for i in ids]); S = corrmat(X); pos = mds(S)
L = np.loadtxt(latent, comments='#'); am = L[:, 3].astype(int) - 1; gm = L[:, 5]; cnt = np.bincount(am[(gm > 0.5) & (am < K)], minlength=K); z1 = np.array([L[(gm > 0.5) & (am == k), 1].mean() if np.any((gm > 0.5) & (am == k)) else 0 for k in range(K)])
cm = plt.get_cmap('viridis'); frames = [cm((z - z1.min()) / max(z1.max() - z1.min(), 1e-9)) for z in z1] + ['tab:red' if i in CLOSED else 'goldenrod' for i in ids]
labels = ['%dk' % (cnt[k] // 1000) for k in range(K)] + ['DEPOSITED %d%s' % (i, ' closed' if i in CLOSED else '') for i in ids]
best = [ids[int(np.argmax(S[k, K:]))] for k in range(K)]; print('our states: best deposited match per state (count): %s' % {i: best.count(i) for i in ids}); print('mean correlation of our states with their best deposited match: %.2f; deposited-deposited nearest: %.2f' % (np.mean([S[k, K + ids.index(best[k])] for k in range(K)]), np.median(np.sort(S[K:, K:] - 9 * np.eye(len(ids)), 1)[:, -1])))
fig, ax = plt.subplots(figsize=(16, 14)); draw(ax, pos, [d.sum(1) for d in Dv] + [Dc[i].sum(1) for i in ids], frames, labels, S, 3, 1.3)
ax.set_title('our 49 states (frames = z1, label = particles) embedded together with the ten deposited conformers (gold / red = closed), by deviation distance at 60-20 A on the coarse grid; 3-NN edges', fontsize=11); plt.tight_layout(); plt.savefig(out + '_joint.png', dpi=80); print('figure ->', out + '_joint.png')
