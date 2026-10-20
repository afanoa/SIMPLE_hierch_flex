#!/usr/bin/env python3
"""Draw the state graph: nodes = states (their deviation-from-mean maps as thumbnails), edges = similarity of the
deviations (kNN graph), layout = spectral embedding refined by a few force-directed steps. Node colour: --gt gives the
dominant true conformer per state (toy); --closedness colours by the mean sub-class closedness of the state's particles
(real data, needs --meta); otherwise the mean z1. usage: draw_state_graph.py CELLS_DIR LATENT OUT.png [--knn 4] [--gt gt.txt] [--meta meta.npz --closedness table.txt]"""
import sys, glob, re, argparse, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
ap = argparse.ArgumentParser(); ap.add_argument('cells'); ap.add_argument('latent'); ap.add_argument('out'); ap.add_argument('--knn', type=int, default=4); ap.add_argument('--gt', default=''); ap.add_argument('--meta', default=''); ap.add_argument('--closedness', default=''); ap.add_argument('--title', default='')
ap.add_argument('--lo', type=float, default=60.); ap.add_argument('--hi', type=float, default=20.); a = ap.parse_args()
def read_mrc(p):
    h = np.fromfile(p, dtype=np.int32, count=256); nx, ny, nz = h[0], h[1], h[2]; return np.fromfile(p, dtype=np.float32, offset=1024 + h[23]).reshape(nz, ny, nx).astype(float), h
fs = sorted(glob.glob(a.cells + '/cell_*.mrc')) + sorted(glob.glob(a.cells + '/state_*.mrc')); V = []; hdr = None
for f in fs: v, hdr = read_mrc(f); V.append(v)
V = np.array(V); K = len(V); n = V.shape[1]; smpd = float(np.frombuffer(hdr[10:13].tobytes(), np.float32)[0]) / n   # coarse-grid pixel size
mean = V.mean(0); Dv = V - mean; zz, yy, xx = np.indices((n, n, n)); r = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2)
F = np.array([np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(d))) for d in Dv]); band = (r >= n * smpd / a.lo) & (r <= n * smpd / a.hi) if n * smpd / a.hi > 2 else r >= 1
X = F[:, band]; X = X / np.maximum(np.sqrt(np.sum(np.abs(X) ** 2, 1))[:, None], 1e-30); S = np.real(X @ np.conj(X).T); np.fill_diagonal(S, 0)
# kNN graph (symmetrised), weights = similarity
W = np.zeros((K, K))
for i in range(K):
    for j in np.argsort(-S[i])[:a.knn]: W[i, j] = max(S[i, j], 0); W[j, i] = max(W[j, i], W[i, j])
# layout: classical MDS of the deviation distances d^2 = 2(1 - S) (keeps the continuum), 2 dimensions
D2 = np.maximum(2 * (1 - (np.real(X @ np.conj(X).T))), 0); J = np.eye(K) - 1. / K; B = -0.5 * J @ D2 @ J; w_, v_ = np.linalg.eigh(B); pos = v_[:, -2:][:, ::-1] * np.sqrt(np.maximum(w_[-2:][::-1], 0)); pos = (pos - pos.mean(0)) / pos.std(0)
L = np.loadtxt(a.latent, comments='#'); am = L[:, 3].astype(int) - 1; gm = L[:, 5]; on = gm > 0.5; cnt = np.bincount(am[on & (am < K)], minlength=K)
col = None; lab = None
if a.gt:
    gt = {}
    for l in open(a.gt):
        if l.strip() and not l.startswith('#'): w = l.split(); gt[int(w[0])] = int(w[1])
    conf = np.array([gt.get(int(p), 0) for p in L[:, 0]]); cs = sorted(set(conf[conf > 0])); cmap = plt.get_cmap('tab10'); col = []; lab = []
    for k in range(K):
        m = on & (am == k) & (conf > 0)
        if m.sum() == 0: col.append((0.8, 0.8, 0.8, 1)); lab.append(''); continue
        bc = np.bincount(conf[m]); c = bc.argmax(); col.append(cmap(cs.index(c) % 10)); lab.append('%d\n%.0f%%' % (c, 100 * bc.max() / m.sum()))
elif a.closedness and a.meta:
    M = np.load(a.meta); flags = M['flags'].astype(int); ref = np.where(flags == 1)[0]; posr = {int(i): j for j, i in enumerate(ref)}; C = np.loadtxt(a.closedness); cl = dict(zip(C[:, 0].astype(int), C[:, 4]))
    rows = np.array([posr.get(int(p) - 1, -1) for p in L[:, 0]]); sub = np.where(rows >= 0, M['lab692'][np.maximum(rows, 0)], 0); cv = np.array([cl.get(int(s_), np.nan) for s_ in sub])
    val = np.array([np.nanmean(cv[on & (am == k)]) if np.any(on & (am == k)) else np.nan for k in range(K)]); cmap = plt.get_cmap('coolwarm'); col = [cmap(np.clip((v + 0.3) / 0.9, 0, 1)) if np.isfinite(v) else (0.8, 0.8, 0.8, 1) for v in val]; lab = ['%+.2f' % v if np.isfinite(v) else '' for v in val]
else:
    z1 = L[:, 1]; val = np.array([z1[on & (am == k)].mean() if np.any(on & (am == k)) else np.nan for k in range(K)]); cmap = plt.get_cmap('viridis'); col = [cmap((v - np.nanmin(val)) / max(np.nanmax(val) - np.nanmin(val), 1e-9)) if np.isfinite(v) else (0.8, 0.8, 0.8, 1) for v in val]; lab = ['' for _ in range(K)]
fig, ax = plt.subplots(figsize=(16, 14)); ax.set_xticks([]); ax.set_yticks([])
for i in range(K):
    for j in range(i + 1, K):
        if W[i, j] > 0: ax.plot([pos[i, 0], pos[j, 0]], [pos[i, 1], pos[j, 1]], color='k', alpha=0.25 + 0.6 * max(W[i, j] - 0.8, 0) / 0.2, lw=0.6 + 2.5 * max(W[i, j] - 0.8, 0) / 0.2, zorder=1)
vmax = np.abs(Dv.sum(1)).max()
for k in range(K):
    img = Dv[k].sum(1); vk = np.abs(img).max(); ob = OffsetImage(img, cmap='RdBu_r', zoom=1.3 if n < 60 else 0.55, norm=matplotlib.colors.Normalize(-vk, vk)); ab = AnnotationBbox(ob, pos[k], frameon=True, bboxprops=dict(edgecolor=col[k], lw=4), zorder=3); ax.add_artist(ab)
    ax.text(pos[k, 0], pos[k, 1] - 0.16, '%s  %dk' % (lab[k], cnt[k] // 1000), ha='center', va='top', fontsize=8, zorder=4, bbox=dict(facecolor='white', alpha=0.7, lw=0, pad=1))
ax.set_xlim(pos[:, 0].min() - 0.4, pos[:, 0].max() + 0.4); ax.set_ylim(pos[:, 1].min() - 0.4, pos[:, 1].max() + 0.4)
ax.set_title(a.title or 'state graph: nodes = state deviation maps (projection), edges = %d-NN similarity at %g-%g A, frame colour = %s' % (a.knn, a.lo, a.hi, 'dominant conformer (label: conformer, purity)' if a.gt else ('mean sub-class closedness' if a.closedness else 'mean z1')), fontsize=11)
plt.tight_layout(); plt.savefig(a.out, dpi=80); print('figure ->', a.out, '| similarity: median nearest-neighbour %.2f, median overall %.2f' % (np.median(np.sort(S, 1)[:, -1]), np.median(S[np.triu_indices(K, 1)])))
