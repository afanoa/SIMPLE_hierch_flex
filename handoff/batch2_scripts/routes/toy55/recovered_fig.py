#!/usr/bin/env python3
"""Figure: the ten deposited RNase P conformers recovered from the real data. One column per conformer: our full-box
map (projection, 8 A), our map's deviation from the set mean and the deposited conformer's deviation from the
consensus (signed projections, 15 A, same colour scale), along z and along y. Column titles carry the own-conformer
correlation at 60-20 A, the FSC resolution and the closed-arm occupancy.
usage: recovered_fig.py MAPS_DIR OUT.png"""
import sys, glob, numpy as np
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/routes/toy55'); from model_class import read_mrc, Shells, fftc
from scipy.ndimage import gaussian_filter
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
maps, out = sys.argv[1], sys.argv[2]; W = '/mnt/beegfs/elmlund/afan/routes/toy55'; SMPD = 1.464; ids = [1, 10, 18, 27, 28, 37, 46, 48, 49, 51]
S = Shells(160, SMPD); msk = S.soft_mask(100 / SMPD); sh = S.sh; n = 160; zz, yy, xx = np.indices((n, n, n)); rr = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2)
def lp(v, res): return np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v * msk))) * (rr <= n * SMPD / res)))))
cons = read_mrc(W + '/maps160/toy_consensus.mrc'); C = {i: read_mrc(glob.glob(W + '/maps160/docked/conf%02d_*.mrc' % i)[0]) for i in ids}
V = [read_mrc('%s/recvol_state%02d.mrc' % (maps, k + 1)) for k in range(10)]; Vm = np.mean(V, 0)
def bc(x, y, lo, hi):
    m = S.band(lo, hi); return float(np.real(np.sum(np.conj(x[m]) * y[m])) / max(np.sqrt(np.sum(np.abs(x[m]) ** 2) * np.sum(np.abs(y[m]) ** 2)), 1e-30))
def fsc_res(x, y):
    a_, b_ = fftc(x * msk), fftc(y * msk); num = np.bincount(sh.ravel(), np.real(np.conj(a_) * b_).ravel()); d1 = np.bincount(sh.ravel(), np.abs(a_.ravel()) ** 2); d2 = np.bincount(sh.ravel(), np.abs(b_.ravel()) ** 2)
    f = num / np.maximum(np.sqrt(d1 * d2), 1e-30); s = 1
    while s < 80 and f[s] > 0.143: s += 1
    return n * SMPD / s
g = lambda v: gaussian_filter(v, 2.5); L18, Lc = g(C[18]), g(cons); dens18 = L18 > L18.mean() + 2.0 * L18.std(); diff = L18 - Lc
arm = (diff > diff.mean() + 3 * diff.std()) & dens18; core = dens18 & (np.abs(diff) < 0.5 * diff.std()); ref = L18[arm].mean() / L18[core].mean()
T = {i: fftc((C[i] - cons) * msk) for i in ids}
rows = [[], [], [], [], []]; titles = []
for k, i in enumerate(ids):
    own = bc(fftc((V[k] - Vm) * msk), T[i], 60, 20)
    try: res = fsc_res(read_mrc('%s/recvol_state%02d_even.mrc' % (maps, k + 1)), read_mrc('%s/recvol_state%02d_odd.mrc' % (maps, k + 1)))
    except Exception: res = float('nan')
    occ = (g(V[k])[arm].mean() / g(V[k])[core].mean()) / ref; titles.append('conf %d\ncorr %+.2f | %.1f A | arm %.2f' % (i, own, res, occ))
    dv = lp(V[k] - Vm, 15.); dc = lp(C[i] - cons, 15.); rows[0].append(lp(V[k], 8.).sum(0)); rows[1].append(dv.sum(0)); rows[2].append(dc.sum(0)); rows[3].append(dv.sum(1)); rows[4].append(dc.sum(1))
labels = ['our map, proj z (8 A)', 'ours: map - set mean, proj z (15 A)', 'deposited: conformer - consensus, proj z', 'ours: map - set mean, proj y', 'deposited: conformer - consensus, proj y']
fig, axes = plt.subplots(5, 10, figsize=(24, 12.5), squeeze=False)
for r in range(5):
    imgs = rows[r]
    if r == 0:
        for j, img in enumerate(imgs): axes[r, j].imshow(img, cmap='gray', vmin=0, vmax=np.percentile(img, 99.7))
    else:
        vmax = np.percentile(np.abs(np.array(imgs)), 99.5)   # one symmetric scale per row (model-derived conformers and reconstructions differ in amplitude)
        for j, img in enumerate(imgs): axes[r, j].imshow(img, cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    for j in range(10):
        axes[r, j].set_xticks([]); axes[r, j].set_yticks([])
        if r == 0: axes[r, j].set_title(titles[j], fontsize=9)
    axes[r, 0].set_ylabel(labels[r], fontsize=8)
plt.suptitle('RNase P: the ten deposited conformers recovered from the real data (family 1 of the klem tree, refined poses; corr = own-conformer deviation correlation at 60-20 A, set mean subtracted)', fontsize=11)
plt.tight_layout(); plt.savefig(out, dpi=100); print('figure ->', out)
