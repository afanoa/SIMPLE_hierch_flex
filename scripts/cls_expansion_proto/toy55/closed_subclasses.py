#!/usr/bin/env python3
"""Where do the CLOSED sub-classes sit in the klem latent? Score every cls_expansion sub-class average for closedness
against the deposited closed (18) vs open (28) conformers rendered as pose-mean projections at the members' refined
poses (klem's calibrated convention: conv 0, shift sign -1, scale 1, |CTF| stacks -> CTF-free averages), then compare
the score with the sub-class mean position along the main latent axis (subclass_axis.txt from the full run).
usage: closed_subclasses.py META CAVGS DEPOSITED_DIR CONSENSUS AXIS_TXT OUT_PREFIX [--lo 40 --hi 14 --nmem 300]"""
import sys, os, glob, argparse, numpy as np
from scipy.ndimage import map_coordinates
from scipy.stats import spearmanr
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'subclass_pool')); sys.path.insert(0, '/mnt/beegfs/elmlund/afan/routes/pvbag/subclass_pool'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from view_tangent_toy import euler2m, Band, fft2c, sample_plane
ap = argparse.ArgumentParser(); ap.add_argument('meta'); ap.add_argument('cavgs'); ap.add_argument('dep'); ap.add_argument('cons'); ap.add_argument('axis'); ap.add_argument('out')
ap.add_argument('--lo', type=float, default=40.); ap.add_argument('--hi', type=float, default=14.); ap.add_argument('--nmem', type=int, default=300); ap.add_argument('--closed', type=int, default=18); ap.add_argument('--open', type=int, default=28)
a = ap.parse_args(); SMPD = 1.464; rng = np.random.default_rng(0)
def read_mrc(p):
    h = np.fromfile(p, dtype=np.int32, count=256); nx, ny, nz = h[0], h[1], h[2]; return np.fromfile(p, dtype=np.float32, offset=1024 + h[23]).reshape(nz, ny, nx).astype(float)
def fftc(v): return np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v)))
def slice_batch(VF, Rs, kx, ky):   # conv 0
    c = VF.shape[0] // 2; hk = np.stack([kx, ky, np.zeros_like(kx)], 1); loc = np.einsum('kj,mji->mki', hk, Rs); coords = [loc[..., 2].ravel() + c, loc[..., 1].ravel() + c, loc[..., 0].ravel() + c]
    return (map_coordinates(VF.real, coords, order=1, mode='nearest') + 1j * map_coordinates(VF.imag, coords, order=1, mode='nearest')).reshape(len(Rs), -1)
M = np.load(a.meta); flags = M['flags'].astype(int); ref = np.where(flags == 1)[0]; lab = M['lab692']; par = M['parent692']
e1, e2, e3, x, y = (M[k][ref] for k in ('e1', 'e2', 'e3', 'x', 'y'))
box = 160; kmin = int(np.floor(box * SMPD / a.lo)); kmax = int(np.floor(box * SMPD / a.hi)); Bd = Band(box, kmax); inb = Bd.r >= kmin; kx, ky, wq = Bd.kx[inb], Bd.ky[inb], Bd.wq[inb]
zz, yy, xx = np.indices((box, box, box)); r = np.hypot(np.hypot(xx - box / 2, yy - box / 2), zz - box / 2); msk = 1. / (1. + np.exp((r - 100 / SMPD) / 2.))
Vc = fftc(read_mrc(a.cons) * msk); V18 = fftc(read_mrc(glob.glob('%s/conf%02d_*.mrc' % (a.dep, a.closed))[0]) * msk); V28 = fftc(read_mrc(glob.glob('%s/conf%02d_*.mrc' % (a.dep, a.open))[0]) * msk)
h = np.fromfile(a.cavgs, np.int32, 256); nimg = int(h[2]); stk = np.memmap(a.cavgs, np.float32, 'r', offset=1024 + int(h[23]), shape=(nimg, box, box))
if a.axis != 'none':
    ax_ = np.loadtxt(a.axis); t_of = {int(s): (m, se, n) for s, p, n, m, se in ax_}
else: t_of = {int(s): (0., 0., 0) for s in np.unique(lab) if s > 0}   # no latent axis: closedness only
subs = [s for s in np.unique(lab) if s > 0 and s <= nimg and int(s) in t_of]; print('%d sub-classes with a latent mean; band %g-%g A (%d coefficients); %d members rendered per sub-class' % (len(subs), a.lo, a.hi, len(kx), a.nmem), flush=True)
res = []
for j, s in enumerate(subs):
    mem = np.where(lab == s)[0]
    if len(mem) > a.nmem: mem = rng.choice(mem, a.nmem, replace=False)
    Rs = np.array([euler2m([-e1[i], e2[i], -e3[i]]) for i in mem]); ph = np.exp(-2j * np.pi * (kx[None, :] * (-x[mem])[:, None] + ky[None, :] * (-y[mem])[:, None]) / box)   # shift sign -1, scale 1
    Pc = (slice_batch(Vc, Rs, kx, ky) * ph).mean(0); P18 = (slice_batch(V18, Rs, kx, ky) * ph).mean(0); P28 = (slice_batch(V28, Rs, kx, ky) * ph).mean(0)
    A = sample_plane(fft2c(np.array(stk[int(s) - 1], float)), kx, ky)
    fit = np.sum(wq * np.real(np.conj(Pc) * A)) / max(np.sum(wq * np.abs(Pc) ** 2), 1e-30); resid = A - fit * Pc; d = P18 - P28
    cc_avg = np.sum(wq * np.real(np.conj(Pc) * A)) / np.sqrt(np.sum(wq * np.abs(Pc) ** 2) * np.sum(wq * np.abs(A) ** 2))
    score = np.sum(wq * np.real(np.conj(d) * resid)) / np.sqrt(np.sum(wq * np.abs(d) ** 2) * np.sum(wq * np.abs(resid) ** 2))
    # sign-free control: the same with the residual against a random other direction is not available; keep the raw agreement of the average with the consensus render as a quality check
    res.append((int(s), int(par[np.where(lab == s)[0][0]]), len(np.where(lab == s)[0]), cc_avg, score, *t_of[int(s)]))
    if j % 100 == 0: print('  %d / %d' % (j, len(subs)), flush=True)
res = np.array(res); s_, p_, n_, ccav, score, mt, se, nlat = res.T
print('average vs pose-mean consensus render: ncc median %.2f (band %g-%g A): the frame/convention check (should be well above 0)' % (np.median(ccav), a.lo, a.hi))
good = ccav > np.percentile(ccav, 10)
print('closedness score (corr of the average residual with the closed-minus-open pose-mean difference): median %.3f, sd %.3f, 5th/95th pct %.2f / %.2f' % (np.median(score), score.std(), np.percentile(score, 5), np.percentile(score, 95)))
print('Spearman(closedness, latent mean t): %+.3f over %d sub-classes (negative t = the bin-1 end)' % (spearmanr(score[good], mt[good])[0], good.sum()))
# in bins of the latent mean: mean closedness
edges = np.percentile(mt[good], np.linspace(0, 100, 6)); b = np.clip(np.searchsorted(edges, mt[good], side='right') - 1, 0, 4)
print('mean closedness by quintile of the latent mean (bin-1 end -> bin-7 end): ' + ' '.join('%.3f' % score[good][b == k].mean() for k in range(5)))
o = np.argsort(-score); print('10 most CLOSED-scoring sub-classes (sub, parent, n, score, latent t):'); [print('  sub %5d parent %3d n %4d score %+.3f t %+.2f' % (s_[k], p_[k], n_[k], score[k], mt[k])) for k in o[:10]]
print('10 most OPEN-scoring:'); [print('  sub %5d parent %3d n %4d score %+.3f t %+.2f' % (s_[k], p_[k], n_[k], score[k], mt[k])) for k in o[-10:][::-1]]
np.savetxt(a.out + '_closedness.txt', res, fmt='%d %d %d %.4f %.4f %.4f %.4f %d', header='subclass parent n ncc_avg_vs_consensus closedness latent_t sem n_latent')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
fig = plt.figure(figsize=(18, 9)); gs = fig.add_gridspec(4, 9)
ax = fig.add_subplot(gs[:, :3]); ax.scatter(mt[good], score[good], s=6, alpha=0.5); ax.set_xlabel('sub-class mean position on the main latent axis (bin-1 end negative)'); ax.set_ylabel('closedness of the sub-class average'); ax.set_title('Spearman %+.2f' % spearmanr(score[good], mt[good])[0]); ax.axhline(0, color='k', lw=0.5)
for row, (title, idx) in enumerate((('most closed-scoring averages', o[:6]), ('most open-scoring averages', o[-6:][::-1]), ('latent bin-1 end (most negative t)', np.argsort(mt)[:6]), ('latent bin-7 end (most positive t)', np.argsort(-mt)[:6]))):
    for col, k in enumerate(idx):
        axx = fig.add_subplot(gs[row, 3 + col]); im = gaussian_filter(np.array(stk[int(s_[k]) - 1], float), 1.5); axx.imshow(im, cmap='gray'); axx.set_xticks([]); axx.set_yticks([])
        axx.set_title('sub %d  c %+.2f  t %+.1f' % (s_[k], score[k], mt[k]), fontsize=8)
        if col == 0: axx.set_ylabel(title, fontsize=8)
plt.tight_layout(); plt.savefig(a.out + '_closedness.png', dpi=100); print('figure ->', a.out + '_closedness.png')
