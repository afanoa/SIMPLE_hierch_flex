#!/usr/bin/env python3
"""land_runcmp.py: do the per-run pole-to-pole paths (land_runpath.sbatch) give the same occupancy landscape? Each run's path is oriented
so that s = 0 is its closed-leaning pole (arm.txt: higher closed-site occupancy). Reports, between runs and against each run's own
halves: the Bhattacharyya coefficient and L1 distance of the 31-bin occupancies p(s), the mass on the closed side / middle / open end,
the Spearman correlation of the particles' posterior mean positions, and the agreement of their sides of the sparsest point.
usage: land_runcmp.py DIR1,DIR2,... name1,name2,... [--fig PNG]   (DIR = rnasep_klem/landscape/runpath_<run>)"""
import sys, json, numpy as np
from scipy.stats import spearmanr
dirs = sys.argv[1].split(','); names = sys.argv[2].split(','); fig = sys.argv[sys.argv.index('--fig') + 1] if '--fig' in sys.argv else ''; XL = sys.argv[sys.argv.index('--xlabel') + 1] if '--xlabel' in sys.argv else 's (closed pole 0 -> open pole 1)'; TT = sys.argv[sys.argv.index('--title') + 1] if '--title' in sys.argv else "r1 along each run's own pole-to-pole path (lam 0.001; halves dotted/dashed)"
R = {}
for d, nm in zip(dirs, names):
    arm = [l.split() for l in open(d + '/arm.txt') if not l.startswith('#')]; cl = [float(x[1]) for x in arm]; flip = (cl[0] < cl[1]) and '--noflip' not in sys.argv
    J = json.load(open(d + '/rw_lmarg_split_0.001/reweight_real.json')); Z = np.load(d + '/rw_lmarg_split_0.001/posterior.npz')
    p = np.array(J['p']); h = [np.array(x) for x in J['halves']]; sm = Z['smean'].astype(float)
    if flip: p, h, sm = p[::-1], [x[::-1] for x in h], 1. - sm
    R[nm] = dict(p=p / p.sum(), h=[x / x.sum() for x in h], s=dict(zip(Z['rows'].tolist(), sm.tolist())), flip=flip)
    import os
    if not os.path.isfile(d + '/pc1.txt'): pc = None
    else: pc = np.loadtxt(d + '/pc1.txt')
    rr = Z['rows'].tolist()
    if pc is not None:
      pcs = dict(zip((pc[:, 0].astype(int) - 1).tolist(), pc[:, 1].tolist()))
      a_ = np.array([R[nm]['s'][r] for r in rr if r in pcs]); b_ = np.array([pcs[r] for r in rr if r in pcs]); rho = spearmanr(a_, b_)[0]
      print('%-6s: Spearman(posterior mean s, own latent PC1) %+.3f; posterior mean s quantiles 10/25/50/75/90 %% %s' % (nm, rho, ' '.join('%.2f' % x for x in np.quantile(a_, [0.1, 0.25, 0.5, 0.75, 0.9]))))
    K = len(p); sg = np.arange(K) / (K - 1.); G = -np.log(np.maximum(p, 1e-12)); inner = slice(3, K - 3); sp = 3 + int(np.argmax(G[inner]))
    R[nm]['sp'] = sg[sp]
    print('%-6s (path %s): closed side s<=0.3 %.3f | middle %.3f | open end s>=0.83 %.3f | sparsest interior point s %.2f (G %.2f above the deepest) | end bins %.3f / %.3f'
          % (nm, 'flipped' if flip else 'as is', p[sg <= 0.3].sum(), p[(sg > 0.3) & (sg < 0.83)].sum(), p[sg >= 0.83].sum(), sg[sp], G[sp] - G.min(), p[0], p[-1]))
bc = lambda a, b: float(np.sqrt(a * b).sum())
for nm in names: print('%-6s halves: BC %.4f, L1 %.3f' % (nm, bc(*R[nm]['h']), np.abs(R[nm]['h'][0] - R[nm]['h'][1]).sum()))
for i, x in enumerate(names):
    for y in names[i + 1:]:
        com = sorted(set(R[x]['s']) & set(R[y]['s'])); a = np.array([R[x]['s'][c] for c in com]); b = np.array([R[y]['s'][c] for c in com])
        cut = 0.5 * (R[x]['sp'] + R[y]['sp']); agree = ((a > cut) == (b > cut)).mean()
        print('%s vs %s: p(s) BC %.4f, L1 %.3f | particles %d: Spearman of posterior mean s %.3f, same side of s %.2f %.3f' % (x, y, bc(R[x]['p'], R[y]['p']),
              np.abs(R[x]['p'] - R[y]['p']).sum(), len(com), spearmanr(a, b)[0], cut, agree))
if fig:
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    f, ax = plt.subplots(1, 2, figsize=(10, 3.4)); col = dict(zip(names, ['#2a6fdb', '#d9822b', '#3a9d5d', '#8e5bd0']))
    for nm in names:
        K = len(R[nm]['p']); sg = np.arange(K) / (K - 1.); ax[0].plot(sg, R[nm]['p'], color=col[nm], lw=1.8, label=nm)
        for hh, ls in zip(R[nm]['h'], (':', '--')): ax[0].plot(sg, hh, color=col[nm], lw=0.8, ls=ls)
        G = -np.log(np.maximum(R[nm]['p'], 1e-12)); ax[1].plot(sg, G - G.min(), color=col[nm], lw=1.8)
    ax[0].set_xlabel(XL); ax[0].set_ylabel('occupancy p(s)'); ax[0].legend(frameon=False, fontsize=8)
    ax[1].set_xlabel('s'); ax[1].set_ylabel('G = -ln p (kT, apparent)'); f.suptitle(TT, fontsize=9)
    f.tight_layout(); f.savefig(fig, dpi=150); print('figure', fig)
