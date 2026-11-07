#!/usr/bin/env python3
"""Do a node's states differ in their view distributions (problem 3: a cause of molecule-wide reproducible differences that a
per-shell regression cannot remove: each state's map carries its own anisotropic resolution)? Per state of a labelling and a
print_project_field dump of the poses the maps were reconstructed at: view histogram over NDIR Fibonacci directions (projection
direction from e1, e2; antipodes not folded), Bhattacharyya overlap with the node's histogram (1 = same views) and the KL divergence,
defocus mean; the same for a random split of the same sizes (the sampling floor).
usage: p3_views.py DUMP LABELS [--name N --ndir 100 --seed 1]"""
import re, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('dump'); ap.add_argument('labels'); ap.add_argument('--name', default='set'); ap.add_argument('--ndir', type=int, default=100); ap.add_argument('--seed', type=int, default=1)
a = ap.parse_args()
pat = {k: re.compile(r'(?:^| )%s=([-0-9.eE+]+)' % k) for k in ('e1', 'e2', 'dfx', 'dfy')}
X = np.array([[float(pat[k].search(l).group(1)) for k in ('e1', 'e2', 'dfx', 'dfy')] for l in open(a.dump) if l.strip()])
lab = np.loadtxt(a.labels, dtype=int); assert len(lab) == len(X), (len(lab), len(X)); on = lab > 0; lab, X = lab[on], X[on]
e1, e2 = np.deg2rad(X[:, 0]), np.deg2rad(X[:, 1]); d = np.c_[np.sin(e2) * np.cos(e1), np.sin(e2) * np.sin(e1), np.cos(e2)]
i = np.arange(a.ndir) + 0.5; z = 1 - 2 * i / a.ndir; ph = np.pi * (1 + 5 ** 0.5) * i; F = np.c_[np.sqrt(1 - z ** 2) * np.cos(ph), np.sqrt(1 - z ** 2) * np.sin(ph), z]
b = (d @ F.T).argmax(1); hn = np.bincount(b, minlength=a.ndir) / len(b); df = (X[:, 2] + X[:, 3]) / 2
def stats(l):
    out = []
    for s in np.unique(l):
        m = l == s; h = np.bincount(b[m], minlength=a.ndir) / m.sum(); bc = np.sqrt(h * hn).sum(); kl = (h[h > 0] * np.log(h[h > 0] / hn[h > 0])).sum()
        out.append((s, m.sum(), bc, kl, df[m].mean()))
    return out
rng = np.random.default_rng(a.seed); st = stats(lab); sr = stats(rng.permutation(lab))
print('%s: %d particles, %d states; view overlap with the node (Bhattacharyya) min / median %.3f / %.3f (random split %.3f / %.3f), KL max %.3f (random %.3f); defocus means %.2f .. %.2f um (random %.2f .. %.2f)'
      % (a.name, len(lab), len(st), min(r[2] for r in st), np.median([r[2] for r in st]), min(r[2] for r in sr), np.median([r[2] for r in sr]), max(r[3] for r in st), max(r[3] for r in sr),
         min(r[4] for r in st), max(r[4] for r in st), min(r[4] for r in sr), max(r[4] for r in sr)))
print('   per state (n, overlap, KL, defocus): ' + ' '.join('%d:%d/%.3f/%.3f/%.2f' % r for r in st))
