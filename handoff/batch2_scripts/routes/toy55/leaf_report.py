#!/usr/bin/env python3
"""Per-leaf truth report of an rflex run (synthetic): for every leaf of leaves.txt (name map even odd n labels_file state [junk]) its
size, junk flag, class shares (apo / holo / swing / junk), top conformers and the swing mean +- sd; then NMI(conformer; leaf) and
purity with and without the flagged junk leaves.
usage: leaf_report.py RUNDIR GT"""
import sys, numpy as np
R, gt = sys.argv[1], sys.argv[2]
L = [l.split() for l in open(R + '/leaves.txt') if l.strip()]; G = np.loadtxt(gt, comments='#'); n = int(G[:, 0].max())
conf = np.full(n, -1); conf[G[:, 0].astype(int) - 1] = G[:, 1].astype(int); lab = np.zeros(n, int); jk = []
def cls(c): return np.where(c == 0, 'junk', np.where(c >= 200, 'swing', np.where(c >= 101, 'holo', 'apo')))
for i, l in enumerate(L, 1):
    s = np.loadtxt(l[5], dtype=int); lab[s == int(l[6])] = i; jk.append(len(l) > 7 and l[7] == 'junk')
print('%-12s %7s %4s | %5s %5s %5s %5s | %-30s | swing' % ('leaf', 'n', 'junk', 'apo', 'holo', 'swing', 'junk', 'top conformers'))
for i, l in enumerate(L, 1):
    c = conf[lab == i]; k = cls(c); cc = np.bincount(c[c > 0]) if (c > 0).any() else np.zeros(1, int); top = np.argsort(-cc)[:3]
    sw = c[(c >= 200) & (c <= 230)] - 200
    print('%-12s %7d %4s | %5.2f %5.2f %5.2f %5.2f | %-30s | %s' % (l[0], len(c), 'J' if jk[i - 1] else '', (k == 'apo').mean(), (k == 'holo').mean(), (k == 'swing').mean(), (k == 'junk').mean(),
          ' '.join('%d:%.2f' % (t, cc[t] / len(c)) for t in top if cc[t] > 0), '%.1f+-%.1f' % (sw.mean(), sw.std()) if len(sw) > 50 else ''))
def nmi(x, y):
    _, xi = np.unique(x, return_inverse=True); _, yi = np.unique(y, return_inverse=True); J = np.zeros((xi.max() + 1, yi.max() + 1)); np.add.at(J, (xi, yi), 1); J /= J.sum()
    px, py = J.sum(1), J.sum(0); nz = J > 0; I = (J[nz] * np.log(J[nz] / np.outer(px, py)[nz])).sum(); return I / max(np.sqrt((px * np.log(px)).sum() * (py * np.log(py)).sum()), 1e-12)
on = lab > 0; real = np.isin(lab, [i for i in range(1, len(L) + 1) if not jk[i - 1]])
jl = lab.copy(); jl[np.isin(lab, [i for i in range(1, len(L) + 1) if jk[i - 1]])] = -1   # all flagged junk leaves as one label
print('leaves %d (%d junk): NMI(conformer; leaf) %.3f; flagged junk leaves as one label %.3f; non-junk leaves only %.3f' % (len(L), sum(jk), nmi(conf[on], lab[on]), nmi(conf[on], jl[on]), nmi(conf[real], lab[real])))
tj = conf == 0; fj = on & ~real
print('junk flag: recall %.3f (%d of %d true junk), precision %.3f (%d flagged)' % ((tj & fj).sum() / max(tj.sum(), 1), (tj & fj).sum(), tj.sum(), (tj & fj).sum() / max(fj.sum(), 1), fj.sum()))
