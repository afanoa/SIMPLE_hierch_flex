#!/usr/bin/env python3
"""Apo/holo x open/closed readout of a klem board on the "ah" toy. Truth from gt_particles.txt (conformer >= HOLO0 = holo)
and closed_ids.txt. Per argmax cell: holo fraction and closed fraction; separability of the two factors: NMI of the
cell labelling with each factor, the share of holo particles in holo-majority cells (and of closed in closed-majority
cells), the 2x2 occupancy of the best cells, and, when the latent carries a layer mass (column 7), the AUC of the
layer mass for holo vs apo and for closed vs open. usage: score_ah.py klem_latent.txt gt_particles.txt closed_ids.txt [--holo0 101] [--tag T]"""
import sys, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('latent'); ap.add_argument('gt'); ap.add_argument('closed'); ap.add_argument('--holo0', type=int, default=101); ap.add_argument('--tag', default=''); a = ap.parse_args()
gt = {}
for l in open(a.gt):
    if l.strip() and not l.startswith('#'): w = l.split(); gt[int(w[0])] = int(w[1])
cl = set(int(v) for v in open(a.closed).read().strip().split(','))
L = np.loadtxt(a.latent, comments='#'); pind = L[:, 0].astype(int); am = L[:, 3].astype(int) - 1; gm = L[:, 5]; lay = L[:, 6] if L.shape[1] > 6 else None
conf = np.array([gt.get(p, 0) for p in pind]); on = (gm > 0.5) & (conf > 0); holo = conf >= a.holo0; closed = np.array([c in cl for c in conf])
def nmi(x, y):
    xs, xi = np.unique(x, return_inverse=True); ys, yi = np.unique(y, return_inverse=True); J = np.zeros((len(xs), len(ys))); np.add.at(J, (xi, yi), 1); J /= J.sum()
    px, py = J.sum(1), J.sum(0); nz = J > 0; I = np.sum(J[nz] * np.log(J[nz] / (px[:, None] * py[None, :])[nz])); hx = -np.sum(px[px > 0] * np.log(px[px > 0])); hy = -np.sum(py[py > 0] * np.log(py[py > 0]))
    return I / max(np.sqrt(hx * hy), 1e-12)
def auc(score, pos):
    r = np.argsort(np.argsort(score)) + 1.; npos = pos.sum(); nneg = len(pos) - npos
    return float((r[pos].sum() - npos * (npos + 1) / 2) / max(npos * nneg, 1)) if npos and nneg else float('nan')
cells = np.unique(am[on]); T = a.tag + ': ' if a.tag else ''
print('%s%d on-grid conformer particles, %d used cells; holo share %.2f, closed share %.2f' % (T, on.sum(), len(cells), holo[on].mean(), closed[on].mean()))
for name, fac in (('holo', holo), ('closed', closed)):
    fr = np.array([fac[on & (am == c)].mean() for c in cells]); nn = np.array([(on & (am == c)).sum() for c in cells]); maj = fr > 0.5
    cover = fac[on & np.isin(am, cells[maj])].sum() / max(fac[on].sum(), 1); pure = fac[on & np.isin(am, cells[maj])].mean() if maj.any() else float('nan')
    print('  %-6s: cell fractions min %.2f median %.2f max %.2f | %d majority cells hold %.2f of the %s particles at purity %.2f | NMI(cell, %s) %.3f%s' % (name, fr.min(), np.median(fr), fr.max(), maj.sum(), cover, name, pure, name, nmi(am[on], fac[on]), '' if lay is None else ' | layer-mass AUC %.3f' % auc(lay[on], fac[on])))
q = {(h, c): (on & (holo == h) & (closed == c)).sum() for h in (0, 1) for c in (0, 1)}
best = {}
for (h, c), nq in q.items():
    m = on & (holo == h) & (closed == c)
    if m.sum() == 0: continue
    fr = np.array([(m & (am == k)).sum() / max((on & (am == k)).sum(), 1) for k in cells]); k = int(np.argmax(fr)); best[(h, c)] = (cells[k], fr[k], (m & (am == cells[k])).sum() / nq)
print('  2x2 (apo/holo x open/closed): best cell purity and recall per quadrant: ' + ' | '.join('%s-%s: cell %d purity %.2f recall %.2f (n %d)' % ('holo' if h else 'apo', 'closed' if c else 'open', v[0], v[1], v[2], q[(h, c)]) for (h, c), v in sorted(best.items())))
print('  NMI(cell, conformer) %.3f; NMI(cell, holo x closed) %.3f' % (nmi(am[on], conf[on]), nmi(am[on], holo[on] * 2 + closed[on])))
