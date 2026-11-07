#!/usr/bin/env python3
"""Per-particle readouts of a flex_pca latent against the toy truth (10 conformers + compositional variants + junk).
usage: ptcl_latent_eval.py FLEXDIR GT --states 1,10,.. --axis a1,a2,.. [--sublab ptcl2D.txt] [--subset states.txt]
                           [--out coord.txt] [--tag name] [--nmax 0]
GT = gt_particles.txt (particle conformer closed e1 e2; conformer 0 = junk). Readouts:
  axis      : Spearman of every PC and of a cross-fitted ridge combination with the truth conformer axis (ordering)
  knn       : cross-validated 15-NN accuracy over the conformers (chance = largest conformer share), balanced accuracy
  fisher    : pairwise d' between conformers in the latent (median / min), closed-vs-rest d' and AUC
  junk      : AUC of the latent norm for junk detection
  pooled    : the same at the 2D sub-class level (mean latent per sub-class, --sublab = ptcl2D table with class=)
  states    : delivered flex_pca states (flex_pca_state_weights.txt argmax) vs conformers (Hungarian), closed recall
OUT: 'pind coord_pc1 coord_sup' = the unsupervised PC1 (window seed) and the cross-fitted supervised axis (oracle seed)."""
import sys, os, argparse
import numpy as np
from scipy.stats import spearmanr
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/routes/pvbag/subclass_pool'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pool_place import read_embedding


def read_gt(fn):
    gt = {}
    for l in open(fn):
        if l.strip() and not l.startswith('#'): w = l.split(); gt[int(w[0])] = (int(w[1]), int(w[2]))
    return gt


def hungarian_acc(pred, truth):
    pu, pi = np.unique(pred, return_inverse=True); tu, ti = np.unique(truth, return_inverse=True)
    C = np.zeros((len(pu), len(tu))); np.add.at(C, (pi, ti), 1); r, c = linear_sum_assignment(-C)
    return C[r, c].sum() / len(pred), {int(pu[i]): int(tu[j]) for i, j in zip(r, c)}, C


def knn_cv(z, y, k=15, nmax=60000, rng=None):
    idx = np.arange(len(y)) if len(y) <= nmax else rng.choice(len(y), nmax, replace=False)
    z, y = z[idx], y[idx]; h = np.arange(len(y)) % 2; acc = []; bal = []
    for a in (0, 1):
        tr, te = h == a, h != a; kk = int(min(k, tr.sum())); tree = cKDTree(z[tr]); _, nn = tree.query(z[te], k=kk); lab = y[tr][np.atleast_2d(nn.T).T if kk > 1 else nn[:, None]]
        pred = np.array([np.bincount(r).argmax() for r in lab]); acc.append(np.mean(pred == y[te]))
        bal.append(np.mean([np.mean(pred[y[te] == c] == c) for c in np.unique(y[te])]))
    return float(np.mean(acc)), float(np.mean(bal))


def fisher(z, y):
    cls = np.unique(y); mu = np.array([z[y == c].mean(0) for c in cls]); W = sum(np.cov((z[y == c] - mu[i]).T, bias=True) * np.sum(y == c) for i, c in enumerate(cls)) / len(y)
    W = np.atleast_2d(W) + 1e-6 * np.trace(np.atleast_2d(W)) / z.shape[1] * np.eye(z.shape[1]); Wi = np.linalg.inv(W)
    d = np.array([[np.sqrt(max((mu[i] - mu[j]) @ Wi @ (mu[i] - mu[j]), 0)) for j in range(len(cls))] for i in range(len(cls))])
    iu = np.triu_indices(len(cls), 1); return d[iu], cls, d


def auc(score, pos):
    r = np.argsort(np.argsort(score)) + 1; npos = pos.sum(); nneg = len(pos) - npos
    return float((r[pos].sum() - npos * (npos + 1) / 2) / max(npos * nneg, 1))


def lda_dir(z, pos):
    mu1, mu0 = z[pos].mean(0), z[~pos].mean(0); W = np.cov(z[pos].T, bias=True) * pos.sum() + np.cov(z[~pos].T, bias=True) * (~pos).sum(); W = np.atleast_2d(W) / len(z)
    W += 1e-6 * np.trace(W) / z.shape[1] * np.eye(z.shape[1]); w = np.linalg.solve(W, mu1 - mu0); return w / np.linalg.norm(w)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('flexdir'); ap.add_argument('gt'); ap.add_argument('--states', required=True); ap.add_argument('--axis', required=True)
    ap.add_argument('--sublab', default=''); ap.add_argument('--subset', default=''); ap.add_argument('--out', default=''); ap.add_argument('--tag', default='latent'); ap.add_argument('--nmax', type=int, default=60000)
    a = ap.parse_args(); rng = np.random.default_rng(0); T = a.tag
    states = [int(s) for s in a.states.split(',')]; axis = dict(zip(states, [float(v) for v in a.axis.split(',')]))
    _, pinds, z, _ = read_embedding(a.flexdir + '/flex_pca_embedding.bin'); gt = read_gt(a.gt)
    if a.subset:
        sub = np.loadtxt(a.subset).astype(int); keep = np.array([sub[p - 1] > 0 for p in pinds]); pinds, z = pinds[keep], z[keep]
    conf = np.array([gt.get(int(p), (0, 0))[0] for p in pinds]); closed = np.array([gt.get(int(p), (0, 0))[1] for p in pinds]) == 1
    good = conf > 0; zs = (z - z.mean(0)) / np.maximum(z.std(0), 1e-12); nc = z.shape[1]
    print('%s: %d particles embedded (%d conformer, %d junk), %d latent dims; eigen-spread of PCs %s' % (T, len(pinds), good.sum(), (~good).sum(), nc, np.array2string(z.std(0)[:6] / z.std(0)[0], precision=2)))
    # ---- junk detection by latent norm
    if (~good).sum() > 0: print('%s junk: AUC of |z| for junk vs conformer particles %.3f' % (T, auc(np.linalg.norm(zs, axis=1), ~good)))
    zg, cg, clg = zs[good], conf[good], closed[good]; ax = np.array([axis[c] for c in cg])
    # ---- ordering: Spearman per PC, cross-fitted ridge combination
    rho = [spearmanr(zg[:, q], ax)[0] for q in range(min(nc, 8))]
    h = np.arange(len(cg)) % 2; sup = np.zeros(len(cg))
    for s in (0, 1):
        tr, te = h == s, h != s; X = np.c_[zg[tr], np.ones(tr.sum())]; w = np.linalg.solve(X.T @ X + 1e-3 * np.eye(X.shape[1]), X.T @ ax[tr]); sup[te] = np.c_[zg[te], np.ones(te.sum())] @ w
    rs = spearmanr(sup, ax)[0]
    print('%s axis: Spearman per PC %s | cross-fitted combination %.3f' % (T, ' '.join('%+.2f' % r for r in rho), rs))
    ladder = sorted(((zg[cg == c, 0].mean(), c) for c in np.unique(cg)))
    print('%s ladder (mean PC1 per conformer, ascending): %s; truth axis in that order: %s' % (T, ' '.join('%d' % c for _, c in ladder), ' '.join('%+.2f' % axis[c] for _, c in ladder)))
    # ---- separability
    acc, bal = knn_cv(zg, cg, nmax=a.nmax, rng=rng); share = np.bincount(cg).max() / len(cg)
    dpair, cls, dmat = fisher(zg, cg)
    print("%s knn: 15-NN accuracy %.3f (largest conformer share %.3f), balanced %.3f | pairwise d' median %.2f min %.2f max %.2f" % (T, acc, share, bal, np.median(dpair), dpair.min(), dpair.max()))
    if clg.sum() > 0:
        w = lda_dir(zg, clg); sc = zg @ w; print("%s closed (%d = %.3f): LDA d' %.2f, AUC %.3f; mean PC1 closed vs open %+.2f vs %+.2f" % (T, clg.sum(), clg.mean(), abs(sc[clg].mean() - sc[~clg].mean()) / np.sqrt(0.5 * (sc[clg].var() + sc[~clg].var())), auc(sc, clg), zg[clg, 0].mean(), zg[~clg, 0].mean()))
    for c in np.unique(cg):
        i = list(cls).index(c); print("%s   conformer %2d (n %6d, axis %+.2f): mean PC1 %+.2f, nearest other conformer d' %.2f" % (T, c, (cg == c).sum(), axis[c], zg[cg == c, 0].mean(), np.min(np.delete(dmat[i], i))))
    # ---- delivered states
    sw = a.flexdir + '/flex_pca_state_weights.txt'
    if os.path.exists(sw):
        rows = np.loadtxt(sw, comments='#'); rows = np.atleast_2d(rows); lab = dict(zip(rows[:, 0].astype(int), rows[:, 1:].argmax(1) + 1))
        pl = np.array([lab.get(int(p), 0) for p in pinds]); m = good & (pl > 0)
        if m.sum() > 0:
            hacc, mp, C = hungarian_acc(pl[m], conf[m]); sizes = np.bincount(pl[m]); cs = closed[m]
            best = max(((np.sum((pl[m] == k) & cs), k) for k in np.unique(pl[m])))
            print('%s states: %d delivered; accuracy vs conformers (Hungarian) %.3f; sizes %s; closed recall %.2f at purity %.2f (state %d)' %
                  (T, len(np.unique(pl[m])), hacc, ' '.join('%d:%d' % (k, sizes[k]) for k in np.unique(pl[m])), best[0] / max(cs.sum(), 1), best[0] / max(np.sum(pl[m] == best[1]), 1), best[1]))
    # ---- pooled by 2D sub-class
    if a.sublab:
        sc = {}
        for i, l in enumerate(open(a.sublab)):
            kv = dict(t.split('=', 1) for t in l.split() if '=' in t); c = int(float(kv.get('class', 0)))
            if c > 0: sc[i + 1] = c
        sl = np.array([sc.get(int(p), 0) for p in pinds]); m = good & (sl > 0); su, si = np.unique(sl[m], return_inverse=True); pop = np.bincount(si)
        ok = pop >= 20; zb = np.array([zs[m][si == j].mean(0) for j in range(len(su))]); axb = np.array([ax[m[good]][si == j].mean() if False else np.mean([axis[c] for c in conf[m][si == j]]) for j in range(len(su))])
        maj = np.array([np.bincount(conf[m][si == j]).argmax() for j in range(len(su))]); pur = np.array([np.mean(conf[m][si == j] == maj[j]) for j in range(len(su))])
        zb, axb, maj, pur, pop = zb[ok], axb[ok], maj[ok], pur[ok], pop[ok]
        rho_b = [spearmanr(zb[:, q], axb)[0] for q in range(min(nc, 4))]
        accb, balb = knn_cv(zb, maj, k=5, nmax=a.nmax, rng=rng); dpb, _, _ = fisher(zb, maj) if len(np.unique(maj)) > 1 else (np.array([0.]), None, None)
        print("%s pooled (%d sub-classes >= 20 members, mean purity %.2f): Spearman per PC %s | 5-NN accuracy vs majority conformer %.3f (balanced %.3f) | pairwise d' median %.2f" %
              (T, len(zb), pur.mean(), ' '.join('%+.2f' % r for r in rho_b), accb, balb, np.median(dpb)))
    # ---- coordinate output: PC1 (unsupervised) and the cross-fitted supervised axis (oracle); junk gets PC1 and 0
    if a.out:
        supall = np.zeros(len(pinds)); supall[good] = sup
        if (~good).any(): supall[~good] = rng.choice(sup, (~good).sum(), replace=True)   # junk: a coordinate drawn from the conformer particles' distribution (no pile-up at 0)
        with open(a.out, 'w') as f:
            f.write('# pind coord_pc1 coord_sup\n')
            for p, z1, s in zip(pinds, zs[:, 0], supall): f.write('%d %.5f %.5f\n' % (p, z1, s))
        print('%s coordinates -> %s' % (T, a.out))


if __name__ == '__main__':
    main()
