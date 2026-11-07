#!/usr/bin/env python3
import argparse, glob, json, os, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--truth', required=True); ap.add_argument('--rows', required=True); ap.add_argument('--json', default='')
for k in ('flex', 'ref', 'em', 'labels'): ap.add_argument('--' + k, nargs='*', default=[])
a = ap.parse_args()
hd = [l.split()[1:] for l in open(a.truth) if l.startswith('#')]; names = hd[0][3:]; kinds = hd[1][1:]; T = np.genfromtxt(a.truth, comments='#'); X = T[:, 3:]; NM = len(T)
use = (np.loadtxt(a.rows) > 0) & (T[:, 2] == 0); grp = [n.split('_')[0] for n in names]; binv = [k == 'd' for k in kinds]
def coords(d):
    f = sorted(glob.glob(d + '/*_flex_pca/flex_pca_coordinates.txt'), key=os.path.getmtime)[-1]; C = np.loadtxt(f, comments='#'); r = C[:, 0].astype(int) - 1
    Z = np.full((NM, C.shape[1] - 5), np.nan); Z[r] = C[:, 5:]; eo = np.full(NM, -1); eo[r] = C[:, 1]; lab = np.zeros(NM, int); lab[r] = C[:, 2]; return Z, eo, lab
def auc(s, y):
    rk = np.argsort(np.argsort(s)) + 1.; n1 = y.sum(); n0 = len(y) - n1; return (rk[y].sum() - n1 * (n1 + 1) / 2) / max(n1 * n0, 1)
def cv(F, v, eo, b):
    p = np.zeros(len(v))
    for f in (0, 1):
        tr, te = eo != f, eo == f; w = np.linalg.lstsq(np.c_[np.ones(tr.sum()), F[tr]], v[tr], rcond=None)[0]; p[te] = np.c_[np.ones(te.sum()), F[te]] @ w
    return auc(p, v > 0.5) if b else 1 - ((v - p) ** 2).sum() / ((v - v.mean()) ** 2).sum()
def eta2(lab, v):
    m = v.mean(); return sum((lab == k).sum() * (v[lab == k].mean() - m) ** 2 for k in np.unique(lab)) / (len(v) * v.var())
quad = lambda Z: np.c_[Z, np.stack([Z[:, i] * Z[:, j] for i in range(Z.shape[1]) for j in range(i, Z.shape[1])], 1)]
one = lambda s, v, b: max(auc(s, v > 0.5), 1 - auc(s, v > 0.5)) if b else np.corrcoef(s, v)[0, 1] ** 2
res = {}; ZS = []; EO = None
def score_flex(nm, d, kind):
    global EO
    Z, eo, lab = coords(d); ok = use & np.isfinite(Z[:, 0]) & (eo >= 0); out = {}
    if kind == 'site': ZS.append(Z); EO = eo
    for j, vn in enumerate(names):
        m = ok & np.isfinite(X[:, j]); v = X[m, j]
        if len(np.unique(v)) < 2: continue
        out[vn] = dict(lin=cv(Z[m], v, eo[m], binv[j]), quad=cv(quad(Z[m]), v, eo[m], binv[j]), single=max(one(Z[m, k], v, binv[j]) for k in range(Z.shape[1])), states=eta2(lab[m], v), n=int(m.sum()))
    res[nm] = dict(kind=kind, vars=out)
for kind, ents in (('site', a.flex), ('ref', a.ref)):
    for e in ents: score_flex(*e.split('=', 1), kind)
if len(ZS) > 1:
    Z = np.concatenate(ZS, 1); ok = use & np.all(np.isfinite(Z), 1) & (EO >= 0); out = {}
    for j, vn in enumerate(names):
        m = ok & np.isfinite(X[:, j]); v = X[m, j]
        if len(np.unique(v)) > 1: out[vn] = dict(lin=cv(Z[m], v, EO[m], binv[j]), n=int(m.sum()))
    res['union_of_sites'] = dict(kind='union', vars=out)
for e in a.em:
    nm, d = e.split('=', 1); out = {}
    for src, f in (('posterior', d + '/post_master.txt'), ('motion', d + '/motion_coord.txt')):
        if not os.path.exists(f): continue
        p = np.loadtxt(f); ok = use & (p >= 0)
        for j, vn in enumerate(names):
            m = ok & np.isfinite(X[:, j]); v = X[m, j]
            if len(np.unique(v)) > 1: out.setdefault(vn, {})[src] = one(p[m], v, binv[j])
    res[nm] = dict(kind='em', vars={vn: dict(lin=max(q.values()), **q) for vn, q in out.items()})
for e in a.labels:
    nm, f = e.split('=', 1)
    if not os.path.exists(f): continue
    lab = np.loadtxt(f, dtype=int); ok = use & (lab > 0); out = {}
    for j, vn in enumerate(names):
        m = ok & np.isfinite(X[:, j]); v = X[m, j]
        if len(np.unique(v)) > 1: out[vn] = dict(lin=eta2(lab[m], v), classes=int(len(np.unique(lab[m]))))
    res[nm] = dict(kind='labels', vars=out)
print('per variable: flex = 2-fold (eo) cross-validated linear readout of the truth from the coordinates (%s), em = the posterior / motion coordinate alone, labels = eta2 of the classes'
      % ('AUC' if all(binv) else 'r^2' if not any(binv) else 'AUC for binary, r^2 for continuous'))
print('%-22s %-6s ' % ('', 'kind') + ' '.join('%7s' % n[:7] for n in names))
for nm, r in res.items(): print('%-22s %-6s ' % (nm[:22], r['kind']) + ' '.join('%7.3f' % r['vars'][n]['lin'] if n in r['vars'] else '      -' for n in names))
print('\nflex sites: own variable group (best linear readout), quadratic readout, best single coordinate, eta2 of the flex states; leakage = best variable of another group')
for nm, r in res.items():
    if r['kind'] not in ('site', 'ref') or not r['vars']: continue
    V = r['vars']; best = max(V, key=lambda n: V[n]['lin'] - (0.5 if binv[names.index(n)] else 0)); g = grp[names.index(best)]
    oth = [n for n in V if grp[names.index(n)] != g]; lk = max(oth, key=lambda n: V[n]['lin'] - (0.5 if binv[names.index(n)] else 0)) if oth else None
    print('%-22s own %s: %s | leakage %s' % (nm, g, '; '.join('%s lin %.3f quad %.3f single %.3f states %.3f' % (n, V[n]['lin'], V[n]['quad'], V[n]['single'], V[n]['states']) for n in V if grp[names.index(n)] == g),
                                          '%s %.3f' % (lk, V[lk]['lin']) if lk else '-'))
if a.json: json.dump(res, open(a.json, 'w'), indent=1)
