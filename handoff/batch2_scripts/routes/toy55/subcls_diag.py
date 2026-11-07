#!/usr/bin/env python3
"""What is the reproducible sub-class structure made of? (diagnostic of a subcls_embed3 run). Per sub-class: the pooled latent of
its even and odd members (raw / bin-centred) and attribute means over all its members: composition (truth: conformer shares, swing
share and mean, junk share), view offset inside the bin, in-plane angle (cos, sin of 2 e3 and of e3), shift x, y, |s|, defocus,
astigmatism (dfx - dfy, cos/sin 2 angast), the flex residual and mean residual, member count. Cross-half R2: a ridge map from an
attribute set to the EVEN-half means, scored on the ODD-half means (and vice versa), over sub-classes; = share of the reproducible
sub-class variance the set explains. Also the per-component split-half correlation and which attribute each top reliability
direction follows.
usage: subcls_diag.py EMBED_DIR ORIS NROWS GT [--codes swing:200:230] [--mode bincentred]"""
import argparse, os, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('d'); ap.add_argument('oris'); ap.add_argument('nrows', type=int); ap.add_argument('gt')
ap.add_argument('--codes', default='swing:200:230'); ap.add_argument('--mode', default='bincentred'); ap.add_argument('--minh', type=int, default=5)
a = ap.parse_args(); E = np.load(os.path.join(a.d, 'subcls_embed.npz')); rng = np.random.default_rng(0)
p, eo, sid, nS, sbin, Z, res, mres = E['p'], E['eo'], E['sid'], E['nS'], E['sbin'], E['Z'], E['res'], E['mres']; S = len(nS); N = len(p); binn = sbin[sid]
keys = ['e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast']; A = np.zeros((a.nrows, len(keys))); i = 0
for ln in open(a.oris):
    if 'e1=' not in ln: continue
    kv = dict(t.split('=', 1) for t in ln.split() if '=' in t); A[i] = [float(kv.get(k, 0)) for k in keys]; i += 1
A = A[p - 1]; e1, e2, e3, x, y, dfx, dfy, ang = A.T
G = np.loadtxt(a.gt, comments='#', ndmin=2); conf_of = np.zeros(a.nrows + 1, int); conf_of[G[:, 0].astype(int)] = G[:, 1].astype(int); conf = conf_of[p]
def gm(X, g=sid, G=S, w=None):
    X = X.reshape(len(g), -1); w = np.ones(len(g)) if w is None else w; c = np.bincount(g, w, minlength=G)
    return np.stack([np.bincount(g, w * X[:, j], minlength=G) for j in range(X.shape[1])], 1) / np.maximum(c, 1e-12)[:, None]
if a.mode == 'bincentred': ub, bi = np.unique(binn, return_inverse=True); Z = Z - gm(Z, bi, len(ub))[bi]
Me = gm(Z, w=(eo == 0).astype(float)); Mo = gm(Z, w=(eo == 1).astype(float)); ne = np.bincount(sid, eo == 0, minlength=S); no = np.bincount(sid, eo == 1, minlength=S)
ok = (ne >= a.minh) & (no >= a.minh)
# attribute sets (sub-class means over all members)
v = np.c_[np.sin(np.radians(e2)) * np.cos(np.radians(e1)), np.sin(np.radians(e2)) * np.sin(np.radians(e1)), np.cos(np.radians(e2))]
ub, bi = np.unique(binn, return_inverse=True); cen = gm(v, bi, len(ub)); cen /= np.linalg.norm(cen, axis=1, keepdims=True); c = cen[bi]
t1 = np.cross(c, [0, 0, 1.]); t1 /= np.maximum(np.linalg.norm(t1, axis=1, keepdims=True), 1e-9); t2 = np.cross(c, t1)
off = np.c_[np.sum((v - c) * t1, 1), np.sum((v - c) * t2, 1)]
r3 = np.radians(e3); ra = np.radians(ang)
sets = {'view offset in bin': off, 'in-plane e3': np.c_[np.cos(r3), np.sin(r3), np.cos(2 * r3), np.sin(2 * r3)], 'shift': np.c_[x, y, np.hypot(x, y)],
        'defocus+astig': np.c_[(dfx + dfy) / 2, dfx - dfy, np.cos(2 * ra), np.sin(2 * ra)], 'flex residuals': np.c_[res, mres, res / np.maximum(mres, 1e-12)],
        'members': np.log(np.bincount(sid, minlength=S))[sid]}
uc, cc = np.unique(conf, return_counts=True); top = uc[np.argsort(-cc)][:12]
comp = np.stack([(conf == t).astype(float) for t in top], 1); lo, hi = [int(t) for t in a.codes.split(':')[1:]]; isw = (conf >= lo) & (conf <= hi)
sets['TRUTH composition'] = np.c_[comp, isw, np.where(isw, conf - lo, 0.), conf == 0]
SM = {k: gm(X) for k, X in sets.items()}; SM['ALL nuisance'] = np.concatenate([SM[k] for k in sets if not k.startswith('TRUTH')], 1); SM['nuisance + composition'] = np.concatenate([SM[k] for k in SM if k != 'ALL nuisance'], 1)
def xr2(X, alpha=1e-2):  # cross-half R2 (fit even -> score odd and vice versa), over ok sub-classes, all latent dims together
    X = X[ok]; X = (X - X.mean(0)) / (X.std(0) + 1e-12); out = []
    for T_, V_ in ((Me[ok], Mo[ok]), (Mo[ok], Me[ok])):
        idx = rng.permutation(len(X)); f = idx % 5; pred = np.zeros_like(V_)
        for k in range(5):
            tr, te = f != k, f == k; mu = T_[tr].mean(0); mx = X[tr].mean(0); Xa = X[tr] - mx
            W = np.linalg.solve(Xa.T @ Xa + alpha * len(Xa) * np.eye(X.shape[1]), Xa.T @ (T_[tr] - mu)); pred[te] = (X[te] - mx) @ W + mu
        out.append(1 - ((V_ - pred) ** 2).sum() / ((V_ - V_.mean(0)) ** 2).sum())
    return np.mean(out)
rep = 1 - ((Me[ok] - Mo[ok]) ** 2).sum() / (((Me[ok] - Me[ok].mean(0)) ** 2).sum() + ((Mo[ok] - Mo[ok].mean(0)) ** 2).sum()) * 2 / 2
cx = [np.corrcoef(Me[ok, j], Mo[ok, j])[0, 1] for j in range(Z.shape[1])]
print('[%s] %d sub-classes; per latent component split-half corr of the sub-class means: %s' % (a.mode, ok.sum(), ' '.join('%.2f' % c_ for c_ in cx)))
print('cross-half R2 of the sub-class mean latent (all %d dims) explained by each attribute set (ceiling = split-half reliability of the means):' % Z.shape[1])
for k in SM: print('  %-24s %.3f' % (k, xr2(SM[k])))
for j in range(min(4, Z.shape[1])):
    yv = (Me[ok, j] + Mo[ok, j]) / 2; line = '  latent z%d: |corr| of its sub-class means with' % (j + 1)
    for k in [k for k in sets]:
        Xs = SM[k][ok]; best = max(abs(np.corrcoef(Xs[:, q], yv)[0, 1]) for q in range(Xs.shape[1]) if Xs[:, q].std() > 0); line += ' %s %.2f;' % (k, best)
    print(line)
