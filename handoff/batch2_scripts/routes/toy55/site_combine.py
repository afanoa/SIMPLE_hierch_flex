#!/usr/bin/env python3
"""Combine per-site classes (site_size.py: top_star.txt = class-2 membership per master row, post_master.txt) into joint labels:
c_i = (c_i1 .. c_iS); independence of every site pair = mutual information of the two binary labels against NPERM row permutations
(and the same on the posteriors' rank correlation); populations of the 2^S combinations; label files for SIMPLE (OUT/combo_labels.txt:
combination id 1 + sum_s c_is 2^(s-1) for combinations >= MINC rows, 0 elsewhere; OUT/site_s_labels.txt: 1 = class 1, 2 = class 2 for
the per-site local maps). With --truth: per site the true state of class 2 (by overlap), precision / recall / NMI against the true site,
and the combination coverage: true combinations >= COMBO_MIN rows that own a joint class >= PURE pure.
usage: site_combine.py --spec SPEC --nmaster N --out DIR [--rows ROWS --minc 1000 --truth sites_truth.txt]  SPEC lines: name top_star post_master [truth_column]"""
import os, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--spec', required=True); ap.add_argument('--nmaster', type=int, required=True); ap.add_argument('--out', required=True); ap.add_argument('--rows', default='')
ap.add_argument('--minc', type=int, default=1000); ap.add_argument('--truth', default=''); ap.add_argument('--nperm', type=int, default=200); ap.add_argument('--combo_min', type=int, default=5000)
ap.add_argument('--pure', type=float, default=0.7); ap.add_argument('--seed', type=int, default=0); a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); rng = np.random.default_rng(a.seed)
S = [l.split() for l in open(a.spec) if l.strip() and not l.startswith('#')]; names = [s[0] for s in S]
L = np.array([np.loadtxt(s[1], dtype=int) > 0 for s in S]); P = np.array([np.loadtxt(s[2]) for s in S]); inn = (P >= 0).all(0)
if a.rows: inn &= np.loadtxt(a.rows, dtype=int) > 0
print('%d sites over %d rows (in every site EM): class-2 sizes %s' % (len(S), inn.sum(), ' '.join('%s %d' % (nm, (L[j] & inn).sum()) for j, nm in enumerate(names))))
def mi(x, y):
    J = np.zeros((2, 2)); np.add.at(J, (x.astype(int), y.astype(int)), 1); J /= J.sum(); px, py = J.sum(1), J.sum(0); nz = J > 0; return float((J[nz] * np.log(J[nz] / np.outer(px, py)[nz])).sum())
def rankc(x, y):
    rx = np.argsort(np.argsort(x)); ry = np.argsort(np.argsort(y)); return float(np.corrcoef(rx, ry)[0, 1])
for i in range(len(S)):
    for j in range(i + 1, len(S)):
        x, y = L[i][inn], L[j][inn]; m0 = mi(x, y); null = np.array([mi(x, rng.permutation(y)) for _ in range(a.nperm)])
        rc = rankc(P[i][inn], P[j][inn]); rn = np.array([rankc(P[i][inn], rng.permutation(P[j][inn])) for _ in range(min(a.nperm, 50))])
        exp_ = x.mean() * y.mean(); obs = (x & y).mean()
        print('independence %s x %s: MI %.2e nats (null mean %.2e, 99th %.2e, z %.1f) | both class 2: %.0f observed vs %.0f expected | posterior rank corr %+.4f (null sd %.4f)' %
              (names[i], names[j], m0, null.mean(), np.percentile(null, 99), (m0 - null.mean()) / max(null.std(), 1e-30), obs * inn.sum(), exp_ * inn.sum(), rc, rn.std()))
cid = np.zeros(a.nmaster, int); cid[inn] = 1 + sum((L[j][inn].astype(int) << j) for j in range(len(S))); cnt = np.bincount(cid, minlength=2 ** len(S) + 1)
keep = [k for k in range(1, 2 ** len(S) + 1) if cnt[k] >= a.minc]; lab = np.where(np.isin(cid, keep), cid, 0)
remap = {k: i + 1 for i, k in enumerate(keep)}; np.savetxt(a.out + '/combo_labels.txt', np.array([remap.get(v, 0) for v in lab]), fmt='%d')
with open(a.out + '/combos.txt', 'w') as f:
    f.write('# state combination_id rows %s (1 = class 2)\n' % ' '.join(names))
    for k in keep: f.write('%d %d %d %s\n' % (remap[k], k, cnt[k], ' '.join(str(((k - 1) >> j) & 1) for j in range(len(S)))))
print('combinations (>= %d rows -> reconstructed): %s' % (a.minc, '; '.join('%s: %d%s' % ('/'.join(str(((k - 1) >> j) & 1) for j in range(len(S))), cnt[k], '' if k in keep else ' (too few)') for k in range(1, 2 ** len(S) + 1) if cnt[k] > 0)))
for j, nm in enumerate(names): np.savetxt('%s/site_%s_labels.txt' % (a.out, nm), np.where(inn, 1 + L[j].astype(int), 0), fmt='%d')
if a.truth:
    T = np.genfromtxt(a.truth, comments='#'); tn = [l.split()[1:] for l in open(a.truth) if l.startswith('#')][0][3:]; X = T[:, 3:]; combo = T[:, 1].astype(int); ok = inn & (T[:, 2] == 0)
    def nmi(x, y):
        _, xi = np.unique(x, return_inverse=True); _, yi = np.unique(y, return_inverse=True); J = np.zeros((xi.max() + 1, yi.max() + 1)); np.add.at(J, (xi, yi), 1); J /= J.sum()
        px, py = J.sum(1), J.sum(0); nz = J > 0; I = (J[nz] * np.log(J[nz] / np.outer(px, py)[nz])).sum(); hx, hy = -(px * np.log(px)).sum(), -(py * np.log(py)).sum(); return I / max(np.sqrt(hx * hy), 1e-12)
    print('truth: %d rows in the site EMs are true junk (excluded from the scores)' % (inn & (T[:, 2] > 0)).sum())
    for j, s in enumerate(S):
        cols = [int(s[3])] if len(s) > 3 else list(range(X.shape[1]))
        best = max(((nmi(L[j][ok], X[ok, t] > 0), t) for t in cols)); t = best[1]; tv = X[ok, t] > 0; cl = L[j][ok]
        st = 1 if (cl & tv).sum() >= (cl & ~tv).sum() else 0; tt = tv if st == 1 else ~tv
        print('site %s ~ true %s state %d: class 2 %d rows, precision %.3f, recall %.3f (true count %d), NMI %.3f' % (s[0], tn[t], st, cl.sum(), tt[cl].mean(), (tt & cl).sum() / max(tt.sum(), 1), tt.sum(), best[0]))
    big = [k for k in keep if (cid[ok] == k).sum() >= 1]; cov = []
    for c_ in np.unique(combo[ok]):
        m = ok & (combo == c_); pur = max([((cid == k) & m).sum() / max(((cid == k) & ok).sum(), 1) for k in big] or [0.]); kk = max(big, key=lambda k: ((cid == k) & m).sum() / max(((cid == k) & ok).sum(), 1)) if big else 0
        print('true combination %d (%s): %6d rows, best joint class %s purity %.2f' % (c_, ' '.join('%s=%d' % (tn[t], X[np.where(combo == c_)[0][0], t]) for t in range(X.shape[1])), m.sum(), '/'.join(str(((kk - 1) >> j) & 1) for j in range(len(S))) if kk else '-', pur))
        if m.sum() >= a.combo_min: cov.append(pur >= a.pure)
    print('COMBINATION COVERAGE: %d of %d true combinations with >= %d rows own a joint class >= %.0f %% pure' % (sum(cov), len(cov), a.combo_min, 100 * a.pure))
