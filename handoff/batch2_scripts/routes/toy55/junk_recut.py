#!/usr/bin/env python3
"""usage: junk_recut.py SCORES [--labels L --k1 1 --truth T --tcol junk --jtype COL --guard COLS --guardtail COLS --tag NAME --out FLAGS]
Re-cut saved junk scores (junk_score.py's junk_scores.txt: amplitude a, DEN) with heteroscedastic mixtures, a_i ~ N(mu_c, s_c^2 + 1/DEN_i):
xd2 = two components, junk = the lower one (junk_score.py); xd3 = junk + a genuine core + a broad genuine tail, genuine = the heaviest
component, junk = every other component >= 2 sd below it (weight >= 0.5 %). Junk = summed junk posterior >= 0.5."""
import argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('scores'); ap.add_argument('--labels', default=''); ap.add_argument('--k1', type=int, default=0)
for k in ('truth', 'tcol', 'jtype', 'guard', 'guardtail', 'tag', 'out'): ap.add_argument('--' + k, default='junk' if k == 'tcol' else '')
a = ap.parse_args(); S = np.loadtxt(a.scores); node = np.isfinite(S[:, 1]); amp = S[node, 1]; var = 1. / np.maximum(S[node, 6], 1e-12); mv = np.median(var)
def xd(mu, s2, pi, iters=500):
    for _ in range(iters):
        T = s2[None] + var[:, None]; L = pi[None] * np.exp(-0.5 * (amp[:, None] - mu[None]) ** 2 / T) / np.sqrt(T); R = L / np.maximum(L.sum(1, keepdims=True), 1e-300)
        bm = mu[None] + s2[None] / T * (amp[:, None] - mu[None]); Bv = s2[None] * var[:, None] / T; w = np.maximum(R.sum(0), 1e-12)
        pi = R.mean(0); mu = (R * bm).sum(0) / w; s2 = np.maximum((R * ((bm - mu[None]) ** 2 + Bv)).sum(0) / w, 1e-8)
    return mu, s2, pi, R
p2, p5, p50 = np.percentile(amp, [2, 5, 50]); vv = amp.var()
L = np.loadtxt(a.labels, dtype=int)[node] if a.labels else np.ones(node.sum(), int)
if a.k1: L = (L > 0).astype(int)
if a.truth:
    hd = [l.split()[1:] for l in open(a.truth) if l.startswith('#')][0]; T = np.genfromtxt(a.truth, comments='#'); tr = T[:, 0].astype(int) - 1
    def col(cn): v = np.full(len(S), np.nan); v[tr] = T[:, hd.index(cn)]; return v[node]
    tj = col(a.tcol) > 0
out = {}
for nm in ('xd2', 'xd3'):
    if nm == 'xd2':
        mu, s2, pi, R = xd(np.array([p5, p50]), np.full(2, 0.25 * vv), np.array([0.2, 0.8])); j = int(np.argmin(mu))
        J = [j] if pi[j] >= 0.005 and abs(mu[1] - mu[0]) / np.sqrt(s2.mean() + mv) >= 2. else []
    else:
        mu, s2, pi, R = xd(np.array([p2, p50, p50]), np.array([0.1, 0.25, 1.]) * vv, np.array([0.05, 0.75, 0.2])); g = int(np.argmax(pi))
        J = [c for c in range(3) if c != g and pi[c] >= 0.005 and (mu[g] - mu[c]) / np.sqrt(0.5 * (s2[g] + s2[c]) + mv) >= 2.]
    F = R[:, J].sum(1) >= 0.5 if J else np.zeros(len(amp), bool); out[nm] = F
    print('%s %s: components (mean / intrinsic sd / weight) %s | junk components %s | flagged %d of %d (%.4f)' % (a.tag, nm, '  '.join('%.2f/%.2f/%.3f' % (m, np.sqrt(s), w) for m, s, w in sorted(zip(mu, s2, pi))),
          [round(float(mu[c]), 2) for c in J] or 'none', F.sum(), len(F), F.mean()))
    if a.labels: print('  per state: ' + ' '.join('%d: %d/%d (%.4f)' % (s, (F & (L == s)).sum(), (L == s).sum(), (F & (L == s)).sum() / max((L == s).sum(), 1)) for s in np.unique(L)))
    if a.truth:
        print('  vs truth: precision %.3f recall %.3f; false-junk on genuine %.4f' % ((F & tj).sum() / max(F.sum(), 1), (F & tj).sum() / max(tj.sum(), 1), (F & ~tj).sum() / max((~tj).sum(), 1))
              + ''.join('; recall junk type %g %.3f' % (v, (F & tj & (col(a.jtype) == v)).sum() / max((tj & (col(a.jtype) == v)).sum(), 1)) for v in (np.unique(col(a.jtype)[tj]) if a.jtype else [])))
        for gc_ in [x for x in a.guard.split(',') if x]:
            gc = col(gc_); print('  guard %s: ' % gc_ + ', '.join('= %g false-junk %.4f' % (v, (F & ~tj & (gc == v)).sum() / max((~tj & (gc == v)).sum(), 1)) for v in np.unique(gc[~tj])))
        tails = [x for x in a.guardtail.split(',') if x]
        if tails:
            ext = np.zeros(len(F), bool)
            for gc_ in tails: gc = col(gc_); lo, hi = np.nanpercentile(gc[~tj], [2.5, 97.5]); ext |= ~tj & ((gc <= lo) | (gc >= hi))
            print('  guard: genuine with any of %s in its outer 5 %% (%d): false-junk %.4f' % (','.join(tails), ext.sum(), (F & ext).sum() / max(ext.sum(), 1)))
if a.out:
    f = np.full(len(S), -1); f[node] = out['xd3']; np.savetxt(a.out, f, fmt='%d')
