#!/usr/bin/env python3
"""usage: cond_resid.py TRUTH TREE_DIR FFLEX_EMB MOTION:CAND_FILE:KEPT_SITE:PROMOTED_FILE ..
Two tests for weak motions, on genuine particles; angles standardised over all genuine particles.
(a) Conditional: the best candidate posterior (CAND_FILE, relative to TREE_DIR) binned into 5 genuine quantiles within its own node level.
  - Per bin: rows, genuine, junk share, mean / sd of the angle.
  - Bin Spearman and largest step.
  - Span ratio against the oracle WITHIN the level and against the population oracle at the same n.
  - Coverage: per node, the best abs(Spearman) of any candidate posterior (node_*/site, node_*/fsite_*) with the angle inside that
    node's level. Nodes >= 0.5; the highest such nodes (no qualifying ancestor) are disjoint branches; their union covers some share of
    all genuine particles.
(b) Residualised: the population-wide score (the kept site KEPT_SITE's posterior and the promoted PROMOTED_FILE, each as a logit) minus
  its least-squares fit, on genuine particles, by:
  - the logits of the population-wide kept sites whose best angle is another angle at abs(rho) >= 0.5;
  - also with the root fflex z components added.
  Each score is scored with population bins: abs(rho), bin Spearman, span ratio vs the population oracle, largest step, 10k tails,
  junk per bin."""
import sys, os, glob, numpy as np
TF, D, FE = sys.argv[1:4]; ENT = [a.split(':') for a in sys.argv[4:]]
hd = [l.split()[1:] for l in open(TF) if l.startswith('#')][0]; T = np.genfromtxt(TF, comments='#'); tr = T[:, 0].astype(int) - 1; NM = tr.max() + 1
ang = hd[3:]; X = np.full((NM, len(ang)), np.nan); X[tr] = T[:, 3:3 + len(ang)]
junk = np.zeros(NM, bool); junk[tr] = T[:, hd.index('junk')] > 0; G = ~junk & ~np.isnan(X[:, 0]); Z = np.nan_to_num((X - X[G].mean(0)) / X[G].std(0))
def ocum(v): return np.cumsum(np.sort(v)[::-1]) / np.arange(1, len(v) + 1)
EH = [ocum(Z[G, c]) for c in range(len(ang))]; EL = [ocum(-Z[G, c]) for c in range(len(ang))]
def rank(v): r = np.empty(len(v)); r[np.argsort(v, kind='stable')] = np.arange(len(v)); return r
def spear(x, y): return float(np.corrcoef(rank(np.asarray(x, float)), rank(np.asarray(y, float)))[0, 1])
def readpost(f):
    Q = np.loadtxt(f, ndmin=2); Q = Q.T if Q.shape[1] > Q.shape[0] else Q; return Q[:, 0].copy()
def logit(p): p = np.clip(p, 1e-6, 1 - 1e-6); return np.log(p / (1 - p))
def bins(sc, lev, c, orc=None):
    gl = lev & G; s = sc + np.random.default_rng(7).uniform(0., 1e-9 * max(1., np.abs(sc[gl]).max()), NM); sg = np.sort(s[gl]); n = len(sg)
    cuts = [sg[int(n * q) - 1] for q in (0.2, 0.4, 0.6, 0.8)]; b = np.searchsorted(np.array(cuts), s, side='left'); r = spear(sc[gl], Z[gl, c]); sgn = 1. if r > 0 else -1.
    mu, sd, ng, nr, js = [], [], [], [], []
    for k in range(5):
        m = lev & (b == k); g = m & G; mu.append(sgn * Z[g, c].mean()); sd.append(Z[g, c].std()); ng.append(g.sum()); nr.append(m.sum()); js.append((m & junk).sum() / max(m.sum(), 1))
    mu = np.array(mu); span = mu[-1] - mu[0]; hiE, loE = (EH[c], EL[c]) if sgn > 0 else (EL[c], EH[c])
    pop = span / (hiE[ng[-1] - 1] + loE[ng[0] - 1]); out = dict(rho=abs(r), bs=spear(np.arange(5), mu), span=span, pop=pop, step=np.diff(mu).max(), mu=mu, sd=sd, ng=ng, nr=nr, js=js)
    if orc is not None: lh, ll = (ocum(sgn * Z[gl, c]), ocum(-sgn * Z[gl, c])); out['lev'] = span / (lh[ng[-1] - 1] + ll[ng[0] - 1])
    o = np.where(lev)[0][np.argsort(s[lev], kind='stable')]; t = []
    for end, rows in (('high', o[-10000:]), ('low', o[:10000])):
        g = rows[G[rows]]; e = (1 if end == 'high' else -1) * sgn * Z[g, c].mean(); w = 'high' if (end == 'high') == (sgn > 0) else 'low'
        t.append('%.2f' % (e / (EH[c] if w == 'high' else EL[c])[len(g) - 1]))
    out['tails'] = '/'.join(t); return out
def fmt(o): return ('abs rho %.3f | bins Spearman %.2f, span %.2f, ratio vs population %.2f%s, max step %.2f | 10k tails %s | means %s | junk %s'
                    % (o['rho'], o['bs'], o['span'], o['pop'], ', vs level %.2f' % o['lev'] if 'lev' in o else '', o['step'], o['tails'],
                       ' '.join('%+.2f' % x for x in o['mu']), ' '.join('%.2f' % x for x in o['js'])))
def node_of(f): return f.split('/node_')[1].split('/')[0]
cand = sorted(glob.glob(D + '/node_*/site/site_post_*.txt') + glob.glob(D + '/node_*/fsite_*/site_post_*.txt')); CR = {}
for f in cand:
    try: sc = np.array(open(f).read().split(), float)
    except Exception: continue
    if len(sc) != NM: continue
    gl = (sc >= 0) & G
    if gl.sum() < 1000: continue
    rs = rank(sc[gl]); CR[f] = (node_of(f), gl, {c: abs(float(np.corrcoef(rs, rank(Z[gl, c]))[0, 1])) for c in range(len(ang))})
print('read %d candidate posteriors' % len(CR), flush=True)
kept = []
for w in [l.split() for l in open(D + '/factorial/table.txt') if l.strip() and not l.startswith('#')]:
    p = readpost(w[3]); gl = (p >= 0) & G; rh = [abs(spear(p[gl], Z[gl, c])) for c in range(len(ang))]; kept.append((w[0], p, gl.sum() == G.sum(), rh))
b = open(FE, 'rb').read(); ver, n, nc, box = (int(v) for v in np.frombuffer(b[8:24], '<i4')); pi = np.frombuffer(b[28:28 + 4 * n], '<i4').astype(int) - 1
zf = np.zeros((NM, nc)); zf[pi] = np.frombuffer(b[28 + 4 * n:28 + 4 * n + 8 * n * nc], '<f8').reshape(nc, n).T; zf = (zf - zf[G].mean(0)) / zf[G].std(0)
for mot, cf, ks, pf in ENT:
    c = ang.index(mot); print('\n######## %s' % mot)
    sc = readpost(D + '/' + cf); lev = sc >= 0; o = bins(sc, lev, c, orc=True)
    print('(a) conditional, %s (node %s): level %d rows, %d genuine' % (cf, node_of(D + '/' + cf), lev.sum(), (lev & G).sum()))
    print('    bins: genuine %s, rows %s, sd %s' % (' '.join('%d' % x for x in o['ng']), ' '.join('%d' % x for x in o['nr']), ' '.join('%.2f' % x for x in o['sd'])))
    print('    ' + fmt(o))
    per = {}
    for f, (nd, gl, rr) in CR.items():
        if rr[c] > per.get(nd, (0., None, ''))[0]: per[nd] = (rr[c], gl, f.split('/1_hetero_tree/')[1])
    q = {nd: v for nd, v in per.items() if v[0] >= 0.5}
    top = [nd for nd in q if not any(nd != a and nd.startswith(a) for a in q)]; cov = np.any([q[nd][1] for nd in top], axis=0) if top else np.zeros(NM, bool)
    print('    coverage: %d nodes with a candidate >= 0.5 on %s; %d disjoint branches (highest such nodes) cover %.3f of genuine: %s'
          % (len(q), mot, len(top), cov.sum() / G.sum(), ', '.join('%s %.2f (%d)' % (nd, q[nd][0], q[nd][1].sum()) for nd in sorted(top, key=lambda x: -q[x][1].sum()))))
    regs = [(nm, logit(p)) for nm, p, full, rh in kept if full and int(np.argmax(rh)) != c and max(rh) >= 0.5]
    print('(b) residualised; regressors = %d population-wide kept sites on other angles: %s' % (len(regs), ' '.join('%s(%s)' % (nm, ang[int(np.argmax(rh))]) for nm, p, full, rh in kept if full and int(np.argmax(rh)) != c and max(rh) >= 0.5)))
    for lab, y in (('kept site %s' % ks, logit(readpost(D + '/factorial/post/%s.txt' % ks))), ('promoted %s' % pf.split('/')[-2] + '/' + pf.split('/')[-1], logit(readpost(pf)))):
        print('    %-40s original : %s' % (lab, fmt(bins(y, np.ones(NM, bool), c))))
        for rl, R in ([('- other kept sites', np.column_stack([r for _, r in regs]))] if regs else []) + [('- other kept sites - fflex z', np.column_stack([r for _, r in regs] + [zf]))]:
            A = np.column_stack([np.ones(NM), R]); beta = np.linalg.lstsq(A[G], y[G], rcond=None)[0]; res = y - A @ beta
            r2 = 1 - ((y[G] - A[G] @ beta) ** 2).sum() / ((y[G] - y[G].mean()) ** 2).sum()
            print('    %-40s %-28s (R^2 %.2f): %s' % ('', rl, r2, fmt(bins(res, np.ones(NM, bool), c))))
