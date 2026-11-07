#!/usr/bin/env python3
"""usage: census.py --sim H --poses PA[,PB] --cache C --rows ROWS --mask M --mean MAP --coords flex_pca_coordinates.txt --truth FILE --tcols COL[,COL] --out DIR [--mra --mrb | --full 1]"""
import argparse, os, subprocess, sys, threading, time, numpy as np
from concurrent.futures import ThreadPoolExecutor
ap = argparse.ArgumentParser()
for k, d in (('sim', None), ('poses', None), ('cache', None), ('rows', None), ('mask', None), ('mean', None), ('coords', None), ('truth', None), ('tcols', None), ('out', None),
             ('mra', ''), ('mrb', ''), ('full', '0'), ('ctfabs', '0'), ('pars', '4'), ('iters', '10'), ('em_lim', '15'), ('nt', '15000'), ('ndir', '6'), ('nmid', '3'), ('nrand', '3'),
             ('small', '2500'), ('same', '0.99'), ('loose', '0.95'), ('disjoint', '0.1'), ('hold', '0.2'), ('topn', '2500')):
    ap.add_argument('--' + k, default=d, required=d is None)
a = ap.parse_args(); T = os.path.dirname(os.path.abspath(__file__)); PY = sys.executable; O = a.out; os.makedirs(O + '/starts', exist_ok=True); LOG = open(O + '/census.log', 'a'); _ll = threading.Lock()
def log(*s):
    m = '[%s census] %s' % (time.strftime('%F %T'), ' '.join(str(x) for x in s))
    with _ll: print(m, flush=True); LOG.write(m + '\n'); LOG.flush()
def run(cmd, out): return subprocess.run('%s > %s 2>&1' % (cmd, out), shell=True).returncode == 0
use = np.loadtxt(a.rows, dtype=int) > 0; NM = len(use); FULL = a.full == '1'
if FULL:
    MR = O + '/identity_rows.txt'
    if not os.path.exists(MR): np.savetxt(MR, np.arange(1, NM + 1), fmt='%d')
    MRA = MRB = MR; SA = SB = a.sim; PA = PB = a.poses; HALF = np.arange(NM) % 2
else:
    MRA, MRB = a.mra or a.sim + '_A/meta_rows.txt', a.mrb or a.sim + '_B/meta_rows.txt'; SA, SB = a.sim + '_A', a.sim + '_B'; PA, PB = a.poses.split(',')
C = np.loadtxt(a.coords, comments='#'); r0 = C[:, 0].astype(int) - 1; k = use[r0]; r0, Z = r0[k], C[k, 5:]; n = len(r0); NT = int(min(int(a.nt), 0.15 * n))
ST = {}
def add(nm, ix): s = np.zeros(NM, bool); s[r0[ix]] = True; ST[nm] = s
for d in range(min(int(a.ndir), Z.shape[1])):
    o = np.argsort(Z[:, d]); add('z%dlo' % (d + 1), o[:NT]); add('z%dhi' % (d + 1), o[::-1][:NT])
for d in range(min(int(a.nmid), Z.shape[1])):
    o = np.argsort(Z[:, d]); m = n // 2; add('z%dmid' % (d + 1), o[m - NT // 2:m + NT // 2])
perm = np.random.default_rng(7).permutation(n)
for i in range(int(a.nrand)): add('rand%d' % (i + 1), perm[i * NT:(i + 1) * NT])
ns = int(min(int(a.small), NT // 2)); o = np.argsort(Z[:, 0]); add('z1lo_small', o[:ns]); add('z1hi_small', o[::-1][:ns])
log('%d rows, latent %s (%d dims), %d starts of %d (small %d): %s' % (n, a.coords, Z.shape[1], len(ST), NT, ns, ' '.join(ST)))
_es = threading.Semaphore(int(a.pars))
def em(nm):
    D = '%s/%s' % (O, nm); os.makedirs(D, exist_ok=True); s = ST[nm]
    if not os.path.exists(D + '/curve.txt'):
        for h, mr in (('A', MRA), ('B', MRB)):
            rr = np.loadtxt(mr).astype(int) - 1; ok = use[rr] & ((HALF[rr] == (0 if h == 'A' else 1)) if FULL else True)
            np.savetxt('%s/start_%s.txt' % (D, h), np.where(ok, 1 + (s & use)[rr].astype(int), 0), fmt='%d')
        with _es:
            ok = run('%s -u %s/minor_em3.py --sim %s,%s --poses %s,%s --labels %s/start_A.txt,%s/start_B.txt --cache %s --ctfabs %s --cons_ref %s --mask %s --soft 1 --iters %s --em_lim %s --post 1 --nomaps 1 --out %s/em'
                     % (PY, T, SA, SB, PA, PB, D, D, a.cache, a.ctfabs, a.mean, a.mask, a.iters, a.em_lim, D), D + '/em.log') and \
                 run('%s -u %s/site_size.py --em %s/em --start %s/start_A.txt,%s/start_B.txt --sim %s,%s --poses %s,%s --mr %s,%s --nmaster %d --mask %s --out %s --cache %s --ctfabs %s --em_lim %s'
                     % (PY, T, D, D, D, SA, SB, PA, PB, MRA, MRB, NM, a.mask, D, a.cache, a.ctfabs, a.em_lim), D + '/size.log')
        if not ok: log('EM %s FAILED' % nm); return None
    cv = np.loadtxt(D + '/curve.txt', comments='#', ndmin=2); nst = int(open(D + '/nstar.txt').read()); hold = float(cv[np.argmin(np.abs(cv[:, 0] - nst)), 2])
    return dict(post=np.loadtxt(D + '/post_master.txt'), nstar=nst, hold=hold)
with ThreadPoolExecutor(len(ST)) as ex: E = {k_: v for k_, v in zip(ST, ex.map(em, list(ST))) if v}
nm = list(E); m = use & np.all([E[x]['post'] >= 0 for x in nm], 0); R = {x: np.argsort(np.argsort(-E[x]['post'][m])) for x in nm}
S = np.array([[np.corrcoef(R[x], R[y])[0, 1] for y in nm] for x in nm])
ov = lambda x, y: (ST[x] & ST[y]).sum() / max(min(ST[x].sum(), ST[y].sum()), 1)
hd = [l.split()[1:] for l in open(a.truth) if l.startswith('#')][0]; TT = np.genfromtxt(a.truth, comments='#', invalid_raise=False, usecols=range(len([c for c in hd if not c.startswith('(')])))
tr = {c: np.full(NM, np.nan) for c in a.tcols.split(',')}
for c in tr: j = hd.index(c); tr[c][TT[:, 0].astype(int) - 1] = TT[:, j]
if 'junk' in hd: jk = np.zeros(NM, bool); jk[TT[:, 0].astype(int) - 1] = TT[:, hd.index('junk')] > 0
else: jk = np.zeros(NM, bool)
def truth_summary(x):
    p = np.where(m & ~jk, E[x]['post'], -1.); top = np.argsort(-p)[:int(a.topn)]; out = []
    for c, v in tr.items():
        ok = m & ~jk & np.isfinite(v); vv = v[top][np.isfinite(v[top])]
        if len(np.unique(v[ok])) <= 4: out.append('%s %s' % (c, ' '.join('%g:%.2f(%.2fx)' % (u, (vv == u).mean(), (vv == u).mean() / max((v[ok] == u).mean(), 1e-9)) for u in np.unique(v[ok]))))
        else: q = np.percentile(v[ok], [10, 50, 90]); out.append('%s %.1f+-%.1f (all %.1f; p10/p90 %.1f/%.1f)' % (c, vv.mean(), vv.std(), v[ok].mean(), q[0], q[2]))
    return '; '.join(out)
for thr in (float(a.same), float(a.loose)):
    g = list(range(len(nm)))
    for i in range(len(nm)):
        for j in range(i + 1, len(nm)):
            if S[i, j] >= thr: gi, gj = g[i], g[j]; g = [gi if x == gj else x for x in g]
    print('\n==== fixed points at Spearman >= %.2f: %d groups over %d endpoints (attractor = >= 2 disjoint starts [overlap <= %.2f] and hold-out >= %.2f)' % (thr, len(set(g)), len(nm), float(a.disjoint), float(a.hold)))
    for G in sorted(set(g), key=lambda G: -sum(x == G for x in g)):
        mem = [nm[i] for i in range(len(nm)) if g[i] == G]; ix = [nm.index(x) for x in mem]; rep = max(mem, key=lambda x: E[x]['hold'])
        pairs = any(ov(x, y) <= float(a.disjoint) for x in mem for y in mem if x < y); att = pairs and E[rep]['hold'] >= float(a.hold)
        print('  %s group of %d (min Spearman %.3f): %s | N* %s | hold-out %s | truth top-%s: %s' % ('ATTRACTOR' if att else 'fixed pt ', len(mem), S[np.ix_(ix, ix)].min(), ' '.join(mem),
              ' '.join(str(E[x]['nstar']) for x in mem), ' '.join('%.2f' % E[x]['hold'] for x in mem), a.topn, truth_summary(rep)))
print('\n==== start overlaps (fraction of the smaller start) for starts sharing a group at >= %.2f are in census.log' % float(a.loose))
for i in range(len(nm)):
    for j in range(i + 1, len(nm)):
        if S[i, j] >= float(a.loose): log('%s ~ %s Spearman %.3f, start overlap %.2f' % (nm[i], nm[j], S[i, j], ov(nm[i], nm[j])))
