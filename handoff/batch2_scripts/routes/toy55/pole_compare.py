#!/usr/bin/env python3
"""usage: pole_compare.py --rows ROWS --mean MEAN_RECONSTRUCT3D_DIR --runs NAME=POLE_DIRS_OUT .. [--lp 15 --mskdiam 232 --smpd 1.464 --arm ARM_TXT --n 2500]"""
import argparse, os, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--rows', required=True); ap.add_argument('--mean', required=True); ap.add_argument('--runs', nargs='+', required=True)
ap.add_argument('--lp', type=float, default=15.); ap.add_argument('--mskdiam', type=float, default=232.); ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--arm', default='')
ap.add_argument('--n', type=int, default=2500); a = ap.parse_args()
use = np.loadtxt(a.rows, dtype=int) > 0
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
me, mo = rd(a.mean + '/recvol_state01_even_unfil.mrc'), rd(a.mean + '/recvol_state01_odd_unfil.mrc'); n = me.shape[0]
f = np.fft.fftfreq(n, d=a.smpd); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); G = np.exp(-0.5 * (kx ** 2 + ky ** 2 + kz ** 2) * (a.lp / 1.2) ** 2)
lp = lambda v: np.real(np.fft.ifftn(np.fft.fftn(v) * G)); c = np.arange(n) - n // 2; X, Y, W = np.meshgrid(c, c, c, indexing='ij')
mf = lp(0.5 * (me + mo)); sup = (mf > 0.2 * mf.max()) & (X ** 2 + Y ** 2 + W ** 2 <= (a.mskdiam / 2 / a.smpd) ** 2); me, mo = lp(me), lp(mo)
def dev(e, o):
    out = []
    for x, m in ((lp(e), me), (lp(o), mo)): s = (x[sup] * m[sup]).sum() / (x[sup] ** 2).sum(); out.append(s * x - m)
    return out[0], out[1], 0.5 * (out[0] + out[1])
cr = lambda x, y: float(np.corrcoef(x[sup], y[sup])[0, 1])
def rel(de, do): r = cr(de, do); return max(2 * r / (1 + r), 1e-6)
arm = {}
if a.arm and os.path.exists(a.arm): arm = {os.path.basename(w[0])[:-4]: (float(w[1]), float(w[2])) for w in (l.split() for l in open(a.arm)) if w and not w[0].startswith('#')}
P = []
for e in a.runs:
    run, D = e.split('=', 1)
    if not os.path.exists(D + '/poles.txt'): print('%s: no poles.txt' % run); continue
    sp = {w[0]: (w[1], w[2]) for w in (l.split() for l in open(D + '/spec_poles.txt')) if len(w) == 3} if os.path.exists(D + '/spec_poles.txt') else {}
    for w in (l.split() for l in open(D + '/poles.txt')):
        if not w or w[0].startswith('#'): continue
        p = np.loadtxt(w[5]); mp = sp.get('%s_top%d' % (w[0], a.n))
        dv = dev(*[rd(x) for x in mp]) if mp else None
        P.append(dict(run=run, name=w[0], d=int(w[1]), side=w[2], nstar=int(w[3]), hold=float(w[4]), post=p, dev=dv, rel=rel(dv[0], dv[1]) if dv else np.nan, rd=arm.get('%s_%s_top%d' % (run, w[0], a.n))))
runs = sorted(set(q['run'] for q in P)); print('poles per run (directions kept; pole: N*, hold-out cross-half corr, top-%d readout closed / open, half-map reliability of its deviation):' % a.n)
for r in runs:
    Q = [q for q in P if q['run'] == r]
    print('  %-8s %d directions (%s), %d poles: %s' % (r, len(set(q['d'] for q in Q)), ' '.join(map(str, sorted(set(q['d'] for q in Q)))), len(Q),
          '; '.join('%s N* %d hold %.2f %s rel %.2f' % (q['name'], q['nstar'], q['hold'], '%.2f/%.2f' % q['rd'] if q['rd'] else '-', q['rel']) for q in Q)))
nm = ['%s:%s' % (q['run'], q['name']) for q in P]; m = use & np.all([q['post'] >= 0 for q in P], 0)
Rk = [np.argsort(np.argsort(-q['post'][m])) for q in P]; S = np.array([[np.corrcoef(x, y)[0, 1] for y in Rk] for x in Rk])
def groups(th):
    g = list(range(len(P)))
    for i in range(len(P)):
        for j in range(i + 1, len(P)):
            if S[i, j] >= th: gi, gj = g[i], g[j]; g = [gi if x == gj else x for x in g]
    return [[i for i in range(len(P)) if g[i] == G] for G in sorted(set(g), key=lambda G: -g.count(G))]
for th in (0.99, 0.95, 0.9):
    print('\ngroups at Spearman >= %.2f:' % th)
    for G in groups(th):
        rn = sorted(set(P[i]['run'] for i in G)); ix = np.ix_(G, G)
        mr = [cr(P[i]['dev'][2], P[j]['dev'][2]) / np.sqrt(P[i]['rel'] * P[j]['rel']) for i in G for j in G if i < j and P[i]['run'] != P[j]['run'] and P[i]['dev'] and P[j]['dev']]
        print('  %d poles from %d / %d runs (min Spearman %.3f): %s | N* %s | map match across runs (fraction of the half-map noise ceiling) %s'
              % (len(G), len(rn), len(runs), S[ix].min(), ' '.join(nm[i] for i in G), ' '.join(str(P[i]['nstar']) for i in G), '%.2f-%.2f' % (min(mr), max(mr)) if mr else '-'))
print('\nbest cross-run match of every pole (Spearman, map match as a fraction of the noise ceiling):')
for i, q in enumerate(P):
    best = []
    for r in runs:
        if r == q['run']: continue
        J = [j for j, x in enumerate(P) if x['run'] == r]
        if not J: continue
        j = max(J, key=lambda j: S[i, j]); mm = cr(q['dev'][2], P[j]['dev'][2]) / np.sqrt(q['rel'] * P[j]['rel']) if q['dev'] and P[j]['dev'] else np.nan
        best.append('%s:%s %.3f / %.2f' % (r, P[j]['name'], S[i, j], mm))
    print('  %-16s %s' % (nm[i], ' | '.join(best)))
