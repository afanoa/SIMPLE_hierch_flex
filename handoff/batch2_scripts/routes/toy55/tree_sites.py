#!/usr/bin/env python3
"""usage: tree_sites.py --tree DIR --out DIR --h HALFPREFIX --cache P [--mra F --mrb F --ctfabs 1 --mskdiam 220 --nodes r,r1 --nsn 3 --npole 15000 --gate 0.8 --tmasks A,B]"""
import argparse, glob, os, subprocess, sys, threading, time, numpy as np
from concurrent.futures import ThreadPoolExecutor
ap = argparse.ArgumentParser()
for k, d in (('tree', None), ('out', None), ('h', None), ('cache', None), ('mra', ''), ('mrb', ''), ('ctfabs', '1'), ('mskdiam', '220'), ('nodes', ''), ('nsn', '3'), ('npole', '15000'),
             ('polefrac', '0.15'), ('gate', '0.8'), ('dice', '0.5'), ('maxmol', '0.3'), ('k', '8'), ('lpf', '15'), ('np', '8'), ('pars', '4'), ('parf', '2'), ('iters', '10'), ('em_lim', '15'), ('tmasks', '')):
    ap.add_argument('--' + k, default=d, required=d is None)
a = ap.parse_args(); T = os.path.dirname(os.path.abspath(__file__)); PY = sys.executable; R, O = a.tree, a.out; os.makedirs(O, exist_ok=True)
MRA = a.mra or a.h + '_A/meta_rows.txt'; MRB = a.mrb or a.h + '_B/meta_rows.txt'; GATE, DICE = float(a.gate), float(a.dice); LOG = open(O + '/tree_sites.log', 'a'); _ll = threading.Lock()
def log(*s):
    m = '[%s tree sites] %s' % (time.strftime('%F %T'), ' '.join(str(x) for x in s))
    with _ll: print(m, flush=True); LOG.write(m + '\n'); LOG.flush()
def run(cmd, out): return subprocess.run('%s > %s 2>&1' % (cmd, out), shell=True).returncode == 0
nodes = sorted([d.split('/node_')[-1] for d in glob.glob(R + '/node_*') if os.path.exists(d + '/node.done')], key=lambda u: (len(u), u))
todo = [u for u in nodes if not a.nodes or u in a.nodes.split(',')]; par = lambda u: u[:-1] if len(u) > 1 else None
_rows = {}
def rows(u):
    if u not in _rows: _rows[u] = np.loadtxt('%s/rows_%s.txt' % (R, u), dtype=int) > 0
    return _rows[u]
NM = len(rows('r'))
def mean(u): return sorted(glob.glob('%s/node_%s/mean/*_reconstruct3D' % (R, u)), key=os.path.getmtime)[-1] + '/recvol_state01.mrc'
_pl = {}
def poses(u):
    D = '%s/node_%s' % (O, u); os.makedirs(D, exist_ok=True)
    with _pl.setdefault(u, threading.Lock()):
        if not os.path.exists(D + '/poses_B.txt'):
            run('simple_exec prg=print_project_field projfile=%s/node_%s/n.simple oritype=ptcl3D mkdir=no 2>/dev/null > %s/pm.txt; for h in A:%s B:%s; do awk "NR == FNR {l[FNR] = \\$0; next} {print l[\\$1 + 0]}" %s/pm.txt ${h#*:} > %s/poses_${h%%%%:*}.tmp && mv %s/poses_${h%%%%:*}.tmp %s/poses_${h%%%%:*}.txt; done; rm -f %s/pm.txt' % (R, u, D, MRA, MRB, D, D, D, D, D), D + '/poses.log')
    return D + '/poses_A.txt', D + '/poses_B.txt'
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n)
def dice(p, q): x, y = rd(p) > 0.5, rd(q) > 0.5; return 2. * (x & y).sum() / max(x.sum() + y.sum(), 1)
def spear(x, y, m): return float(np.corrcoef(np.argsort(np.argsort(x[m])), np.argsort(np.argsort(y[m])))[0, 1])
def discover(u):
    D = '%s/node_%s/disc' % (O, u)
    if not os.path.exists(D + '/sites.txt'):
        V = '%s/node_%s/view' % (O, u); os.makedirs(V, exist_ok=True)
        if not os.path.lexists(V + '/node_' + u): os.symlink('%s/node_%s' % (R, u), V + '/node_' + u)
        run('%s -u %s/discover_sites2.py --tree %s --out %s --mskdiam %s --maxmol %s --dice %s %s' % (PY, T, V, D, a.mskdiam, a.maxmol, a.dice, '--truth ' + a.tmasks if a.tmasks else ''), '%s/node_%s/disc.log' % (O, u))
    S = [l.split() for l in open(D + '/sites.txt') if l.strip() and not l.startswith('#')] if os.path.exists(D + '/sites.txt') else []
    return [dict(name='%s_s%d' % (u, int(s[0])), rep=s[1], repro=float(s[2]), mol=float(s[4]), mask='%s/site_%02d.mrc' % (D, int(s[0]))) for s in S[:int(a.nsn)]]
_fs, _es = threading.Semaphore(int(a.parf)), threading.Semaphore(int(a.pars))
def flex(s, u):
    D = '%s/fits/%s@%s/flex' % (O, s['name'], u); os.makedirs(D, exist_ok=True); c = glob.glob(D + '/*_flex_pca/flex_pca_coordinates.txt')
    if not c:
        with _fs:
            run('cd %s && cp -f %s/node_%s/n.simple f.simple && simple_exec prg=selection projfile=f.simple oritype=ptcl3D infile=%s/rows_%s.txt mkdir=no && simple_exec prg=flex_pca projfile=f.simple vol1=%s mskdiam=%s lp=%s npreimages=%s rec_states=yes min_state_frac=0 umap=no n_probe_iters=2 pcg_mskfile=%s nparts=%s nthr=8'
                % (D, R, u, R, u, mean(u), a.mskdiam, a.lpf, a.k, s['mask'], a.np), D + '/flex.log')
        c = glob.glob(D + '/*_flex_pca/flex_pca_coordinates.txt')
        if not c: log('flex %s at %s FAILED (%s/flex.log)' % (s['name'], u, D)); return None
    return np.loadtxt(sorted(c, key=os.path.getmtime)[-1], comments='#')
def em(mask, u, seed, D):
    os.makedirs(D, exist_ok=True)
    if not os.path.exists(D + '/post_master.txt'):
        use = rows(u); pa, pb = poses(u)
        for h, mr in (('A', MRA), ('B', MRB)):
            r = np.loadtxt(mr).astype(int) - 1; np.savetxt('%s/start_%s.txt' % (D, h), np.where(use[r], 1 + (seed & use)[r].astype(int), 0), fmt='%d')
        with _es:
            ok = run('%s -u %s/minor_em3.py --sim %s_A,%s_B --poses %s,%s --labels %s/start_A.txt,%s/start_B.txt --cache %s --ctfabs %s --cons_ref %s --mask %s --soft 1 --iters %s --em_lim %s --post 1 --nomaps 1 --out %s/em'
                     % (PY, T, a.h, a.h, pa, pb, D, D, a.cache, a.ctfabs, mean(u), mask, a.iters, a.em_lim, D), D + '/em.log') and \
                 run('%s -u %s/site_size.py --em %s/em --start %s/start_A.txt,%s/start_B.txt --sim %s_A,%s_B --poses %s,%s --mr %s,%s --nmaster %d --mask %s --out %s --cache %s --ctfabs %s --em_lim %s'
                     % (PY, T, D, D, D, a.h, a.h, pa, pb, MRA, MRB, NM, mask, D, a.cache, a.ctfabs, a.em_lim), D + '/size.log')
        if not ok: log('EM %s FAILED' % D); return None
    return np.loadtxt(D + '/post_master.txt'), np.loadtxt(D + '/top_star.txt', dtype=int) > 0
def fit(s, u):
    C = flex(s, u)
    if C is None: return None
    r, z = C[:, 0].astype(int) - 1, C[:, 5]; k = int(min(int(a.npole), float(a.polefrac) * len(r))); o = np.argsort(z); lo = np.zeros(NM, bool); hi = np.zeros(NM, bool)
    lo[r[o[:k]]] = True; hi[r[o[::-1][:k]]] = True; B = '%s/fits/%s@%s' % (O, s['name'], u)
    with ThreadPoolExecutor(2) as ex: f1, f2 = ex.submit(em, s['mask'], u, lo, B + '/lo'), ex.submit(em, s['mask'], u, hi, B + '/hi'); El, Eh = f1.result(), f2.result()
    if El is None or Eh is None: return None
    log('site %s at %s: poles from both tails of the site flex z1 (%d rows each): N* %d / %d' % (s['name'], u, k, El[1].sum(), Eh[1].sum())); return dict(lo=El, hi=Eh)
acc, seen, _al = [], [], threading.Lock()
def site(s, u):
    F = fit(s, u)
    if F is None: return
    m = rows(u)
    with _al: anc = [A for A in acc if u.startswith(A['level']) and A['level'] != u and dice(A['mask'], s['mask']) >= DICE]
    for A in anc:
        g = max(abs(spear(A['lo'][0], F['lo'][0], m)), abs(spear(A['lo'][0], F['hi'][0], m)))
        log('site %s at %s vs %s pooled at %s: Spearman on %s particles %.3f' % (s['name'], u, A['name'], A['level'], u, g))
        if g >= GATE:
            with _al: seen.append((s['name'], u, A['name'], g))
            return
    L, FL = u, F
    while par(L) is not None:
        p = par(L); B = '%s/fits/%s@%s' % (O, s['name'], p); El = em(s['mask'], p, FL['lo'][1], B + '/lo')
        if El is None: break
        g = spear(El[0], F['lo'][0], m); log('site %s: pooling test %s -> %s, Spearman on %s particles %.3f (gate %.2f)' % (s['name'], L, p, u, g, GATE))
        if g < GATE: break
        Eh = em(s['mask'], p, FL['hi'][1], B + '/hi')
        if Eh is None: break
        L, FL = p, dict(lo=El, hi=Eh)
    with _al: acc.append(dict(name=s['name'], mask=s['mask'], home=u, level=L, lo=FL['lo'], hi=FL['hi'], rep=s['rep'], repro=s['repro']))
    log('site %s (home %s, %s, %.0f %% of the molecule) -> level %s' % (s['name'], u, s['rep'], 100 * s['mol'], L))
def site_safe(t):
    try: site(*t)
    except Exception as e: log('site %s at %s ERROR %r' % (t[0]['name'], t[1], e))
for dep in sorted(set(len(u) for u in todo)):
    lev = [u for u in todo if len(u) == dep]; cand = {u: discover(u) for u in lev}
    log('depth %d: %s' % (dep - 1, '; '.join('%s: %s' % (u, ' '.join('%s(%s, %.0f%%)' % (s['name'], s['rep'], 100 * s['mol']) for s in cand[u]) or '-') for u in lev)))
    with ThreadPoolExecutor(max(1, sum(len(c) for c in cand.values()))) as ex: list(ex.map(site_safe, [(s, u) for u in lev for s in cand[u]]))
keep = []
for A in sorted(acc, key=lambda A: -A['repro']):
    if any(B['level'] == A['level'] and dice(A['mask'], B['mask']) >= DICE and abs(spear(A['lo'][0], B['lo'][0], rows(A['level']))) >= GATE for B in keep): log('site %s duplicates a kept site at %s' % (A['name'], A['level'])); continue
    keep.append(A)
os.makedirs(O + '/labels', exist_ok=True)
with open(O + '/sites_final.txt', 'w') as f:
    f.write('# name home level representative repro n_lo n_hi labels mask\n')
    for A in keep:
        m = rows(A['level']); l = np.zeros(NM, int); l[A['lo'][1] & m] = 1; l[A['hi'][1] & m & ~A['lo'][1]] = 2; both = A['lo'][1] & A['hi'][1] & m; l[both] = np.where(A['lo'][0][both] >= A['hi'][0][both], 1, 2)
        lf = '%s/labels/%s.txt' % (O, A['name']); np.savetxt(lf, l, fmt='%d')
        f.write('%s %s %s %s %.3f %d %d %s %s\n' % (A['name'], A['home'], A['level'], A['rep'], A['repro'], (l == 1).sum(), (l == 2).sum(), lf, A['mask']))
for t in seen: log('site %s at %s = %s (Spearman %.3f)' % t)
log('%d sites kept: %s' % (len(keep), '; '.join('%s home %s level %s' % (A['name'], A['home'], A['level']) for A in keep)))
if keep and os.path.exists(R + '/final_labels.txt'):
    run('%s %s/hybrid_labels.py %s/final_labels.txt %s %s/hybrid_labels.txt --minn 500' % (PY, T, R, ','.join('%s/labels/%s.txt' % (O, A['name']) for A in keep), O), O + '/hybrid.log'); log(open(O + '/hybrid.log').read().strip())
log('DONE')
