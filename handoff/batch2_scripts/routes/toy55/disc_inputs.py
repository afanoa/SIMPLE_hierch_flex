#!/usr/bin/env python3
"""usage: disc_inputs.py --base SITE_REPOSED_DIR --basepost POST_MASTER --seed SEED --tree DIR --node U --out DIR --h PREFIX --cache P --ports NAME=PORT_NODE_DIR .."""
import argparse, glob, os, subprocess, sys, threading, time, numpy as np
from concurrent.futures import ThreadPoolExecutor
ap = argparse.ArgumentParser()
for k, d in (('base', None), ('basepost', None), ('seed', None), ('tree', None), ('node', 'r1'), ('out', None), ('h', None), ('cache', None), ('mra', ''), ('mrb', ''), ('ctfabs', '1'),
             ('mskdiam', '232'), ('pars', '4'), ('iters', '10'), ('em_lim', '15')):
    ap.add_argument('--' + k, default=d, required=d is None)
ap.add_argument('--ports', nargs='*', default=[]); a = ap.parse_args(); T = os.path.dirname(os.path.abspath(__file__)); PY = sys.executable; O = a.out; os.makedirs(O, exist_ok=True)
MRA = a.mra or a.h + '_A/meta_rows.txt'; MRB = a.mrb or a.h + '_B/meta_rows.txt'; LOG = open(O + '/disc_inputs.log', 'a'); _ll = threading.Lock()
def log(*s):
    m = '[%s disc inputs] %s' % (time.strftime('%F %T'), ' '.join(str(x) for x in s))
    with _ll: print(m, flush=True); LOG.write(m + '\n'); LOG.flush()
def run(cmd, out): return subprocess.run('%s > %s 2>&1' % (cmd, out), shell=True).returncode == 0
B, R, U = a.base, a.tree, a.node; use = np.loadtxt(B + '/rows_use.txt', dtype=int) > 0; NM = len(use); PA, PB = B + '/poses_A.txt', B + '/poses_B.txt'
MEAN = sorted(glob.glob('%s/node_%s/mean/*_reconstruct3D' % (R, U)), key=os.path.getmtime)[-1] + '/recvol_state01.mrc'; BASEM = os.path.realpath(B + '/picked/site_02.mrc')
seed = np.loadtxt(a.seed, dtype=int) > 0; P0 = np.loadtxt(a.basepost)
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n)
def dice(p, q): x, y = rd(p) > 0.5, rd(q) > 0.5; return 2. * (x & y).sum() / max(x.sum() + y.sum(), 1)
def lab_coords(f):
    C = np.loadtxt(f, comments='#'); l = np.zeros(NM, int); l[C[:, 0].astype(int) - 1] = C[:, 2].astype(int); return l
def lab_project(p, D):
    out = D + '/proj_states.txt'
    if not os.path.exists(out) or os.path.getsize(out) == 0: run('(simple_exec prg=print_project_field projfile=%s oritype=ptcl3D mkdir=no 2>/dev/null | awk \'{s = 0; for (i = 1; i <= NF; i++) if ($i ~ /^state=/) {split($i, q, "="); s = q[2] + 0}; print s}\' > %s)' % (p, out), D + '/proj.log')
    return np.loadtxt(out, dtype=int)
def lab_leaves(nd):
    l = np.zeros(NM, int)
    for f in sorted(glob.glob(nd + '/leaf_r_s*_rows.txt')): l[np.loadtxt(f, dtype=int) - 1] = int(f.split('_s')[-1].split('_')[0])
    return l
SN = '%s/node_%s' % (R, U); INP = []
INP += [('script', 'fflex', SN + '/focus/1_flex_pca', lambda D: np.loadtxt(SN + '/fflex_states_raw.txt', dtype=int)),
        ('script', 'cheap', SN + '/cheap/rec/1_reconstruct3D', lambda D: np.loadtxt(SN + '/cheap/children_1.txt', dtype=int)),
        ('script', 'recB', sorted(glob.glob(SN + '/fstageB/*_refine3D'), key=os.path.getmtime)[-1], lambda D: np.loadtxt(SN + '/refined_f.txt', dtype=int))]
for e in a.ports:
    nm, nd = e.split('=', 1)
    INP += [(nm, 'fflex', nd + '/fflex', lambda D, nd=nd: lab_coords(nd + '/fflex/flex_pca_coordinates.txt')), (nm, 'cheap', nd + '/cheap', lambda D, nd=nd: lab_project(nd + '/cheap/c.simple', D)),
            (nm, 'recB', nd + '/recB', lambda D, nd=nd: lab_leaves(nd))]
def discover(t):
    run_, kind, maps, labf = t; D = '%s/%s_%s' % (O, run_, kind); V = D + '/view/node_' + U; os.makedirs(V + '/fstageB', exist_ok=True)
    if not os.path.exists(V + '/refined_f.txt'): l = labf(D); np.savetxt(V + '/refined_f.txt', np.where(use, l, 0), fmt='%d'); open(V + '/node.done', 'w').close()
    if not os.path.lexists(V + '/fstageB/1_refine3D'): os.symlink(maps, V + '/fstageB/1_refine3D')
    if not os.path.exists(D + '/disc/sites.txt'): run('%s %s/discover_sites2.py --tree %s/view --out %s/disc --mskdiam %s --maxmol 0.3 --dice 0.5 --allcands 1' % (PY, T, D, D, a.mskdiam), D + '/disc.log')
    S = [w for w in (l.split() for l in open(D + '/disc/sites.txt')) if w and not w[0].startswith('#')] if os.path.exists(D + '/disc/sites.txt') else []
    C = [w for w in (l.split() for l in open(D + '/disc/candidates.txt')) if w and not w[0].startswith('#')] if os.path.exists(D + '/disc/candidates.txt') else []
    sites = [dict(rank=int(w[0]), rep=w[1], mol=float(w[4]), mask='%s/disc/site_%02d.mrc' % (D, int(w[0]))) for w in S]
    for s in sites: s['dice'] = dice(s['mask'], BASEM)
    best = max(C, key=lambda w: dice(w[-1], BASEM)) if C else None
    log('%s %s: %d sites (%s); best candidate overall %s Dice %.3f' % (run_, kind, len(sites), ' '.join('%d:%s %.0f%% Dice %.2f' % (s['rank'], s['rep'], 100 * s['mol'], s['dice']) for s in sites[:4]),
                                                                      best[1] if best else '-', dice(best[-1], BASEM) if best else 0.))
    return dict(run=run_, kind=kind, D=D, sites=sites, best=(best[1], best[-1], dice(best[-1], BASEM)) if best else None)
with ThreadPoolExecutor(4) as ex: DISC = list(ex.map(discover, INP))
_es = threading.Semaphore(int(a.pars))
def em(t):
    tag, mk = t; D = '%s/em_%s' % (O, tag); os.makedirs(D, exist_ok=True)
    if not os.path.exists(D + '/curve.txt'):
        for h, mr in (('A', MRA), ('B', MRB)):
            rr = np.loadtxt(mr).astype(int) - 1; np.savetxt('%s/start_%s.txt' % (D, h), np.where(use[rr], 1 + (seed & use)[rr].astype(int), 0), fmt='%d')
        with _es:
            ok = run('%s -u %s/minor_em3.py --sim %s_A,%s_B --poses %s,%s --labels %s/start_A.txt,%s/start_B.txt --cache %s --ctfabs %s --cons_ref %s --mask %s --soft 1 --iters %s --em_lim %s --post 1 --nomaps 1 --out %s/em'
                     % (PY, T, a.h, a.h, PA, PB, D, D, a.cache, a.ctfabs, MEAN, mk, a.iters, a.em_lim, D), D + '/em.log') and \
                 run('%s -u %s/site_size.py --em %s/em --start %s/start_A.txt,%s/start_B.txt --sim %s_A,%s_B --poses %s,%s --mr %s,%s --nmaster %d --mask %s --out %s --cache %s --ctfabs %s --em_lim %s'
                     % (PY, T, D, D, D, a.h, a.h, PA, PB, MRA, MRB, NM, mk, D, a.cache, a.ctfabs, a.em_lim), D + '/size.log')
        if not ok: return None
    p = np.loadtxt(D + '/post_master.txt'); m = use & (p >= 0) & (P0 >= 0)
    return float(np.corrcoef(np.argsort(np.argsort(-p[m])), np.argsort(np.argsort(-P0[m])))[0, 1]), int(open(D + '/nstar.txt').read())
jobs = []
for d in DISC:
    if d['sites']: jobs.append(('%s_%s_rule' % (d['run'], d['kind']), d['sites'][0]['mask']))
    if d['best'] and (not d['sites'] or os.path.realpath(d['best'][1]) != os.path.realpath(d['sites'][0]['mask']) and d['best'][2] > d['sites'][0]['dice'] + 0.05): jobs.append(('%s_%s_best' % (d['run'], d['kind']), d['best'][1]))
with ThreadPoolExecutor(len(jobs) or 1) as ex: EMR = dict(zip([j[0] for j in jobs], ex.map(em, jobs)))
print('discovery on pre-stage-A inputs vs recB: rule = site 1 of discover_sites2 (local <= 30 %%, ranked by reproducibility, ties by variance share); EM from the fixed closed start at the shared input poses, Spearman vs the 1f1 baseline pole')
print('%-7s %-6s %-22s %9s %13s %7s | %-22s %9s %13s' % ('run', 'input', 'rule pick (site 1)', 'Dice 1f1', 'pole Spearman', 'N*', 'best candidate', 'Dice 1f1', 'pole Spearman'))
for d in DISC:
    s = d['sites'][0] if d['sites'] else None; r = EMR.get('%s_%s_rule' % (d['run'], d['kind'])); b = EMR.get('%s_%s_best' % (d['run'], d['kind']))
    print('%-7s %-6s %-22s %9s %13s %7s | %-22s %9s %13s' % (d['run'], d['kind'], s['rep'] if s else '-', '%.3f' % s['dice'] if s else '-', '%.3f' % r[0] if r else '-', r[1] if r else '-',
          d['best'][0] if d['best'] else '-', '%.3f' % d['best'][2] if d['best'] else '-', '%.3f' % b[0] if b else ('= rule' if d['best'] else '-')))
