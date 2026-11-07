#!/usr/bin/env python3
"""usage: pole_dirs.py --tree DIR --node U --coords flex_pca_coordinates.txt --out DIR --h HALFPREFIX --cache P [--dmax 5 --npole 15000 --hold 0.2 --same 0.99 --cons C --closed C]"""
import argparse, glob, os, subprocess, sys, threading, time, numpy as np
from concurrent.futures import ThreadPoolExecutor
ap = argparse.ArgumentParser()
for k, d in (('tree', None), ('node', 'r'), ('coords', None), ('out', None), ('h', None), ('cache', None), ('mra', ''), ('mrb', ''), ('ctfabs', '1'), ('mskdiam', '232'), ('dmax', '5'),
             ('npole', '15000'), ('hold', '0.2'), ('same', '0.99'), ('pars', '4'), ('iters', '10'), ('em_lim', '15'), ('cons', ''), ('closed', ''), ('mapn', '2500')):
    ap.add_argument('--' + k, default=d, required=d is None)
a = ap.parse_args(); T = os.path.dirname(os.path.abspath(__file__)); PY = sys.executable; R, U, O = a.tree, a.node, a.out; os.makedirs(O, exist_ok=True)
MRA = a.mra or a.h + '_A/meta_rows.txt'; MRB = a.mrb or a.h + '_B/meta_rows.txt'; SAME, HOLD = float(a.same), float(a.hold); LOG = open(O + '/pole_dirs.log', 'a'); _ll = threading.Lock()
def log(*s):
    m = '[%s pole dirs] %s' % (time.strftime('%F %T'), ' '.join(str(x) for x in s))
    with _ll: print(m, flush=True); LOG.write(m + '\n'); LOG.flush()
def run(cmd, out): return subprocess.run('%s > %s 2>&1' % (cmd, out), shell=True).returncode == 0
use = np.loadtxt('%s/rows_%s.txt' % (R, U), dtype=int) > 0; NM = len(use); NPJ = '%s/node_%s/n.simple' % (R, U)
MEAN = sorted(glob.glob('%s/node_%s/mean/*_reconstruct3D' % (R, U)), key=os.path.getmtime)[-1] + '/recvol_state01.mrc'
PA, PB = O + '/poses_A.txt', O + '/poses_B.txt'
if not os.path.exists(PB):
    run('simple_exec prg=print_project_field projfile=%s oritype=ptcl3D mkdir=no 2>/dev/null > %s/pm.txt; for h in A:%s B:%s; do awk "NR == FNR {l[FNR] = \\$0; next} {print l[\\$1 + 0]}" %s/pm.txt ${h#*:} > %s/poses_${h%%%%:*}.tmp && mv %s/poses_${h%%%%:*}.tmp %s/poses_${h%%%%:*}.txt; done; rm -f %s/pm.txt'
        % (NPJ, O, MRA, MRB, O, O, O, O, O), O + '/poses.log')
C = np.loadtxt(a.coords, comments='#'); r0 = C[:, 0].astype(int) - 1; k = use[r0]; r0, Z = r0[k], C[k, 5:]
log('node %s: %d rows, latent %s (%d coordinates), shared input poses %s' % (U, len(r0), a.coords, Z.shape[1], NPJ))
def rec(D, lab, ns):
    os.makedirs(D, exist_ok=True); g = glob.glob(D + '/*_reconstruct3D/recvol_state01_even_unfil.mrc')
    if not g:
        np.savetxt(D + '/labels.txt', lab, fmt='%d')
        with _rs: run('cd %s && cp -f %s r.simple && simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=labels.txt mkdir=no && simple_exec prg=update_project projfile=r.simple qsys_name=local nparts=2 nthr=8 mkdir=no && simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=%s nstates=%d nparts=2 nthr=8'
                      % (D, NPJ, a.mskdiam, ns), D + '/rec.log')
        g = glob.glob(D + '/*_reconstruct3D/recvol_state01_even_unfil.mrc')
    return os.path.dirname(sorted(g, key=os.path.getmtime)[-1]) if g else None
_es, _rs = threading.Semaphore(int(a.pars)), threading.Semaphore(2)
def em(mask, seed, D):
    os.makedirs(D, exist_ok=True)
    if not os.path.exists(D + '/curve.txt'):
        for h, mr in (('A', MRA), ('B', MRB)):
            rr = np.loadtxt(mr).astype(int) - 1; np.savetxt('%s/start_%s.txt' % (D, h), np.where(use[rr], 1 + (seed & use)[rr].astype(int), 0), fmt='%d')
        with _es:
            ok = run('%s -u %s/minor_em3.py --sim %s_A,%s_B --poses %s,%s --labels %s/start_A.txt,%s/start_B.txt --cache %s --ctfabs %s --cons_ref %s --mask %s --soft 1 --iters %s --em_lim %s --post 1 --nomaps 1 --out %s/em'
                     % (PY, T, a.h, a.h, PA, PB, D, D, a.cache, a.ctfabs, MEAN, mask, a.iters, a.em_lim, D), D + '/em.log') and \
                 run('%s -u %s/site_size.py --em %s/em --start %s/start_A.txt,%s/start_B.txt --sim %s_A,%s_B --poses %s,%s --mr %s,%s --nmaster %d --mask %s --out %s --cache %s --ctfabs %s --em_lim %s'
                     % (PY, T, D, D, D, a.h, a.h, PA, PB, MRA, MRB, NM, mask, D, a.cache, a.ctfabs, a.em_lim), D + '/size.log')
        if not ok: log('EM %s FAILED' % D); return None
    cv = np.loadtxt(D + '/curve.txt', comments='#', ndmin=2); ns = int(open(D + '/nstar.txt').read()); hold = float(cv[np.argmin(np.abs(cv[:, 0] - ns)), 2])
    return dict(post=np.loadtxt(D + '/post_master.txt'), top=np.loadtxt(D + '/top_star.txt', dtype=int) > 0, nstar=ns, hold=hold, dir=D)
def spear(x, y):
    m = use & (x >= 0) & (y >= 0); return float(np.corrcoef(np.argsort(np.argsort(x[m])), np.argsort(np.argsort(y[m])))[0, 1])
def direction(d):
    Dd = '%s/d%d' % (O, d); os.makedirs(Dd, exist_ok=True); o = np.argsort(Z[:, d - 1]); n = int(a.npole); S = {}
    for nm, ix in (('lo', o[:n]), ('hi', o[::-1][:n]), ('lo2', o[:n // 2]), ('hi2', o[::-1][:n // 2])): S[nm] = np.zeros(NM, bool); S[nm][r0[ix]] = True
    lab = np.zeros(NM, int); lab[S['lo']] = 1; lab[S['hi']] = 2; RD = rec(Dd + '/tails', lab, 2)
    if RD is None: log('direction %d: tail maps FAILED' % d); return None
    with open(Dd + '/tails_list.txt', 'w') as f:
        for j, nm in ((1, 'lo'), (2, 'hi')): f.write('%s %s/recvol_state%02d_even_unfil.mrc %s/recvol_state%02d_odd_unfil.mrc %d\n' % (nm, RD, j, RD, j, S[nm].sum()))
    MK = Dd + '/mask.mrc'
    if not os.path.exists(MK): run('%s %s/arm_mask2.py %s --list %s/tails_list.txt --mskdiam %s --frac 0.2 --dil 3 --falloff 6 %s %s' % (PY, T, MK, Dd, a.mskdiam, '--cons ' + a.cons if a.cons else '', '--closed ' + a.closed if a.closed else ''), Dd + '/mask.log')
    if not os.path.exists(MK): log('direction %d: mask FAILED' % d); return None
    ml = [l for l in open(Dd + '/mask.log').read().splitlines() if l.strip() and 'Warning' not in l]; log('direction %d: tail-difference mask: %s' % (d, ml[-1][:170] if ml else '-'))
    with ThreadPoolExecutor(4) as ex: E = dict(zip(('lo', 'hi', 'lo2', 'hi2'), ex.map(lambda nm: em(MK, S[nm], '%s/%s' % (Dd, nm)), ('lo', 'hi', 'lo2', 'hi2'))))
    return dict(d=d, MK=MK, E=E, n=n)
with ThreadPoolExecutor(int(a.dmax)) as ex: DIRS = list(ex.map(direction, range(1, min(int(a.dmax), Z.shape[1]) + 1)))
acc, kept = [], []
for X in DIRS:
    if X is None: continue
    d, MK, E, n = X['d'], X['MK'], X['E'], X['n']; new_here = 0
    for side in ('lo', 'hi'):
        e, e2 = E[side], E[side + '2']
        if e is None or e2 is None: continue
        sa = spear(e['post'], e2['post']); sim = [(A['name'], spear(e['post'], A['post'])) for A in acc]; best = max(sim, key=lambda t: t[1]) if sim else ('-', 0.)
        ok = sa >= SAME and e['hold'] >= HOLD; new = best[1] < SAME
        log('direction %d %s: N* %d, hold-out cross-half corr %.3f, start agreement (15k vs %dk tail) %.3f, closest earlier pole %s %.3f -> %s'
            % (d, side, e['nstar'], e['hold'], n // 2000, sa, best[0], best[1], 'NEW POLE' if ok and new else ('same as ' + best[0] if not new else 'not reproducible')))
        if ok and new: acc.append(dict(name='d%d%s' % (d, side), d=d, side=side, mask=MK, **e)); new_here += 1
    if new_here: kept.append(d)
log('%d directions kept (%s), %d poles: %s' % (len(kept), ' '.join(map(str, kept)), len(acc), ' '.join('%s(N* %d, hold %.2f)' % (A['name'], A['nstar'], A['hold']) for A in acc)))
def pmap(t):
    A, N, tag = t; l = np.zeros(NM, int); p = np.where(use, A['post'], -1.); l[np.argsort(-p)[:N]] = 1; l[p < 0] = 0; RD = rec('%s/maps/%s_%s' % (O, A['name'], tag), l, 1)
    return '%s_%s %s/recvol_state01_even_unfil.mrc %s/recvol_state01_odd_unfil.mrc' % (A['name'], tag, RD, RD) if RD else None
with ThreadPoolExecutor(2) as ex: spec = [x for x in ex.map(pmap, [(A, N, tag) for A in acc for N, tag in ((A['nstar'], 'nstar'), (int(a.mapn), 'top%d' % int(a.mapn)))]) if x]
with open(O + '/poles.txt', 'w') as f:
    f.write('# name direction side nstar holdout post_master top_star mask\n')
    for A in acc: f.write('%s %d %s %d %.4f %s/post_master.txt %s/top_star.txt %s\n' % (A['name'], A['d'], A['side'], A['nstar'], A['hold'], A['dir'], A['dir'], A['mask']))
open(O + '/spec_poles.txt', 'w').write('\n'.join(spec) + '\n'); log('DONE')
