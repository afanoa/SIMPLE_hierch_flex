#!/usr/bin/env python3
"""usage: mask_series.py --base SITE_REPOSED_DIR --seed SEED --tree DIR --node U --out DIR --h PREFIX --cache P --ports NAME=NODE_DIR .. [--mskdiam 232 --cons C --closed C]"""
import argparse, glob, os, subprocess, sys, threading, time, numpy as np
from concurrent.futures import ThreadPoolExecutor
from scipy import ndimage
ap = argparse.ArgumentParser()
for k, d in (('base', None), ('seed', None), ('tree', None), ('node', 'r1'), ('out', None), ('h', None), ('cache', None), ('mra', ''), ('mrb', ''), ('ctfabs', '1'), ('mskdiam', '232'),
             ('pars', '4'), ('iters', '10'), ('em_lim', '15'), ('cons', ''), ('closed', ''), ('mapn', '2500'), ('fall', '6')):
    ap.add_argument('--' + k, default=d, required=d is None)
ap.add_argument('--ports', nargs='*', default=[]); a = ap.parse_args(); T = os.path.dirname(os.path.abspath(__file__)); PY = sys.executable; O = a.out; os.makedirs(O + '/masks', exist_ok=True)
MRA = a.mra or a.h + '_A/meta_rows.txt'; MRB = a.mrb or a.h + '_B/meta_rows.txt'; LOG = open(O + '/mask_series.log', 'a'); _ll = threading.Lock()
def log(*s):
    m = '[%s mask series] %s' % (time.strftime('%F %T'), ' '.join(str(x) for x in s))
    with _ll: print(m, flush=True); LOG.write(m + '\n'); LOG.flush()
def run(cmd, out): return subprocess.run('%s > %s 2>&1' % (cmd, out), shell=True).returncode == 0
B, R, U = a.base, a.tree, a.node; use = np.loadtxt(B + '/rows_use.txt', dtype=int) > 0; NM = len(use); PA, PB = B + '/poses_A.txt', B + '/poses_B.txt'
NPJ = '%s/node_%s/n.simple' % (R, U); MEAN = sorted(glob.glob('%s/node_%s/mean/*_reconstruct3D' % (R, U)), key=os.path.getmtime)[-1] + '/recvol_state01.mrc'
seed = np.loadtxt(a.seed, dtype=int) > 0; BASEM = os.path.realpath(B + '/picked/site_02.mrc')
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n)
def wr(p, x, ref):
    with open(ref, 'rb') as f: h = f.read(1024)
    hh = np.frombuffer(h, np.int32).copy(); hh[23] = 0; hf = hh.view(np.float32); hf[19:22] = x.min(), x.max(), x.mean()
    with open(p, 'wb') as f: f.write(hh.tobytes()); f.write(x.astype(np.float32).tobytes())
def dice(p, q): x, y = rd(p) > 0.5, rd(q) > 0.5; return 2. * (x & y).sum() / max(x.sum() + y.sum(), 1)
M0 = rd(BASEM); n = M0.shape[0]; c = np.arange(n) - n // 2; X, Y, Z = np.meshgrid(c, c, c, indexing='ij'); sph = X ** 2 + Y ** 2 + Z ** 2 <= (float(a.mskdiam) / 2 / 1.464) ** 2
def falloff(core): d = ndimage.distance_transform_edt(~core); return np.where(core, 1., np.exp(-d ** 2 / (2 * float(a.fall) ** 2))) * sph
masks = [('1f1_base', BASEM)]
core0 = M0 > 0.999
for k in (1, 2):
    for tag, op in (('dil', ndimage.binary_dilation), ('ero', ndimage.binary_erosion)):
        p = '%s/masks/1f1_%s%d.mrc' % (O, tag, k)
        if not os.path.exists(p): wr(p, falloff(op(core0, iterations=k)), BASEM)
        masks.append(('1f1_%s%d' % (tag, k), p))
def disc(tag, view, frac):
    D = '%s/disc_%s' % (O, tag)
    if not os.path.exists(D + '/candidates.txt'): run('%s %s/discover_sites2.py --tree %s --out %s --mskdiam %s --frac %s --maxmol 1 --dice 0.5 --allcands 1' % (PY, T, view, D, a.mskdiam, frac), D + '.log')
    C = [w for w in (l.split() for l in open(D + '/candidates.txt')) if w and not w[0].startswith('#') and ':states:' in w[1]] if os.path.exists(D + '/candidates.txt') else []
    if not C: log('discovery %s found no state candidates' % tag); return None
    best = max(C, key=lambda w: dice(w[-1], BASEM)); log('discovery %s (frac %s): best match to 1f1 = %s, Dice %.3f' % (tag, frac, best[1], dice(best[-1], BASEM))); return best[-1]
V = '%s/view_script' % O; os.makedirs(V, exist_ok=True)
if not os.path.lexists(V + '/node_' + U): os.symlink('%s/node_%s' % (R, U), V + '/node_' + U)
for fr in ('0.4', '0.6'):
    p = disc('script_frac' + fr, V, fr)
    if p: masks.append(('1f1_frac%s' % fr, p))
masks.append(('script_nodemask', '%s/node_%s/rsmask.mrc' % (R, U)))
for e in a.ports:
    nm, nd = e.split('=', 1); masks.append(('%s_nodemask' % nm, nd + '/mask.mrc'))
    Vp = '%s/view_%s/node_%s' % (O, nm, U); os.makedirs(Vp + '/fstageB', exist_ok=True)
    if not os.path.exists(Vp + '/refined_f.txt'):
        lab = np.zeros(NM, int)
        for f in sorted(glob.glob(nd + '/leaf_r_s*_rows.txt')): lab[np.loadtxt(f, dtype=int) - 1] = int(f.split('_s')[-1].split('_')[0])
        np.savetxt(Vp + '/refined_f.txt', lab, fmt='%d'); open(Vp + '/node.done', 'w').close()
    if not os.path.lexists(Vp + '/fstageB/1_refine3D'): os.symlink(nd + '/recB', Vp + '/fstageB/1_refine3D')
    p = disc(nm, os.path.dirname(Vp), '0.2')
    if p: masks.append(('%s_1f1' % nm, p))
_es, _rs = threading.Semaphore(int(a.pars)), threading.Semaphore(2)
def em(t):
    nm, mk = t; D = '%s/em_%s' % (O, nm); os.makedirs(D, exist_ok=True)
    if not os.path.exists(D + '/curve.txt'):
        for h, mr in (('A', MRA), ('B', MRB)):
            rr = np.loadtxt(mr).astype(int) - 1; np.savetxt('%s/start_%s.txt' % (D, h), np.where(use[rr], 1 + (seed & use)[rr].astype(int), 0), fmt='%d')
        with _es:
            ok = run('%s -u %s/minor_em3.py --sim %s_A,%s_B --poses %s,%s --labels %s/start_A.txt,%s/start_B.txt --cache %s --ctfabs %s --cons_ref %s --mask %s --soft 1 --iters %s --em_lim %s --post 1 --nomaps 1 --out %s/em'
                     % (PY, T, a.h, a.h, PA, PB, D, D, a.cache, a.ctfabs, MEAN, mk, a.iters, a.em_lim, D), D + '/em.log') and \
                 run('%s -u %s/site_size.py --em %s/em --start %s/start_A.txt,%s/start_B.txt --sim %s_A,%s_B --poses %s,%s --mr %s,%s --nmaster %d --mask %s --out %s --cache %s --ctfabs %s --em_lim %s'
                     % (PY, T, D, D, D, a.h, a.h, PA, PB, MRA, MRB, NM, mk, D, a.cache, a.ctfabs, a.em_lim), D + '/size.log')
        if not ok: log('EM %s FAILED' % nm); return None
    p = np.loadtxt(D + '/post_master.txt'); ns = int(open(D + '/nstar.txt').read()); N = int(a.mapn); l = np.zeros(NM, int); q = np.where(use, p, -1.); l[np.argsort(-q)[:N]] = 1; l[q < 0] = 0
    R3 = '%s/rec_%s' % (O, nm); os.makedirs(R3, exist_ok=True)
    if not glob.glob(R3 + '/*_reconstruct3D/recvol_state01_even_unfil.mrc'):
        np.savetxt(R3 + '/labels.txt', l, fmt='%d')
        with _rs: run('cd %s && cp -f %s r.simple && simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=labels.txt mkdir=no && simple_exec prg=update_project projfile=r.simple qsys_name=local nparts=2 nthr=8 mkdir=no && simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=%s nstates=1 nparts=2 nthr=8' % (R3, NPJ, a.mskdiam), R3 + '/rec.log')
    g = glob.glob(R3 + '/*_reconstruct3D'); return dict(post=p, nstar=ns, rec=sorted(g)[-1] if g else None)
with ThreadPoolExecutor(len(masks)) as ex: RES = dict(zip([m[0] for m in masks], ex.map(em, masks)))
spec = ['%s %s/recvol_state01_even_unfil.mrc %s/recvol_state01_odd_unfil.mrc' % (nm, RES[nm]['rec'], RES[nm]['rec']) for nm, _ in masks if RES.get(nm) and RES[nm]['rec']]
open(O + '/spec.txt', 'w').write('\n'.join(spec) + '\n'); Vw = O + '/view'
run('%s %s/make_view.py %s %s/spec.txt --box 160' % (PY, T, Vw, O), O + '/view.log')
run('%s %s/arm_readout.py %s %s %s/arm.txt %s' % (PY, T, a.cons, a.closed, O, ' '.join('%s/%s.mrc' % (Vw, s.split()[0]) for s in spec)), O + '/arm.log')
arm = {os.path.basename(w[0])[:-4]: (float(w[1]), float(w[2])) for w in (l.split() for l in open(O + '/arm.txt')) if w and not w[0].startswith('#')} if os.path.exists(O + '/arm.txt') else {}
b = RES['1f1_base']['post']; m = use & (b >= 0); rb = np.argsort(np.argsort(-b[m]))
def jac(x, y, N): A, Bq = set(np.argsort(-x)[:N]), set(np.argsort(-y)[:N]); return len(A & Bq) / len(A | Bq)
rows = []
for nm, mk in masks:
    r = RES.get(nm)
    if not r: continue
    x = r['post']; mm = m & (x >= 0); sp = float(np.corrcoef(np.argsort(np.argsort(-b[mm])), np.argsort(np.argsort(-x[mm])))[0, 1])
    bb, xx = np.where(m, b, -1.), np.where(m, x, -1.); M = rd(mk)
    rows.append((dice(mk, BASEM), nm, float((M > 0.5).mean() / sph.mean()), sp, [jac(bb, xx, N) for N in (1000, 2500, 5000, 10000)], r['nstar'], arm.get(nm)))
print('closed-pole mask sensitivity: fixed start %s, shared input poses %s, same EM and sizing; only the site mask varies (baseline 1f1 = %s)' % (a.seed, NPJ, BASEM))
print('%-18s %6s %6s %8s %24s %7s %14s' % ('mask', 'Dice', 'vol', 'Spearman', 'Jaccard 1k/2.5k/5k/10k', 'N*', 'top-2.5k c/o'))
for d, nm, vol, sp, jc, ns, rr in sorted(rows, key=lambda t: -t[0]):
    print('%-18s %6.3f %5.1f%% %8.3f %24s %7d %14s' % (nm, d, 100 * vol, sp, ' / '.join('%.2f' % j for j in jc), ns, '%.2f / %.2f' % rr if rr else '-'))
