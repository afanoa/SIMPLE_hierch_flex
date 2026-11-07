#!/usr/bin/env python3
"""usage: half_poles.py --tree DIR --node U --out DIR --h HALFPREFIX --cache P [--nsn 3 --npole 7500 --mapn 1250 --cons C --closed C]"""
import argparse, glob, os, subprocess, sys, threading, time, numpy as np
from concurrent.futures import ThreadPoolExecutor
ap = argparse.ArgumentParser()
for k, d in (('tree', None), ('node', 'r1'), ('out', None), ('h', None), ('cache', None), ('mra', ''), ('mrb', ''), ('ctfabs', '1'), ('mskdiam', '232'), ('nsn', '3'), ('npole', '7500'),
             ('polefrac', '0.15'), ('k', '8'), ('lpf', '15'), ('np', '8'), ('pars', '4'), ('iters', '10'), ('em_lim', '15'), ('mapn', '1250'), ('cons', ''), ('closed', ''), ('lp', '15')):
    ap.add_argument('--' + k, default=d, required=d is None)
a = ap.parse_args(); T = os.path.dirname(os.path.abspath(__file__)); PY = sys.executable; R, U, O = a.tree, a.node, a.out; os.makedirs(O, exist_ok=True)
MRA = a.mra or a.h + '_A/meta_rows.txt'; MRB = a.mrb or a.h + '_B/meta_rows.txt'; LOG = open(O + '/half_poles.log', 'a'); _ll = threading.Lock()
def log(*s):
    m = '[%s half poles] %s' % (time.strftime('%F %T'), ' '.join(str(x) for x in s))
    with _ll: print(m, flush=True); LOG.write(m + '\n'); LOG.flush()
def run(cmd, out): return subprocess.run('%s > %s 2>&1' % (cmd, out), shell=True).returncode == 0
use = np.loadtxt('%s/rows_%s.txt' % (R, U), dtype=int) > 0; NM = len(use); NPJ = '%s/node_%s/n.simple' % (R, U); idx = np.arange(NM)
HALF = {'h1': use & ((idx // 2) % 2 == 0), 'h2': use & ((idx // 2) % 2 == 1)}; ST = np.loadtxt('%s/node_%s/refined_f.txt' % (R, U), dtype=int)
PA, PB = O + '/poses_A.txt', O + '/poses_B.txt'
if not os.path.exists(PB):
    run('simple_exec prg=print_project_field projfile=%s oritype=ptcl3D mkdir=no 2>/dev/null > %s/pm.txt; for h in A:%s B:%s; do awk "NR == FNR {l[FNR] = \\$0; next} {print l[\\$1 + 0]}" %s/pm.txt ${h#*:} > %s/poses_${h%%%%:*}.tmp && mv %s/poses_${h%%%%:*}.tmp %s/poses_${h%%%%:*}.txt; done; rm -f %s/pm.txt'
        % (NPJ, O, MRA, MRB, O, O, O, O, O), O + '/poses.log')
_rs, _fs, _es = threading.Semaphore(2), threading.Semaphore(2), threading.Semaphore(int(a.pars))
def rec(D, lab, ns):
    os.makedirs(D, exist_ok=True); g = glob.glob(D + '/*_reconstruct3D/recvol_state01_even_unfil.mrc')
    if not g:
        np.savetxt(D + '/labels.txt', lab, fmt='%d')
        with _rs: run('cd %s && cp -f %s r.simple && simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=labels.txt mkdir=no && simple_exec prg=update_project projfile=r.simple qsys_name=local nparts=2 nthr=8 mkdir=no && simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=%s nstates=%d nparts=2 nthr=8' % (D, NPJ, a.mskdiam, ns), D + '/rec.log')
        g = glob.glob(D + '/*_reconstruct3D/recvol_state01_even_unfil.mrc')
    return os.path.dirname(sorted(g, key=os.path.getmtime)[-1]) if g else None
def em(rows, mean, mask, seed, D):
    os.makedirs(D, exist_ok=True)
    if not os.path.exists(D + '/curve.txt'):
        for h, mr in (('A', MRA), ('B', MRB)):
            rr = np.loadtxt(mr).astype(int) - 1; np.savetxt('%s/start_%s.txt' % (D, h), np.where(rows[rr], 1 + (seed & rows)[rr].astype(int), 0), fmt='%d')
        with _es:
            ok = run('%s -u %s/minor_em3.py --sim %s_A,%s_B --poses %s,%s --labels %s/start_A.txt,%s/start_B.txt --cache %s --ctfabs %s --cons_ref %s --mask %s --soft 1 --iters %s --em_lim %s --post 1 --nomaps 1 --out %s/em'
                     % (PY, T, a.h, a.h, PA, PB, D, D, a.cache, a.ctfabs, mean, mask, a.iters, a.em_lim, D), D + '/em.log') and \
                 run('%s -u %s/site_size.py --em %s/em --start %s/start_A.txt,%s/start_B.txt --sim %s_A,%s_B --poses %s,%s --mr %s,%s --nmaster %d --mask %s --out %s --cache %s --ctfabs %s --em_lim %s'
                     % (PY, T, D, D, D, a.h, a.h, PA, PB, MRA, MRB, NM, mask, D, a.cache, a.ctfabs, a.em_lim), D + '/size.log')
        if not ok: log('EM %s FAILED' % D); return None
    return dict(post=np.loadtxt(D + '/post_master.txt'), nstar=int(open(D + '/nstar.txt').read()))
def half(h):
    rows = HALF[h]; D = '%s/%s' % (O, h); os.makedirs(D, exist_ok=True); np.savetxt(D + '/rows.txt', rows.astype(int), fmt='%d')
    MD = rec(D + '/mean', rows.astype(int), 1); SD = rec(D + '/states', np.where(rows, ST, 0), int(ST.max()))
    if MD is None or SD is None: log('%s: maps FAILED' % h); return None
    mean = MD + '/recvol_state01.mrc'; V = D + '/view/node_' + U; os.makedirs(V + '/fstageB', exist_ok=True); np.savetxt(V + '/refined_f.txt', np.where(rows, ST, 0), fmt='%d'); open(V + '/node.done', 'w').close()
    if not os.path.lexists(V + '/fstageB/1_refine3D'): os.symlink(SD, V + '/fstageB/1_refine3D')
    if not os.path.exists(D + '/disc/sites.txt'): run('%s %s/discover_sites2.py --tree %s/view --out %s/disc --mskdiam %s --maxmol 0.3 --dice 0.5' % (PY, T, D, D, a.mskdiam), D + '/disc.log')
    S = [w for w in (l.split() for l in open(D + '/disc/sites.txt')) if w and not w[0].startswith('#')][:int(a.nsn)] if os.path.exists(D + '/disc/sites.txt') else []
    log('%s: %d particles, sites %s' % (h, rows.sum(), ' '.join('%s(%s, %.0f%%)' % (w[0], w[1], 100 * float(w[4])) for w in S)))
    def site(w):
        j = int(w[0]); mk = '%s/disc/site_%02d.mrc' % (D, j); F = '%s/s%d/flex' % (D, j); os.makedirs(F, exist_ok=True); c = glob.glob(F + '/*_flex_pca/flex_pca_coordinates.txt')
        if not c:
            np.savetxt(F + '/rows.txt', rows.astype(int), fmt='%d')
            with _fs: run('cd %s && cp -f %s f.simple && simple_exec prg=selection projfile=f.simple oritype=ptcl3D infile=rows.txt mkdir=no && simple_exec prg=update_project projfile=f.simple qsys_name=slurm nparts=%s nthr=8 mkdir=no && simple_exec prg=flex_pca projfile=f.simple vol1=%s mskdiam=%s lp=%s npreimages=%s rec_states=no min_state_frac=0 umap=no n_probe_iters=2 pcg_mskfile=%s nparts=%s nthr=8'
                          % (F, NPJ, a.np, mean, a.mskdiam, a.lpf, a.k, mk, a.np), F + '/flex.log')
            c = glob.glob(F + '/*_flex_pca/flex_pca_coordinates.txt')
            if not c: log('%s site %d flex FAILED' % (h, j)); return []
        C = np.loadtxt(sorted(c, key=os.path.getmtime)[-1], comments='#'); r, z = C[:, 0].astype(int) - 1, C[:, 5]; k = int(min(int(a.npole), float(a.polefrac) * len(r))); o = np.argsort(z); out = []
        for side, ix in (('lo', o[:k]), ('hi', o[::-1][:k])):
            s = np.zeros(NM, bool); s[r[ix]] = True; e = em(rows, mean, mk, s, '%s/s%d/%s' % (D, j, side))
            if e is None: continue
            l = np.zeros(NM, int); p = np.where(rows, e['post'], -1.); l[np.argsort(-p)[:int(a.mapn)]] = 1; l[p < 0] = 0; RD = rec('%s/s%d/%s_map' % (D, j, side), l, 1)
            out.append(dict(half=h, site=j, rep=w[1], mask=mk, side=side, nstar=e['nstar'], rec=RD, mean=MD))
        log('%s site %d (%s): poles %s' % (h, j, w[1], ' '.join('%s N* %d' % (q['side'], q['nstar']) for q in out))); return out
    with ThreadPoolExecutor(len(S) or 1) as ex: return [q for L in ex.map(site, S) for q in L]
with ThreadPoolExecutor(2) as ex: P = dict(zip(('h1', 'h2'), ex.map(half, ('h1', 'h2'))))
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
spec = [('%s_s%d_%s' % (q['half'], q['site'], q['side']), q['rec']) for h in ('h1', 'h2') for q in (P[h] or []) if q['rec']]
open(O + '/spec.txt', 'w').write(''.join('%s %s/recvol_state01_even_unfil.mrc %s/recvol_state01_odd_unfil.mrc\n' % (n_, d_, d_) for n_, d_ in spec))
run('%s %s/make_view.py %s/view %s/spec.txt --box 160' % (PY, T, O, O), O + '/view.log')
run('%s %s/arm_readout.py %s %s %s/arm.txt %s' % (PY, T, a.cons, a.closed, O, ' '.join('%s/view/%s.mrc' % (O, n_) for n_, _ in spec)), O + '/arm.log')
arm = {os.path.basename(w[0])[:-4]: (float(w[1]), float(w[2])) for w in (l.split() for l in open(O + '/arm.txt')) if w and not w[0].startswith('#')} if os.path.exists(O + '/arm.txt') else {}
n = 160; f = np.fft.fftfreq(n, d=1.464); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); G = np.exp(-0.5 * (kx ** 2 + ky ** 2 + kz ** 2) * (float(a.lp) / 1.2) ** 2)
lp = lambda v: np.real(np.fft.ifftn(np.fft.fftn(v) * G)); c = np.arange(n) - n // 2; X, Y, W = np.meshgrid(c, c, c, indexing='ij')
def dev(q):
    me, mo = lp(rd(q['mean'] + '/recvol_state01_even_unfil.mrc')), lp(rd(q['mean'] + '/recvol_state01_odd_unfil.mrc')); out = []
    for x, m in ((lp(rd(q['rec'] + '/recvol_state01_even_unfil.mrc')), me), (lp(rd(q['rec'] + '/recvol_state01_odd_unfil.mrc')), mo)): s = (x[SUP] * m[SUP]).sum() / (x[SUP] ** 2).sum(); out.append(s * x - m)
    return out
m0 = P['h1'][0]['mean'] if P['h1'] else None
if m0:
    mf = lp(rd(m0 + '/recvol_state01.mrc')); SUP = (mf > 0.2 * mf.max()) & (X ** 2 + Y ** 2 + W ** 2 <= (float(a.mskdiam) / 2 / 1.464) ** 2)
    cr = lambda x, y: float(np.corrcoef(x[SUP], y[SUP])[0, 1]); D_ = {}
    for h in ('h1', 'h2'):
        for q in P[h] or []:
            if q['rec']: e, o = dev(q); r = cr(e, o); D_[(h, q['site'], q['side'])] = (0.5 * (e + o), max(2 * r / (1 + r), 1e-6))
    print('\nhalf-independent poles (h1 / h2 = disjoint halves of %s, each with its own state maps, discovery, site flex and EMs; maps top-%s at the shared input poses):' % (U, a.mapn))
    for q in P['h1'] or []:
        best = max(((cr(D_[('h1', q['site'], q['side'])][0], D_[('h2', x['site'], x['side'])][0]) / np.sqrt(D_[('h1', q['site'], q['side'])][1] * D_[('h2', x['site'], x['side'])][1]), x) for x in P['h2'] or [] if ('h2', x['site'], x['side']) in D_ and ('h1', q['site'], q['side']) in D_), key=lambda t: t[0], default=(np.nan, None))
        x = best[1]; nm1 = 'h1_s%d_%s' % (q['site'], q['side'])
        if x is None: print('  %-12s no partner' % nm1); continue
        nm2 = 'h2_s%d_%s' % (x['site'], x['side']); dm = (lambda p_, q_: 2. * ((rd(p_) > 0.5) & (rd(q_) > 0.5)).sum() / max((rd(p_) > 0.5).sum() + (rd(q_) > 0.5).sum(), 1))(q['mask'], x['mask'])
        print('  %-12s (%s, N* %d, %s) <-> %-12s (%s, N* %d, %s): mask Dice %.2f, map match %.2f of the half-map noise ceiling'
              % (nm1, q['rep'], q['nstar'], '%.2f/%.2f' % arm[nm1] if nm1 in arm else '-', nm2, x['rep'], x['nstar'], '%.2f/%.2f' % arm[nm2] if nm2 in arm else '-', dm, best[0]))
log('DONE')
