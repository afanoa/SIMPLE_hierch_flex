#!/usr/bin/env python3
"""Sample a run's state landscape with a particle budget: every state of a consolidation (e.g. the zoom leaves) is placed by
its map's deviation from its family's parent inside the family's periphery (scale-matched, low-passed, particle-weighted PCA);
then the smallest group is folded into its nearest group until every group holds >= NMIN particles. Extremes only merge with
near neighbours, so the ends of every motion survive and the groups between them are the intermediates; the number of groups
follows the number of particles. NMIN comes from the states' own resolution inside the periphery vs particle count (half-map FSC 0.143, fit
res = a - b ln N) at TARGET A, unless --nmin is given. No truth, no dataset-specific mask.
Writes OUT/states.txt (rows of the consolidation's states.txt -> group), OUT/groups.txt, OUT/landscape.txt.
usage: landscape_sample.py CONS PARENT_PREFIX PERI_PREFIX OUT [--target 8 --nmin N --lp 12 --var 0.9 --smpd 1.464]"""
import sys, os, re, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc
ap = argparse.ArgumentParser(); ap.add_argument('cons'); ap.add_argument('parent'); ap.add_argument('peri'); ap.add_argument('out')
ap.add_argument('--target', type=float, default=8.); ap.add_argument('--nmin', type=int, default=0); ap.add_argument('--lp', type=float, default=12.)
ap.add_argument('--var', type=float, default=0.9); ap.add_argument('--smpd', type=float, default=1.464)
ap.add_argument('--npc', type=int, default=3, help='sampling space = the leading NPC components (the dominant motions)')
ap.add_argument('--mode', default='fold', help='fold: the smallest group joins its nearest until all >= NMIN | fps: farthest-point seeds (extremes first), nearest-seed groups, then fold')
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
L1 = a.cons + '/fam0/regionall/level1/maps'
st = np.loadtxt(a.cons + '/states.txt', dtype=int); run = os.path.abspath(a.cons)
while run != '/' and not os.path.isdir(run + '/expand'): run = os.path.dirname(run)
famrow = np.zeros(len(st), int)
for F in (1, 2):
    fs = '%s/expand/family%d_subset.txt' % (run, F)
    if os.path.isfile(fs): famrow[np.loadtxt(fs, dtype=int).ravel() - 1] = F
G = {}
for l in open(a.cons + '/groups.txt'):   # a state's family = the family holding most of its particles
    if l.startswith('#') or not l.strip(): continue
    w = l.split(); s = int(w[0]); fam = np.bincount(famrow[st == s], minlength=3)[1:]
    if fam.sum(): G[s] = (int(w[1]), w[2], int(np.argmax(fam)) + 1)
def half(s, h):
    for d in (L1, L1 + '/halves'):
        f = '%s/recvol_state%02d%s.mrc' % (d, s, h)
        if os.path.isfile(f): return f
    return None
n = read_mrc(a.parent.replace('{F}', '1') + '.mrc').shape[0]; x = np.arange(n) - n // 2
rr = np.sqrt(x[:, None, None] ** 2 + x[None, :, None] ** 2 + x[None, None, :] ** 2); sph = np.clip((0.45 * n - rr) / (0.05 * n), 0, 1)
f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); k2 = kx ** 2 + ky ** 2 + kz ** 2
Glp = np.exp(-0.5 * k2 * (a.lp / 1.2) ** 2 / a.smpd ** 2); sh = np.minimum(np.round(np.sqrt(k2) * n).astype(int), n // 2)
def lp(v): return np.real(np.fft.ifftn(np.fft.fftn(v) * Glp))
def fsc_res(e, o):   # FSC 0.143 of two half maps under a soft sphere
    E, O = np.fft.fftn(e * sph), np.fft.fftn(o * sph); num = np.bincount(sh.ravel(), np.real(E * np.conj(O)).ravel(), n // 2 + 1)
    den = np.sqrt(np.bincount(sh.ravel(), (np.abs(E) ** 2).ravel(), n // 2 + 1) * np.bincount(sh.ravel(), (np.abs(O) ** 2).ravel(), n // 2 + 1)); c = num / np.maximum(den, 1e-30)
    below = np.where(c[1:n // 2] < 0.143)[0]; k = (below[0] + 1) if len(below) else n // 2 - 1; return n * a.smpd / max(k, 1)
out_lab = np.zeros_like(st); gid = 0; lines = []; land = []
MK = {F: np.clip(read_mrc(a.peri.replace('{F}', str(F)) + '.mrc'), 0, 1).astype(np.float32) for F in {v[2] for v in G.values()}}
res_all = {s: fsc_res(read_mrc(half(s, '_even')) * MK[G[s][2]], read_mrc(half(s, '_odd')) * MK[G[s][2]]) for s in G if half(s, '_even') and half(s, '_odd')}   # resolution where the states differ
ok = [s for s in res_all if res_all[s] < 40]; N_ = np.array([G[s][0] for s in ok], float); R_ = np.array([res_all[s] for s in ok])
b, a0 = np.polyfit(np.log(N_), R_, 1)   # res = a0 + b ln N (b < 0)
nmin = a.nmin or int(np.clip(np.exp((a.target - a0) / b) if b < 0 else 5000, 1000, 50000))
print('periphery resolution vs particles over %d states: res = %.1f %+.2f ln N -> %d particles for %.1f A (median state %d particles, %.1f A)' % (len(ok), a0, b, nmin, a.target, int(np.median(N_)), float(np.median(R_))))
for F in sorted({v[2] for v in G.values()}):
    P = read_mrc(a.parent.replace('{F}', str(F)) + '.mrc').astype(np.float32); M = MK[F]
    Pl = lp(P); supp = Pl > 0.2 * Pl.max(); ss = [s for s in sorted(G) if G[s][2] == F and os.path.isfile('%s/recvol_state%02d.mrc' % (L1, s))]
    X = []
    for s in ss:
        v = lp(read_mrc('%s/recvol_state%02d.mrc' % (L1, s)).astype(np.float32)); A = np.stack([Pl[supp], np.ones(supp.sum())], 1); c = np.linalg.lstsq(A, v[supp], rcond=None)[0]
        X.append((((v - c[1]) / max(c[0], 1e-12) - Pl) * M)[M > 0.05])
    X = np.array(X); w = np.array([G[s][0] for s in ss], float); mu = (w[:, None] * X).sum(0) / w.sum(); Xc = (X - mu) * np.sqrt(w / w.sum())[:, None]
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False); ev = S ** 2 / (S ** 2).sum(); k = int(np.searchsorted(np.cumsum(ev), a.var) + 1); Z = (X - mu) @ Vt[:k].T
    print('family %d: %d states, %d components for %.0f %% of the variance (first three %s); sampling in the first %d (%s)' % (F, len(ss), k, 100 * a.var, ' '.join('%.2f' % e for e in ev[:3]), a.npc, a.mode))
    Z = Z[:, :a.npc]   # distances along the dominant motions only (the minor components are mostly noise)
    if a.mode == 'fps':   # seeds: the most extreme leaf, then repeatedly the leaf farthest from every seed so far
        K = max(1, int(w.sum() // nmin)); sd = [int(np.argmax(np.sum(Z ** 2, 1)))]; d = np.sum((Z - Z[sd[0]]) ** 2, 1)
        while len(sd) < min(K, len(ss)): j = int(np.argmax(d)); sd.append(j); d = np.minimum(d, np.sum((Z - Z[j]) ** 2, 1))
        lab = np.argmin(((Z[:, None, :] - Z[sd][None]) ** 2).sum(-1), 1); grp = {g: list(np.where(lab == g)[0]) for g in range(len(sd)) if (lab == g).any()}
    else: grp = {i: [i] for i in range(len(ss))}
    cen = {g: (w[m][:, None] * Z[m]).sum(0) / w[m].sum() for g, m in ((g, np.array(v)) for g, v in grp.items())}; cnt = {g: w[np.array(v)].sum() for g, v in grp.items()}
    while len(grp) > 1 and min(cnt.values()) < nmin:
        i = min(cnt, key=cnt.get); j = min((g for g in grp if g != i), key=lambda g: np.sum((cen[g] - cen[i]) ** 2))
        cen[j] = (cnt[i] * cen[i] + cnt[j] * cen[j]) / (cnt[i] + cnt[j]); cnt[j] += cnt[i]; grp[j] += grp.pop(i); cen.pop(i); cnt.pop(i)
    order = sorted(grp, key=lambda g: cen[g][0])   # along the first motion
    for g in order:
        gid += 1; mem = [ss[i] for i in grp[g]]
        for s in mem: out_lab[st == s] = gid
        lines.append('%d %d %s' % (gid, int(cnt[g]), '+'.join(G[s][1] for s in mem))); land.append('%d %d %s %s' % (gid, F, ' '.join('%.3f' % v for v in cen[g][:3]), ','.join(str(s) for s in mem)))
    print('  -> %d groups of >= %d particles; along the first motion: %s' % (len(grp), nmin, ' '.join(str(int(cnt[g])) for g in order)))
np.savetxt(a.out + '/states.txt', out_lab, fmt='%d')
with open(a.out + '/groups.txt', 'w') as f: f.write('# sampled state: particles, members (source states); budget %d particles\n' % nmin); f.write('\n'.join(lines) + '\n')
with open(a.out + '/landscape.txt', 'w') as f: f.write('# group family pc1 pc2 pc3 source_states\n'); f.write('\n'.join(land) + '\n')
print('%d sampled states from %d source states' % (gid, len(G)))
