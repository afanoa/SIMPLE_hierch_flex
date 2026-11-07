#!/usr/bin/env python3
"""Audit of node f692st r1's partitions (92,554 particles). (1) Confounders: NMI of each partition (repro_r1/labels, e3 seeds) with
projection-direction bins (64 equal-area), defocus octiles, the eo half and alignment-score (corr) octiles, each next to a permutation null;
per-state defocus / eo / corr / direction spread of the script's stage B. (2) Decisiveness: cross-half log-likelihoods of every r1 particle
against the plain leaf maps (score_embed_f692st, fixed node poses, band to 15 A): own r1 leaf vs the best other r1 leaf (rank, margin in nats,
share decided by >= 3 nats), and for contrast the family decision (best r1 leaf vs best leaf of the other families).  usage: audit_r1.py"""
import os, numpy as np
K_ = '/mnt/beegfs/elmlund/afan/rnasep_klem'; L = K_ + '/repro_r1/labels'; R = K_ + '/rflex_f692st'; H = '/mnt/beegfs/elmlund/afan/routes/toy55/simrnasep692'
FL = np.loadtxt(R + '/final_labels.txt', dtype=int); NR = len(FL); rows = np.where((FL >= 1) & (FL <= 8))[0]
def nmi(a, b):
    ua, ia = np.unique(a, return_inverse=True); ub, ib = np.unique(b, return_inverse=True)
    C = np.zeros((len(ua), len(ub))); np.add.at(C, (ia, ib), 1); P = C / C.sum(); pa, pb = P.sum(1), P.sum(0); nz = P > 0
    I = (P[nz] * np.log(P[nz] / np.outer(pa, pb)[nz])).sum(); h = lambda p: -(p[p > 0] * np.log(p[p > 0])).sum(); return I / np.sqrt(h(pa) * h(pb))
e1, e2, df, eo, corr = (np.zeros(NR) for _ in range(5)); want = ('e1', 'e2', 'dfx', 'dfy', 'eo', 'corr')
for i, l in enumerate(open(R + '/final_plain/poses_combined.txt')):
    d = dict(t.split('=', 1) for t in l.split() if '=' in t and t.split('=', 1)[0] in want)
    e1[i], e2[i] = float(d.get('e1', 0)), float(d.get('e2', 0)); df[i] = 0.5 * (float(d.get('dfx', 0)) + float(d.get('dfy', 0)))
    eo[i], corr[i] = float(d.get('eo', 0)), float(d.get('corr', 0))
def octile(x): return np.searchsorted(np.quantile(x, np.linspace(0, 1, 9)[1:-1]), x)
ct = np.cos(np.radians(e2[rows])); tb = octile(ct); az = np.zeros(len(rows), int)
for t in range(8):
    m = tb == t; az[m] = octile(e1[rows][m])
CONF = {'direction (64)': tb * 8 + az, 'defocus (8)': octile(df[rows]), 'eo half': eo[rows].astype(int), 'align score (8)': octile(corr[rows])}
parts = [(n, os.path.join(L, n + '.txt')) for n in ('S_flex', 'S_fflex', 'S_A', 'S_B', 'O_fflex', 'O_A', 'P_flex', 'P_fflex', 'P_A', 'P_B')]
parts += [('E3_script_seeds', K_ + '/repro_r1/e3_script/seeds.txt'), ('E3_port1_seeds', K_ + '/repro_r1/e3_port1/seeds.txt')]
rng = np.random.default_rng(0)
print('== (1) confounders: NMI of each r1 partition with nuisance variables (permutation null in brackets)')
print('%-16s' % 'partition' + ''.join('%22s' % k for k in CONF))
for n, f in parts:
    if not os.path.exists(f): continue
    a = np.loadtxt(f, dtype=int)[rows]
    print('%-16s' % n + ''.join('%22s' % ('%.3f [%.3f]' % (nmi(a, v), np.mean([nmi(rng.permutation(a), v) for _ in range(3)]))) for v in CONF.values()))
a = np.loadtxt(L + '/S_B.txt', dtype=int)[rows]
print('-- script stage B per state: particles, mean defocus (um), eo share, mean align score, direction-bin entropy (node %.2f nats)' %
      -(np.bincount(CONF['direction (64)']) / len(rows) * np.log(np.maximum(np.bincount(CONF['direction (64)']) / len(rows), 1e-12))).sum())
for s in np.unique(a):
    m = a == s; p = np.bincount(CONF['direction (64)'][m], minlength=64) / m.sum()
    print('   state %d: %6d  %.3f  %.3f  %.4f  %.2f' % (s, m.sum(), df[rows][m].mean(), eo[rows][m].mean(), corr[rows][m].mean(), -(p[p > 0] * np.log(p[p > 0])).sum()))
print('== (2) decisiveness: cross-half log-likelihoods (nats) of each r1 particle against the 8 r1 plain leaf maps')
own_all, M8, Mfam = [], [], []
for h in 'AB':
    mr = np.loadtxt(H + '_%s/meta_rows.txt' % h).astype(int) - 1; S = np.load(K_ + '/score_embed_f692st/scores_%s.npy' % h).astype(np.float64)
    sel = np.where((FL[mr] >= 1) & (FL[mr] <= 8))[0]; S = S[sel]; own = FL[mr[sel]] - 1
    S8 = S[:, :8]; o = S8[np.arange(len(sel)), own]; oth = S8.copy(); oth[np.arange(len(sel)), own] = -np.inf
    own_all.append(own); M8.append(np.c_[o - oth.max(1), np.sort(S8, 1)[:, -1] - np.sort(S8, 1)[:, -2]]); Mfam.append(S8.max(1) - S[:, 8:].max(1))
own = np.concatenate(own_all); M8 = np.vstack(M8); Mfam = np.concatenate(Mfam)
q = lambda x: ' '.join('%.1f' % v for v in np.quantile(x, [0.1, 0.25, 0.5, 0.75, 0.9]))
print('   %d particles; own leaf best among the 8: %.1f %%; own minus best other (q10 q25 q50 q75 q90): %s' % (len(own), 100 * (M8[:, 0] > 0).mean(), q(M8[:, 0])))
print('   best minus second best (q10..q90): %s; share with a gap >= 3 nats: %.1f %%, >= 1 nat: %.1f %%' % (q(M8[:, 1]), 100 * (M8[:, 1] >= 3).mean(), 100 * (M8[:, 1] >= 1).mean()))
print('   family decision, best r1 leaf minus best leaf elsewhere (q10..q90): %s; share >= 3 nats: %.1f %%' % (q(Mfam), 100 * (Mfam >= 3).mean()))
names = ['r1_s5', 'r1_s1', 'r1_s4', 'r1_s3', 'r1_s2', 'r1_s6', 'r1_s7', 'r1_s8']
for g in range(8):
    m = own == g; print('   %-6s %6d: own best %.1f %%, median margin %.1f nats, gap >= 3 nats %.1f %%' % (names[g], m.sum(), 100 * (M8[m, 0] > 0).mean(), np.median(M8[m, 0]), 100 * (M8[m, 1] >= 3).mean()))
