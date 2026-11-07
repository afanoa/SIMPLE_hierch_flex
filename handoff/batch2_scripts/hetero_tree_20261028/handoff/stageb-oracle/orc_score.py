#!/usr/bin/env python3
"""orc_score.py PREP OUT ARM...: refine=state_eval labels of the oracle arms against the truth. ARM = name:dir:ident, dir holding
after.txt (the arm's ptcl3D dump); ident says how a state names its class: maj (its majority non-junk class), truth or coarse (the
label set the arm's references were built from, prep's truth_map.txt / lab_coarse.txt), or tmap=FILE (rows 'label class n').
Columns: NMI with the truth classes over all particles (junk a class) and over the non-junk ones; ARI; class accuracy (non-junk in
their class's state); group accuracy (non-junk in the state of their coarse group = the stage-B state their class mostly fills);
best-state F1 of junk, closed (18) and holo (101); label agreement with stage B. Then, for the non-junk particles stage B gets
wrong, what each arm does with them: 'merged' = in their class's majority state, which another class owns (a K limit);
'misassigned' = elsewhere."""
import os
import sys
import numpy as np

prep, fout = sys.argv[1:3]
arms = [a.split(':') for a in sys.argv[3:]]
info = np.loadtxt(os.path.join(prep, 'ptcl_info.txt'), ndmin=2)
combo = info[:, 0].astype(int); junk = info[:, 1] > 0; sB = info[:, 5].astype(int)
tmap = np.loadtxt(os.path.join(prep, 'truth_map.txt'), dtype=int, ndmin=2)
classes = tmap[:, 1]
lab_c = np.loadtxt(os.path.join(prep, 'lab_coarse.txt'), dtype=int)
act = sB > 0
nj = act & ~junk
fo = open(fout, 'w')


def P(s=''):
    fo.write(s + '\n'); print(s, flush=True)


def contingency(a, b):
    a = np.unique(a, return_inverse=True)[1]; b = np.unique(b, return_inverse=True)[1]
    C = np.zeros((a.max() + 1, b.max() + 1)); np.add.at(C, (a, b), 1)
    return C


def nmi(a, b):
    C = contingency(a, b); Pm = C / C.sum(); pa, pb = Pm.sum(1), Pm.sum(0); nz = Pm > 0
    mi = (Pm[nz] * np.log(Pm[nz] / np.outer(pa, pb)[nz])).sum()
    ha = -(pa[pa > 0] * np.log(pa[pa > 0])).sum(); hb = -(pb[pb > 0] * np.log(pb[pb > 0])).sum()
    return float(mi / max(np.sqrt(ha * hb), 1e-12))


def ari(a, b):
    C = contingency(a, b); c2 = lambda x: x * (x - 1) / 2.
    sij, sa, sb, nn = c2(C).sum(), c2(C.sum(1)).sum(), c2(C.sum(0)).sum(), c2(C.sum())
    ex = sa * sb / nn
    return float((sij - ex) / (0.5 * (sa + sb) - ex))


def states(fn):
    return np.array([int(float(dict(t.split('=', 1) for t in ln.split() if '=' in t).get('state', '0')))
                     for ln in open(fn) if 'e1=' in ln])


def majority(lab, target, mask):
    m = {}
    for u in np.unique(lab[mask]):
        v = target[mask & (lab == u)]
        m[u] = np.bincount(v - v.min()).argmax() + v.min() if len(v) else -1
    return m


def class_right(lab, ident):
    if ident == 'truth' or ident.startswith('tmap='):
        cl = classes if ident == 'truth' else np.loadtxt(ident[5:], dtype=int, ndmin=2)[:, 1]
        cls = np.where(lab > 0, cl[np.maximum(lab, 1) - 1], -1)
    else:
        mj = majority(lab, combo, nj)
        cls = np.array([mj.get(u, -1) for u in lab])
    return cls == combo


def group_right(lab, ident):
    if ident == 'coarse':
        return lab == lab_c
    mg = majority(lab, lab_c, nj)
    return np.array([mg.get(u, -1) for u in lab]) == lab_c


def f1(lab, c):
    mc = act & (combo == c)
    return max(2 * (mc & (lab == u)).sum() / (mc.sum() + (act & (lab == u)).sum()) for u in np.unique(lab[act]))


# stage B's errors on the non-junk particles
rB = class_right(sB, 'maj')
majB = majority(sB, combo, nj)
own = {c: int(np.bincount(sB[nj & (combo == c)]).argmax()) for c in np.unique(combo[nj])}
merged = nj & ~rB & (sB == np.array([own.get(c, -1) for c in combo]))
misas = nj & ~rB & ~merged
P('# %s: %d particles, %d non-junk; stage B (K %d) class-right %.3f, merged %.3f, misassigned %.3f of non-junk'
  % (prep, act.sum(), nj.sum(), sB.max(), rB[nj].mean(), merged[nj].mean(), misas[nj].mean()))
P('# stage-B state owners (majority non-junk class): %s; classes without a state: %s'
  % ({int(u): int(c) for u, c in majB.items()}, sorted(set(np.unique(combo[nj]).tolist()) - set(majB.values()))))
P('# ceiling at stage B\'s grouping (the truth merged into stage B\'s states, perfectly assigned): NMIall %.4f NMInj %.4f, K %d'
  % (nmi(lab_c[act], combo[act]), nmi(lab_c[nj], combo[nj]), len(np.unique(lab_c[act]))))
P('%-10s %3s %7s %7s %7s %7s %7s %7s %7s %7s %7s | %s' % ('arm', 'K', 'NMIall', 'NMInj', 'ARI', 'accCls', 'accGrp',
  'F1junk', 'F1cl18', 'F1h101', 'agreeB', 'stage-B wrong non-junk now class-right: merged / misassigned; B-right kept'))
rows = []
for name, d, ident in arms:
    fn = os.path.join(d, 'after.txt')
    if not os.path.exists(fn):
        P('%-10s missing %s' % (name, fn)); continue
    lab = states(fn)
    assert len(lab) == len(sB)
    cr = class_right(lab, ident); gr = group_right(lab, ident)
    agree = (lab[act] == sB[act]).mean() if (ident == 'maj' and lab.max() == sB.max()) else float('nan')
    f = {c: f1(lab, c) for c in np.unique(combo[act])}
    line = '%-10s %3d %7.4f %7.4f %7.4f %7.4f %7.4f %7.3f %7.3f %7.3f %7.4f | %.3f / %.3f; %.3f' % (
        name, len(np.unique(lab[act])), nmi(lab[act], combo[act]), nmi(lab[nj], combo[nj]), ari(lab[act], combo[act]),
        cr[nj].mean(), gr[nj].mean(), f.get(0, np.nan), f.get(18, np.nan), f.get(101, np.nan), agree,
        cr[merged].mean() if merged.any() else np.nan, cr[misas].mean() if misas.any() else np.nan, cr[nj & rB].mean())
    P(line)
    rows.append((name, lab, f))
P()
P('# best-state F1 per truth class')
for name, lab, f in rows:
    P('%-10s %s' % (name, ' '.join('%d:%.3f' % (c, v) for c, v in sorted(f.items()))))
P()
snrq = info[:, 3]
qs = np.percentile(snrq[nj], [0, 20, 40, 60, 80, 100])
P('# class accuracy / group accuracy of the non-junk particles per snrq quintile (edges %s)' % ' '.join('%.2f' % v for v in qs))
for name, d, ident in arms:
    fn = os.path.join(d, 'after.txt')
    if not os.path.exists(fn):
        continue
    lab = states(fn); cr = class_right(lab, ident); gr = group_right(lab, ident)
    cells = []
    for lo, hi in zip(qs[:-1], qs[1:]):
        m = nj & (snrq >= lo) & (snrq <= hi)
        cells.append('%.3f/%.3f' % (cr[m].mean(), gr[m].mean()))
    P('%-10s %s' % (name, '  '.join(cells)))
P()
P('# where each truth class goes (state: share), largest two states')
for name, lab, f in rows:
    parts = []
    for c in np.unique(combo[act]):
        m = act & (combo == c); cnt = np.bincount(lab[m]); o = np.argsort(-cnt)[:2]
        parts.append('%d->%s' % (c, ','.join('%d:%.2f' % (u, cnt[u] / m.sum()) for u in o if cnt[u] > 0)))
    P('%-10s %s' % (name, ' '.join(parts)))
fo.close()
