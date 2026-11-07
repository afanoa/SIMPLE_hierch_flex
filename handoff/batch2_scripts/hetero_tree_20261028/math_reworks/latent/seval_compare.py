#!/usr/bin/env python3
"""seval_compare.py: refine=state_eval against the project it started from. Poses must not move (e1 e2 e3 x y); the states are
compared with the input labels (agreement, per-state transfer) and, given a qtoy truth, both labelings with the truth (NMI).
usage: seval_compare.py BEFORE_DUMP AFTER_DUMP [TRUTH_QTOY ROWS_TXT]"""
import sys
import numpy as np


def read(fn):
    return [dict(t.split('=', 1) for t in ln.split() if '=' in t) for ln in open(fn) if 'e1=' in ln]


def col(rows, k):
    return np.array([float(r.get(k, 'nan')) for r in rows])


def nmi(a, b):
    a = np.unique(a, return_inverse=True)[1]
    b = np.unique(b, return_inverse=True)[1]
    C = np.zeros((a.max() + 1, b.max() + 1))
    np.add.at(C, (a, b), 1)
    P = C / len(a)
    pa, pb = P.sum(1), P.sum(0)
    nz = P > 0
    mi = (P[nz] * np.log(P[nz] / np.outer(pa, pb)[nz])).sum()
    ha = -(pa[pa > 0] * np.log(pa[pa > 0])).sum()
    hb = -(pb[pb > 0] * np.log(pb[pb > 0])).sum()
    return float(mi / max(np.sqrt(ha * hb), 1e-12))


b, a = read(sys.argv[1]), read(sys.argv[2])
assert len(a) == len(b), (len(a), len(b))
dpose = max(float(np.nanmax(np.abs(col(a, k) - col(b, k)))) for k in ('e1', 'e2', 'e3', 'x', 'y'))
sb, sa = col(b, 'state').astype(int), col(a, 'state').astype(int)
act = sb > 0
print('particles %d (active %d); largest pose change %.2e; state agreement with the input %.4f'
      % (len(a), act.sum(), dpose, (sa[act] == sb[act]).mean()))
K = max(sb.max(), sa.max())
M = np.zeros((K + 1, K + 1), int)
np.add.at(M, (sb[act], sa[act]), 1)
for s in range(1, K + 1):
    print('  state %d: %6d in, %6d out, kept %.3f' % (s, M[s].sum(), M[:, s].sum(), M[s, s] / max(M[s].sum(), 1)))
if len(sys.argv) > 4:
    t = np.loadtxt(sys.argv[3], ndmin=2)
    rows = np.loadtxt(sys.argv[4], dtype=int, ndmin=1)
    combo = t[rows - 1, 1] if len(rows) == len(a) else t[:len(a), 1]
    print('NMI with the truth combo: input %.4f, state_eval %.4f' % (nmi(sb[act], combo[act]), nmi(sa[act], combo[act])))
