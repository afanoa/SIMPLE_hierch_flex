#!/usr/bin/env python3
"""Per-conformer recovery readout from a deposited_match table over the 55 docked apo maps: a conformer is RECOVERED when
some state matches it with deviation correlation >= CMIN and a margin >= MARGIN over that state's second-best conformer.
usage: recover55.py MATCHES.txt [--cmin 0.4] [--margin 0.05] [--names map_order.txt]"""
import sys, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('matches'); ap.add_argument('--cmin', type=float, default=0.4); ap.add_argument('--margin', type=float, default=0.05); ap.add_argument('--names', default=''); a = ap.parse_args()
hdr = open(a.matches).readline(); ids = [int(x) for x in hdr.split('|')[1].split()[3:]]
rows = [l.split() for l in open(a.matches) if l.strip() and not l.startswith('#')]; C = np.array([[float(x) for x in r[3:]] for r in rows]); names = [l.strip() for l in open(a.names)] if a.names else ['state %s' % r[0] for r in rows]
rec = {}
for s in range(len(rows)):
    o = np.argsort(-C[s]); c1, c2 = C[s, o[0]], C[s, o[1]]
    if c1 >= a.cmin and c1 - c2 >= a.margin:
        k = ids[o[0]]
        if k not in rec or c1 > rec[k][0]: rec[k] = (c1, c1 - c2, names[s])
print('recovered %d of %d deposited conformers (corr >= %.2f, margin >= %.2f): %s' % (len(rec), len(ids), a.cmin, a.margin, ' '.join(str(k) for k in sorted(rec))))
for k in sorted(rec): print('  conformer %2d: corr %.2f margin %.2f  by %s' % (k, rec[k][0], rec[k][1], rec[k][2]))
