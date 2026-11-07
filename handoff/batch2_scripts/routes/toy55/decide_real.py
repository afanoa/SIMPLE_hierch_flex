#!/usr/bin/env python3
"""decide_real.py = rflex6split_decide.py (split-decision session, copied 2026-10-21) for REAL half-stacks: the test's children files live in
each half's row space and are mapped to master rows through MRA / MRB (half row -> master row, 1-based).
rflex6split node decision from the split test (splitbench11 output in TESTDIR): refined states the junk screen set aside (junk_units.txt)
are junk leaves; every other refined state joins the supported group (child of the test) holding most of its particles. Per group:
family (>= MINP particles, depth < MAXD, at least two supported groups: the node splits) or leaves; one supported group = the node is
one leaf. Prints one line per group: kind particles members(+) note
usage: decide_real.py STATES TESTDIR DEPTH MRA MRB [--minp 50000 --maxd 2]"""
import argparse, os, re, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('states'); ap.add_argument('test'); ap.add_argument('depth', type=int); ap.add_argument('mra'); ap.add_argument('mrb')
ap.add_argument('--minp', type=int, default=50000); ap.add_argument('--maxd', type=int, default=2); a = ap.parse_args()
lab = np.loadtxt(a.states, dtype=int); ch = np.zeros(len(lab), int)
for h, mr in (('A', a.mra), ('B', a.mrb)): ch[np.loadtxt(mr).astype(int) - 1] = np.loadtxt(a.test + '/children_%s.txt' % h, dtype=int)
jf = a.test + '/junk_units.txt'; junk = [int(x) for x in open(jf).read().split()] if os.path.isfile(jf) else []
jm = {}
for l in open(a.test + '/bench.log', errors='replace'):
    m = re.match(r'\s*junk screen unit (\d+): projection on the screened consensus (\S+) / (\S+), rel. SSNR (\S+) / (\S+)', l)
    if m: jm[int(m.group(1))] = 'projection %s/%s rel.SSNR %s/%s' % m.groups()[1:]
to = {}
for s in np.unique(lab[lab > 0]):
    if int(s) in junk: continue
    c = ch[lab == s]; c = c[c > 0]; to[int(s)] = int(np.bincount(c).argmax()) if len(c) else 0
big = max(set(to.values()) - {0}, key=lambda g: int((ch == g).sum()), default=1)
to = {s: (g if g > 0 else big) for s, g in to.items()}   # a state with no particles in any child joins the largest group
groups = sorted(set(to.values())); split = len(groups) >= 2
for g in groups:
    mem = sorted(s for s in to if to[s] == g); n = int(np.isin(lab, mem).sum())
    print('family' if (split and n >= a.minp and a.depth < a.maxd) else 'leaves', n, '+'.join(map(str, mem)), 'supported group %d of %d' % (g, len(groups)))
for s in junk: print('junk', int((lab == s).sum()), s, 'junk screen: ' + jm.get(s, ''))
