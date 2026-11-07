#!/usr/bin/env python3
"""decide_all.py: over-split node decision (no merging at the node; merging is left to the final merge-up). Same arguments and output as
decide_real.py: every refined state the split test's junk screen kept is its own group, a family when it holds >= MINP particles and
depth < MAXD, else its own leaf; junk-screen states are junk leaves. Used because splitbench15's post-EM re-test in the beyond-EM band
cannot see low-resolution differences (f692u7: a 0.87 closed-arm state of 8.8k merged into one leaf at 0.56).
usage: decide_all.py STATES TESTDIR DEPTH MRA MRB [--minp 15000 --maxd 4]"""
import argparse, os, re, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('states'); ap.add_argument('test'); ap.add_argument('depth', type=int); ap.add_argument('mra'); ap.add_argument('mrb')
ap.add_argument('--minp', type=int, default=15000); ap.add_argument('--maxd', type=int, default=4); a = ap.parse_args()
lab = np.loadtxt(a.states, dtype=int)
jf = a.test + '/junk_units.txt'; junk = [int(x) for x in open(jf).read().split()] if os.path.isfile(jf) else []
jm = {}
if os.path.isfile(a.test + '/bench.log'):
    for l in open(a.test + '/bench.log', errors='replace'):
        m = re.match(r'\s*junk screen unit (\d+): projection on the screened consensus (\S+) / (\S+), rel. SSNR (\S+) / (\S+)', l)
        if m: jm[int(m.group(1))] = 'projection %s/%s rel.SSNR %s/%s' % m.groups()[1:]
for s in np.unique(lab[lab > 0]):
    s = int(s)
    if s in junk: continue
    n = int((lab == s).sum())
    print('family' if (n >= a.minp and a.depth < a.maxd) else 'leaves', n, s, 'refined state %d (no node merge)' % s)
for s in junk: print('junk', int((lab == s).sum()), s, 'junk screen: ' + jm.get(s, ''))
