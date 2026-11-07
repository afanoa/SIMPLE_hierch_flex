#!/usr/bin/env python3
"""usage: split_maps.py TREE_DIR > MAPLIST
fa_maps10.py map list of every split's stage-B maps: per node with a 'split group' line in the run's tree.log, each stage-B state that
dominates one of its groups, named split_<node>_s<state> (kind other, family split, particles = that group's size), with its halves."""
import os, re, sys
T = sys.argv[1].rstrip('/'); SP = {}
for l in open(os.path.join(os.path.dirname(T), 'tree.log'), errors='replace'):
    m = re.search(r'HETERO_TREE node (\S+): split group (\d+) \((\d+) particles; stage-B states ([^)]*)\)', l)
    if m: st = {int(p.split(':')[0]): int(p.split(':')[1]) for p in m.group(4).split()}; SP.setdefault(m.group(1), {})[max(st, key=st.get)] = int(m.group(3))
for X in sorted(SP, key=lambda s: (len(s), s)):
    for s, n in sorted(SP[X].items()):
        f = ['%s/node_%s/stageB/recvol_state%02d%s.mrc' % (T, X, s, x) for x in ('', '_even', '_odd')]
        if all(os.path.exists(p) for p in f): print('split_%s_s%d other split %d %s' % (X, s, n, ' '.join(f)))
