#!/usr/bin/env python3
"""Does the tree keep a rare population together? Tracer = rows of chosen states of an earlier run (e.g. f692r's closed states). For every
node of an rflex run (queue.txt / level*.txt: node depth parent_project rows): tracer particles in the node and their spread over the
node's refined states (share of the largest state, effective number of states exp(entropy)), and over the node's outcome (families /
leaves from groups.txt). A split that keeps the tracer in one state / one child keeps it together.
usage: tree_tracer.py RUNDIR TRACER_UNITS IDS"""
import sys, os, glob, numpy as np
R, tf, ids = sys.argv[1], sys.argv[2], [int(x) for x in sys.argv[3].split(',')]
tr = np.isin(np.loadtxt(tf, dtype=int), ids); print('tracer: states %s of %s, %d particles' % (sys.argv[3], tf, tr.sum()))
lines = []
for q in [R + '/queue.txt'] + sorted(glob.glob(R + '/level*.txt')):
    if os.path.isfile(q): lines += [l.split() for l in open(q) if l.strip()]
seen = set()
for l in lines:
    u, dep, rows = l[0], int(l[1]), l[3]
    if u in seen or not os.path.isfile(rows): continue
    seen.add(u); D = '%s/node_%s' % (R, u); m = np.loadtxt(rows, dtype=int) > 0; t = tr & m; nt = t.sum()
    out = '%s%s (depth %d): %d particles, tracer %d (%.0f %% of the tracer, %.1f %% of the node)' % ('  ' * dep, u, dep, m.sum(), nt, 100 * nt / max(tr.sum(), 1), 100 * nt / max(m.sum(), 1))
    if os.path.isfile(D + '/refined_states.txt') and nt:
        st = np.loadtxt(D + '/refined_states.txt', dtype=int); c = np.bincount(st[t], minlength=st.max() + 1)[1:]; p = c / c.sum(); p = p[p > 0]
        top = np.argsort(c)[::-1][:4]
        out += ' | states: ' + ' '.join('s%d %.0f%%' % (s + 1, 100 * c[s] / c.sum()) for s in top if c[s]) + ' | largest %.2f, effective number %.1f' % (c.max() / c.sum(), np.exp(-(p * np.log(p)).sum()))
        if os.path.isfile(D + '/groups.txt'):
            parts = []
            for g in open(D + '/groups.txt'):
                w = g.split();
                if len(w) < 3: continue
                mem = [int(x) for x in w[2].split('+')]; parts.append('%s[%s] %.0f%%' % (w[0], w[2], 100 * np.isin(st[t], mem).mean()))
            out += ' | outcome: ' + ', '.join(parts)
    print(out)
