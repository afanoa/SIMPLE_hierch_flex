#!/usr/bin/env python3
"""profile.py SACCT_DUMP: where the 692k tree's wall and core-hours went, from the step jobs' accounting (WorkDir = node/step)."""
import sys, re, collections, datetime as dt
P = lambda s: dt.datetime.strptime(s, '%Y-%m-%dT%H:%M:%S')
jobs = []
for l in open(sys.argv[1]):
    r = l.rstrip('\n').split('|')
    if len(r) < 9 or 'f692_r17d_s1/1_hetero_tree' not in r[8]: continue
    tail = r[8].split('1_hetero_tree/')[1]
    parts = tail.split('/')
    node = parts[0] if parts[0].startswith('node_') else 'end'
    step = parts[1] if len(parts) > 1 and node != 'end' else parts[0]
    step = re.sub(r'^fsite_.*', 'fsite', step)
    depth = len(node) - len('node_r') if node != 'end' else -1
    jobs.append(dict(node=node, step=step, depth=depth, cpus=int(r[3]), ch=int(r[4]) / 3600, s=P(r[5]), e=P(r[6])))
t0 = min(j['s'] for j in jobs); t1 = max(j['e'] for j in jobs)
print('run span %s -> %s = %.2f h; %d step jobs; %.0f core-h' % (t0, t1, (t1 - t0).total_seconds() / 3600, len(jobs), sum(j['ch'] for j in jobs)))
# core-h and node-span by step
by = collections.defaultdict(lambda: [0., 0, set()])
for j in jobs:
    b = by[j['step']]; b[0] += j['ch']; b[1] += 1; b[2].add(j['node'])
# wall each step occupies on its node (first start to last end per node, summed over nodes, and the mean per node)
span = collections.defaultdict(lambda: collections.defaultdict(lambda: [None, None]))
for j in jobs:
    s = span[j['step']][j['node']]
    s[0] = j['s'] if s[0] is None or j['s'] < s[0] else s[0]; s[1] = j['e'] if s[1] is None or j['e'] > s[1] else s[1]
print('\nstep        core-h  share  jobs  nodes  node-span mean (min)  max (min)')
tot = sum(j['ch'] for j in jobs)
for st, (ch, n, nodes) in sorted(by.items(), key=lambda x: -x[1][0]):
    sp = [(v[1] - v[0]).total_seconds() / 60 for v in span[st].values()]
    print('%-10s %7.1f  %5.1f%%  %4d  %5d  %8.1f  %8.1f' % (st, ch, 100 * ch / tot, n, len(nodes), sum(sp) / len(sp), max(sp)))
# per depth
print('\ndepth  nodes  core-h  first start -> last end (h after t0)')
bd = collections.defaultdict(lambda: [set(), 0., None, None])
for j in jobs:
    b = bd[j['depth']]; b[0].add(j['node']); b[1] += j['ch']
    b[2] = j['s'] if b[2] is None or j['s'] < b[2] else b[2]; b[3] = j['e'] if b[3] is None or j['e'] > b[3] else b[3]
for d in sorted(bd):
    b = bd[d]
    print('%5s  %5d  %6.0f  %5.2f -> %5.2f' % ('end' if d < 0 else d, len(b[0]), b[1], (b[2] - t0).total_seconds() / 3600, (b[3] - t0).total_seconds() / 3600))
# per-node wall: the root and the closed lineage
print('\nnode wall (first job start -> last job end), the root, depth-1 nodes and the closed lineage:')
nd = collections.defaultdict(lambda: [None, None, 0.])
for j in jobs:
    b = nd[j['node']]; b[0] = j['s'] if b[0] is None or j['s'] < b[0] else b[0]; b[1] = j['e'] if b[1] is None or j['e'] > b[1] else b[1]; b[2] += j['ch']
for n in ['node_r', 'node_r1', 'node_r2', 'node_r3', 'node_r4', 'node_r13', 'node_r131', 'node_r1311']:
    if n in nd: b = nd[n]; print('  %-10s %6.1f min  %6.0f core-h  (%.2f -> %.2f h)' % (n, (b[1] - b[0]).total_seconds() / 60, b[2], (b[0] - t0).total_seconds() / 3600, (b[1] - t0).total_seconds() / 3600))
# concurrency timeline: cores busy per 10-min bin
bins = collections.defaultdict(float)
for j in jobs:
    t = j['s']
    while t < j['e']:
        nb = min(j['e'], (t + dt.timedelta(minutes=10)).replace(minute=(t.minute // 10) * 10, second=0, microsecond=0) + dt.timedelta(minutes=10) if False else t + dt.timedelta(minutes=10))
        bins[int((t - t0).total_seconds() // 600)] += j['cpus'] * (nb - t).total_seconds() / 600
        t = nb
print('\ncores busy per 10-min bin (hours after start):')
ks = sorted(bins); line = ''
for k in ks:
    line += '%4.0f' % bins[k]
    if (k + 1) % 12 == 0: print('  h%2d-%2d:%s' % (k // 12 * 2, k // 12 * 2 + 2, line)); line = ''
if line: print('  h%2d+  :%s' % (ks[-1] // 12 * 2, line))
idle = [k for k in ks if bins[k] < 40]
print('bins under 40 cores busy: %d of %d (%.0f min)' % (len(idle), len(ks), 10 * len(idle)))
