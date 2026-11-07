#!/usr/bin/env python3
"""Post-hoc readout of junk_em.py runs by junk type: each truth-junk particle is typed blob or broken by matching its CTF triple (dfx, dfy,
angast) in SIM/all.txt to SIM/oris_blob.txt / SIM/oris_broken.txt; per run, the junk-component calls (call > K) are counted per type and per
real conformer, and the junk components' own members are broken down the same way.
usage: junk_types.py SIM RUNDIR [RUNDIR ...]"""
import sys, re, numpy as np
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
def key(l): d = kv(l); return (round(float(d['dfx']), 5), round(float(d['dfy']), 5), round(float(d['angast']), 3))
sim = sys.argv[1]; typ = {}
for t in ('blob', 'broken'):
    for l in open('%s/oris_%s.txt' % (sim, t)):
        if l.strip(): typ[key(l)] = t
G = np.loadtxt(sim + '/gt_particles.txt', comments='#'); conf = np.full(int(G[:, 0].max()), -1); conf[G[:, 0].astype(int) - 1] = G[:, 1].astype(int)
lab = np.array(['c%d' % c for c in conf], dtype=object); rows = [key(l) for l in open(sim + '/all.txt') if l.strip()]
for i, k in enumerate(rows):
    if conf[i] == 0: lab[i] = typ.get(k, 'junk?')
cats = sorted(set(lab)); print('%s: %s' % (sim, ', '.join('%s %d' % (c, (lab == c).sum()) for c in cats)))
for d in sys.argv[2:]:
    L = open(d + '/junk_em.log').read(); m = re.search(r'components: states 1\.\.(\d+), free ([\d,-]+), noise ([\d-]+)', L)
    if not m: print('%s: no components line' % d); continue
    K = int(m.group(1)); fr = [int(x) for x in m.group(2).split(',') if x != '-']; no = [int(m.group(3))] if m.group(3) != '-' else []
    call = np.maximum(np.loadtxt(d + '/calls_A.txt', dtype=int), np.loadtxt(d + '/calls_B.txt', dtype=int)); jc = call > K
    s = '  '.join('%s %.2f' % (c, jc[lab == c].mean()) for c in cats)
    comp = '  '.join('%s%d: %s' % ('free' if g in fr else 'noise', g, ' '.join('%s %d' % (c, ((call == g) & (lab == c)).sum()) for c in cats)) for g in fr + no)
    print('%s: junk-call rate per category: %s | %s' % (d.rstrip('/').split('/')[-1], s, comp))
