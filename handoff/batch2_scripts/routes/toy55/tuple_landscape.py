#!/usr/bin/env python3
"""Factorised landscape of a run: every particle carries the tuple of its region states (argmax-cell membership of each
region board of its family, 0 = no supported state), a conformer = an occupied tuple. Prints per region the state counts
and per family the occupied tuples above --minp with their populations. usage: tuple_landscape.py RUN NROWS [--minp 1000]"""
import sys, os, glob, re, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('run'); ap.add_argument('nrows', type=int); ap.add_argument('--minp', type=int, default=1000); a = ap.parse_args()
for fam in (1, 2):
    regs = sorted(glob.glob('%s/fam%d/region*/level1/clusters.txt' % (a.run, fam))); regs = [r for r in regs if '/zoomregion' not in r]
    labs, names = [], []
    for cf in regs:
        L1 = os.path.dirname(cf); lat = np.loadtxt(L1 + '/../board/klem_latent.txt', comments='#'); pind, am = lat[:, 0].astype(int), lat[:, 3].astype(int); c2s = {}
        for l in open(cf):
            if l.startswith('#') or not l.strip(): continue
            t = l.split(); [c2s.__setitem__(int(c), int(t[0])) for c in t[4:]]
        lab = np.zeros(a.nrows, int); lab[pind - 1] = [c2s.get(c, 0) for c in am]; labs.append(lab); names.append(re.search(r'fam\d/(region[^/]+)/', cf).group(1))
    if not labs: continue
    T = np.stack(labs, 1); inpart = (T > 0).any(1); print('family %d: %d particles in %d regions: %s' % (fam, int(inpart.sum()), len(names), ' '.join('%s(%d states)' % (nm, len(set(T[:, j][T[:, j] > 0]))) for j, nm in enumerate(names))))
    keys, cnt = np.unique(T[inpart], axis=0, return_counts=True); order = np.argsort(-cnt); big = [(tuple(keys[i]), cnt[i]) for i in order if cnt[i] >= a.minp]
    print('  occupied tuples >= %d particles: %d (of %d occupied); top 15: %s' % (a.minp, len(big), len(cnt), '  '.join('%s:%d' % ('-'.join(str(x) for x in k), n) for k, n in big[:15])))
    full = [(k, n) for k, n in big if all(x > 0 for x in k)]; print('  tuples with a state in EVERY region: %d above %d particles, %d particles in total' % (len(full), a.minp, sum(n for _, n in full)))
