#!/usr/bin/env python3
"""Partitions of a coordinate node for the held-out comparison of its leaf models (problem 1): hard = the node's refined states,
bins = NB equal-width bins of the motion coordinate (P .. 100-P percentiles, ends clipped in), rand = random labels 1..NB (a model
without heterogeneity); all on the common rows (hard > 0, finite coordinate). Also the scoring passes' state files: random labels
over all 3*NB states (score_states.txt) and over the NB hard states (score_states_hard.txt) on the common rows, so every state is
populated (prob_neigh geom scores every state in the same pose neighbourhood, the input state only seeds nothing), and eo per row
from a print_project_field dump.
usage: p1_parts.py OUT --hard refined_states.txt --coord motion_coord.txt --dump base_dump.txt [--nb 8 --p 1 --seed 1]"""
import re, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--hard', required=True); ap.add_argument('--coord', required=True)
ap.add_argument('--dump', required=True); ap.add_argument('--nb', type=int, default=8); ap.add_argument('--p', type=float, default=1.); ap.add_argument('--seed', type=int, default=1)
a = ap.parse_args(); rng = np.random.default_rng(a.seed)
hard = np.loadtxt(a.hard, dtype=int); c = np.loadtxt(a.coord); eo = np.array([int(float(re.search(r' eo=([-0-9.eE+]+)', ' ' + l).group(1))) for l in open(a.dump) if l.strip()])
assert len(hard) == len(c) == len(eo), (len(hard), len(c), len(eo))
com = (hard > 0) & np.isfinite(c); nb = a.nb; z = np.zeros(len(hard), int)
lo, hi = np.percentile(c[com], [a.p, 100 - a.p]); e = np.linspace(lo, hi, nb + 1); bins = z.copy(); bins[com] = np.clip(np.searchsorted(e, c[com]) - 1, 0, nb - 1) + 1
hd = np.where(com, hard, 0); rnd = z.copy(); rnd[com] = rng.integers(1, nb + 1, com.sum())
s3 = z.copy(); s3[com] = rng.integers(1, 3 * nb + 1, com.sum()); s1 = z.copy(); s1[com] = rng.integers(1, nb + 1, com.sum())
for nm, x in (('hard', hd), ('bins', bins), ('rand', rnd), ('score_states', s3), ('score_states_hard', s1), ('eo', eo), ('common', com.astype(int))): np.savetxt('%s/%s.txt' % (a.out, nm), x, fmt='%d')
print('%d common rows (even %d, odd %d); hard %s | bins %s | rand %s' % (com.sum(), (com & (eo == 0)).sum(), (com & (eo == 1)).sum(),
      ' '.join(str((hd == k).sum()) for k in range(1, hd.max() + 1)), ' '.join(str((bins == k).sum()) for k in range(1, nb + 1)), ' '.join(str((rnd == k).sum()) for k in range(1, nb + 1))))
