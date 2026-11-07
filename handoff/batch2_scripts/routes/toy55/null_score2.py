#!/usr/bin/env python3
"""usage: null_score2.py NODE_DIR OUT   (null_score.py; a seed dir without site_pole_rule.txt keeps no pole, the fflex fallback)
The score along which a tree node was split, one value per project row (-1 outside the node):
- seeded by a site pole (the seed dir's site_pole_rule.txt keeps one): the mean of that pole's endpoint posteriors
  (site_post_<endpoint>.txt), high = the pole;
- the fflex fallback (no pole kept; seeds/seeds.txt holds the two groups of the node's fflex states): the Fisher discriminant of the
  node's fflex z between the two groups, as a percentile rank over the node, high = group 1.
The seed dir is the node's site/ or, when seed_site.txt names site k > 1, nsite_<node>_s<k> (fsite_ in older trees)."""
import os, sys, numpy as np
ND, OUT = sys.argv[1:3]; ND = ND.rstrip('/'); nid = os.path.basename(ND)[5:]
rows = np.loadtxt(ND + '/rows.txt', dtype=int).ravel() - 1; seeds = np.loadtxt(ND + '/seeds/seeds.txt', dtype=int); NM = len(seeds)
f = ND + '/sites/seed_site.txt'; k = int(open(f).read().split()[0]) if os.path.exists(f) and open(f).read().strip() else 1
sd = ND + '/site' if k <= 1 else next(d for d in (ND + '/nsite_%s_s%d' % (nid, k), ND + '/fsite_%s_s%d' % (nid, k)) if os.path.isdir(d))
rule = [l.split() for l in open(sd + '/site_pole_rule.txt') if l.strip() and not l.startswith('#')] if os.path.exists(sd + '/site_pole_rule.txt') else []; kept = [int(w[0]) for w in rule if w[-1] == 'kept']
poles = {int(w[0]): w[2:] for w in [l.split() for l in open(sd + '/site_poles.txt') if l.strip()]}
s = np.full(NM, -1.); inn = np.zeros(NM, bool); inn[rows] = True
if kept:
    eps = poles[kept[0]]; P = np.mean([np.loadtxt('%s/site_post_%s.txt' % (sd, e)) for e in eps], 0); s[inn] = P[inn]
    print('%s: %d rows; seed dir %s; pole %d kept (endpoints %s, %d seeds): score = the mean of its endpoint posteriors' % (nid, inn.sum(), os.path.basename(sd), kept[0], ' '.join(eps), (seeds == 1).sum()))
else:
    b = open(ND + '/fflex/flex_pca_embedding.bin', 'rb').read(); ver, n, nc, box = (int(v) for v in np.frombuffer(b[8:24], '<i4'))
    pi = np.frombuffer(b[28:28 + 4 * n], '<i4').astype(int) - 1; z = np.frombuffer(b[28 + 4 * n:28 + 4 * n + 8 * n * nc], '<f8').reshape(nc, n).T
    g1, g2 = seeds[pi] == 1, seeds[pi] == 2; mu = z[g1].mean(0) - z[g2].mean(0); Sw = np.cov(z[g1].T) * g1.sum() + np.cov(z[g2].T) * g2.sum()
    wv = np.linalg.solve(Sw / (g1.sum() + g2.sum()) + 1e-9 * np.eye(nc), mu); x = z @ wv; s[pi] = (np.argsort(np.argsort(x)) + 0.5) / n
    auc = (np.argsort(np.argsort(x[g1 | g2]))[g1[g1 | g2]] + 1.).sum(); n1, n2 = g1.sum(), g2.sum(); auc = (auc - n1 * (n1 + 1) / 2) / (n1 * n2)
    print('%s: %d rows; seed dir %s keeps no pole: fflex fallback, groups %d / %d; score = Fisher discriminant of %d fflex z (group AUC %.3f)' % (nid, inn.sum(), os.path.basename(sd), n1, n2, nc, auc))
np.savetxt(OUT, s, fmt='%.6f')
