#!/usr/bin/env python3
"""usage: site_score.py SEED_DIR NODE_DIR OUT   (one candidate site of a node as a score; null_score3.py's pole, for any site)
The site's non-rest poles are its site_poles.txt lines that are not the rest; with several, the one holding the most site seeds. The
score is the mean of that pole's endpoint posteriors (site_post_<endpoint>.txt), one value per project row, high = the pole, -1
outside the node and on the node's junk rows (flag 1 in the seed dir's site_junk.txt). site_pole_rule.txt is not read."""
import os, sys, numpy as np
SD, ND, OUT = (a.rstrip('/') for a in sys.argv[1:4])
rows = np.loadtxt(ND + '/rows.txt', dtype=int).ravel() - 1
real = {int(w[0]): w[2:] for w in (l.split() for l in open(SD + '/site_poles.txt')) if len(w) >= 3 and w[2:] != ['rest']}
if not real: sys.exit('%s: no non-rest pole' % SD)
ss = np.loadtxt(SD + '/site_seeds.txt', dtype=int); NM = len(ss); j = max(real, key=lambda q: (ss == q).sum()); eps = real[j]
P = np.mean([np.loadtxt('%s/site_post_%s.txt' % (SD, e)) for e in eps], 0); s = np.full(NM, -1.); s[rows] = P[rows]; s[s < 0] = -1.
nj = 0
if os.path.exists(SD + '/site_junk.txt'):
    J = np.loadtxt(SD + '/site_junk.txt', comments='#', ndmin=2); jr = J[J[:, 1] == 1, 0].astype(int) - 1; s[jr] = -1.; nj = len(jr)
np.savetxt(OUT, s, fmt='%.6f')
print('%s: %d node rows, %d junk rows out; %d non-rest pole(s) %s; scored pole state %d (endpoints %s, %d site seeds)' % (os.path.basename(SD), len(rows), nj, len(real),
      ' '.join('%d:%s' % (q, '+'.join(real[q])) for q in sorted(real)), j, ' '.join(eps), (ss == j).sum()))
