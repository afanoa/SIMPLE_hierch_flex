#!/usr/bin/env python3
"""Truth purity of consolidated (merged) states of a synthetic run: per merged state the majority conformer, its purity
and size; summary = states, particles per state, mean purity, conformers covered (majority of >= 1 state).
usage: consolidated_purity.py CONS_DIR GT_PARTICLES"""
import sys, numpy as np
cons, gtf = sys.argv[1:3]; st = np.loadtxt(cons + '/states.txt', dtype=int)
conf = {int(w[0]): int(w[1]) for w in (l.split() for l in open(gtf)) if w and not w[0].startswith('#')}; g = np.array([conf.get(i + 1, 0) for i in range(len(st))])
K = st.max(); rows = []
for s in range(1, K + 1):
    m = st == s; n = int(m.sum())
    if n == 0: continue
    cc_ = g[m][g[m] > 0]; c = np.bincount(cc_) if len(cc_) else np.array([1]); maj = int(c.argmax()); rows.append((s, n, maj, c[maj] / max(len(cc_), 1), float((g[m] == 0).mean())))   # purity over conformer particles; junk fraction separately
rows.sort(key=lambda r: -r[1]); pur = np.array([r[3] for r in rows]); sizes = np.array([r[1] for r in rows]); cov = len(set(r[2] for r in rows if r[2] > 0))
print('%d merged states; particles per state median %d (min %d, max %d); mean purity %.2f (size-weighted %.2f); conformers covered %d; states with purity >= 0.5: %d' % (len(rows), np.median(sizes), sizes.min(), sizes.max(), pur.mean(), (pur * sizes).sum() / sizes.sum(), cov, int((pur >= 0.5).sum())))
print('  ' + '  '.join('%d:%d(conf%d %.2f, junk %.2f)' % r for r in rows[:12]))
