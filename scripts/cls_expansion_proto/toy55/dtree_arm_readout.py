#!/usr/bin/env python3
"""Per-arm truth readout of a data-derived tree (dtree_v10) on a hinged-arm synthetic set: for every node of both halves,
each child group's std of every arm angle (theta, phi) relative to the std over the node's particles and over all
particles; a child that localises arm k has a small ratio for arm k. Ends with the best ratio per arm over all children.
usage: dtree_arm_readout.py K_DIR MODE TAG SIM_BASE TRUTH   (half run dirs K_DIR/synth_<MODE>_A|B/TAG/<node>, truth rows
SIM_BASE_A|B/gt_particles.txt)"""
import sys, os, glob, numpy as np
K, mode, tag, simb, trf = sys.argv[1:6]
ang = {int(w[0]): np.array([float(x) for x in w[1:]]) for w in (l.split() for l in open(trf)) if w and not w[0].startswith('#')}
narm = len(next(iter(ang.values()))) // 2; names = ['a%d%s' % (k // 2 + 1, 'th' if k % 2 == 0 else 'ph') for k in range(2 * narm)]
best = np.ones(2 * narm)
for H in 'AB':
    G = np.loadtxt('%s_%s/gt_particles.txt' % (simb, H), comments='#'); conf = G[:, 1].astype(int)
    Aall = np.array([ang.get(c, np.full(2 * narm, np.nan)) for c in conf]); sall = np.nanstd(Aall, 0)
    print('=== half %s: %d rows, %d arms; angle std over all conformer particles %s' % (H, len(conf), narm, ' '.join('%.1f' % v for v in sall)))
    for nd in sorted(glob.glob('%s/synth_%s_%s/%s/r*' % (K, mode, H, tag)), key=lambda p: (p.count('.'), p)):
        f = nd + '/children.txt'
        if not os.path.exists(f): continue
        lab = np.loadtxt(f, dtype=int, ndmin=1); sub = np.loadtxt(nd + '/subset.txt', dtype=int, ndmin=1) - 1
        snode = np.nanstd(Aall[sub], 0); ks = sorted(set(lab[lab > 0].tolist()))
        print('-- node %s: %d particles, %d children (std over the node %s)' % (os.path.basename(nd), len(sub), len(ks), ' '.join('%.1f' % v for v in snode)))
        for c in ks:
            rows = np.where(lab == c)[0]; a = Aall[rows]; a = a[~np.isnan(a[:, 0])]
            if len(a) < 50: print('   child %d: %d particles (%d junk)' % (c, len(rows), len(rows) - len(a))); continue
            r = np.nanstd(a, 0) / sall; k = int(np.argmin(r)); best = np.minimum(best, r)
            print('   child %d: %5d ptcls (junk %4.0f %%)  ratio vs all %s  -> %s %.2f (mean %+.0f deg)' % (c, len(rows), 100. * (len(rows) - len(a)) / len(rows), ' '.join('%.2f' % v for v in r), names[k], r[k], np.mean(a[:, k])))
print('== best std ratio per arm angle over all children: %s' % ' '.join('%s %.2f' % (n, v) for n, v in zip(names, best)))
