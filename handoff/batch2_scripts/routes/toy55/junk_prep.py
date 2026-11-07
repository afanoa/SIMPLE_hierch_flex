#!/usr/bin/env python3
"""usage: junk_prep.py smoke DIR | junk_prep.py simsites_truth SIMDIR CSV SITES_TRUTH OUT
smoke: tiny synthetic two-half stack (noise + a blob in genuine images, 20 % junk without it) to exercise junk_score.py end to end.
simsites_truth: sites_truth.txt plus a junk-type column (1 blob, 2 broken) from merge55.py's seeded permutation of the stacks (conformers, blob, broken)."""
import os, sys, numpy as np
if sys.argv[1] == 'smoke':
    D = sys.argv[2]; os.makedirs(D, exist_ok=True); os.chdir(D); rng = np.random.default_rng(1); n = 64; NM = 1200
    mr = rng.permutation(NM) + 1; jt = np.zeros(NM, int); jt[rng.choice(NM, 240, replace=False)] = rng.integers(1, 3, 240); lab = rng.integers(1, 3, NM)
    x = np.arange(n) - n // 2; blob = np.exp(-(x[:, None] ** 2 + x[None, :] ** 2) / 50.)
    for h, M in (('A', mr[:600]), ('B', mr[600:])):
        os.makedirs('sim_' + h, exist_ok=True); hdr = np.zeros(256, np.int32); hdr[0] = hdr[1] = n; hdr[2] = len(M)
        img = rng.normal(size=(len(M), n, n)).astype(np.float32) + 0.3 * blob[None] * (jt[M - 1] == 0)[:, None, None]
        with open('sim_%s/all.mrcs' % h, 'wb') as f: f.write(hdr.tobytes()); f.write(img.astype(np.float32).tobytes())
        np.savetxt('mr_%s.txt' % h, M, fmt='%d'); np.savetxt('lab_%s.txt' % h, lab[M - 1], fmt='%d')
        with open('poses_%s.txt' % h, 'w') as f:
            for i in range(len(M)): f.write('e1=%.1f e2=%.1f e3=%.1f x=0 y=0 dfx=%.2f dfy=%.2f angast=0\n' % (rng.uniform(0, 360), rng.uniform(0, 180), rng.uniform(0, 360), rng.uniform(1, 3), rng.uniform(1, 3)))
    with open('truth.txt', 'w') as f:
        f.write('# row junk jt C_closed th\n')
        for i in range(NM): f.write('%d %d %d %d %.2f\n' % (i + 1, jt[i] > 0, jt[i], rng.integers(0, 2), rng.normal()))
elif sys.argv[1] == 'simsites_truth':
    import csv
    S, CSV, TF, O = sys.argv[2:6]; names = ['conf%02d' % int(r['conformer']) for r in csv.DictReader(open(CSV))] + ['blob', 'broken']
    nz = [int(np.fromfile('%s/stk_%s.mrcs' % (S, k), np.int32, 3)[2]) for k in names]; src = np.repeat(np.arange(len(names)), nz)
    perm = np.random.default_rng(11).permutation(len(src)); G = np.loadtxt(S + '/gt_particles.txt', comments='#'); conf = G[:, 1].astype(int)
    labs = np.array([int(r['conformer']) for r in csv.DictReader(open(CSV))] + [0, 0]); assert (labs[src[perm]] == conf).all(), 'merge permutation not reproduced'
    off = np.r_[0, np.cumsum(nz)]
    for k in (len(names) - 2, len(names) - 1):   # each junk row's e1 equals the e1 at its place in its own junk stack
        e1 = np.array([float(dict(t.split('=', 1) for t in l.split() if '=' in t)['e1']) for l in open('%s/oris_%s.txt' % (S, names[k])) if l.strip()]); j = np.flatnonzero(src[perm] == k)
        assert np.allclose(e1[perm[j] - off[k]], G[j, 3], atol=1e-3), names[k]
    hd = [l.split()[1:] for l in open(TF) if l.startswith('#')][0]; T = np.loadtxt(TF, comments='#'); assert ((T[:, hd.index('junk')] > 0) == (conf == 0)).all()
    jt = np.where(src[perm] == len(names) - 2, 1, np.where(src[perm] == len(names) - 1, 2, 0))
    np.savetxt(O, np.c_[T, jt], fmt='%g', header=' '.join(hd + ['jtype']), comments='# ')
    print('simsites truth with junk type (merge permutation reproduced): %d blob, %d broken -> %s' % ((jt == 1).sum(), (jt == 2).sum(), O))
