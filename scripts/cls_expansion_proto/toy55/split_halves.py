#!/usr/bin/env python3
"""Split a merged simulation (all.mrcs/all.txt/imp.txt/gt_particles.txt) into two independent halves (alternate rows)
as sim<MODE>_A / sim<MODE>_B, with marker files so run_toy55_scale.sbatch skips simulation and merging.
usage: split_halves.py SIMDIR OUT_A OUT_B"""
import sys, os, numpy as np
S, OA, OB = sys.argv[1:4]
h = np.fromfile(S + '/all.mrcs', dtype=np.int32, count=256); nx, ny, N, nsym = int(h[0]), int(h[1]), int(h[2]), int(h[23])
mm = np.memmap(S + '/all.mrcs', dtype=np.float32, mode='r', offset=1024 + nsym, shape=(N, ny, nx))
for O, par in ((OA, 0), (OB, 1)):
    os.makedirs(O, exist_ok=True); idx = np.arange(par, N, 2); hdr = h.copy(); hdr[2] = len(idx); hdr[9] = len(idx); hdr[23] = 0
    with open(O + '/all.mrcs', 'wb') as f:
        f.write(hdr.tobytes())
        for c0 in range(0, len(idx), 20000): f.write(np.ascontiguousarray(mm[idx[c0:c0 + 20000]]).tobytes())
    for fn in ('all.txt', 'imp.txt'):
        L = [l for l in open('%s/%s' % (S, fn)) if l.strip()]; assert len(L) == N, (fn, len(L), N); open('%s/%s' % (O, fn), 'w').writelines(L[par::2])
    G = [l for l in open(S + '/gt_particles.txt') if l.strip() and not l.startswith('#')]; assert len(G) == N
    with open(O + '/gt_particles.txt', 'w') as g:
        g.write('# particle conformer closed e1 e2 (truth; conformer 0 = junk)\n')
        for j, i in enumerate(idx): g.write('%d %s' % (j + 1, G[i].split(' ', 1)[1]))
    open(O + '/oris_in_blob.txt', 'w').write('# split half: no simulation\n# (marker)\n'); open(O + '/stk_blob.mrcs', 'wb').write(b'\0')
    print('%s: %d particles (rows %d::2 of %s)' % (O, len(idx), par, S))
