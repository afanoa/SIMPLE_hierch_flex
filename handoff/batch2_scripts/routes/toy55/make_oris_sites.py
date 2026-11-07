#!/usr/bin/env python3
"""simsites oritabs: every non-junk particle draws its three site states INDEPENDENTLY (Bernoulli PA / PB / PC from COMBOS = combos.txt
of make_sites_maps.py), so the combination counts are multinomial; poses (e1 e2 e3) and CTF (dfx dfy angast) drawn without replacement
from the refined real particles (meta_all.npz, flags == 1), as make_oris55.py; two junk classes (blob, broken) at JUNK_FRAC each.
Writes OUT/oris_in_conf0k.txt (k = combination id) + oris_in_blob.txt / oris_in_broken.txt.
usage: make_oris_sites.py META COMBOS OUT N JUNK_FRAC SEED"""
import sys, numpy as np
meta, combos, out, N, JF, seed = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6])
rng = np.random.default_rng(seed); M = np.load(meta); ok = np.where(M['flags'] == 1)[0]
pick = rng.choice(ok, N, replace=False); e1, e2, e3, dfx, dfy, ang = (M[k][pick] for k in ('e1', 'e2', 'e3', 'dfx', 'dfy', 'angast'))
C = np.loadtxt(combos, comments='#'); pa = C[C[:, 1] == 1, 4].sum(); pb = C[C[:, 2] == 1, 4].sum(); pc = C[C[:, 3] == 1, 4].sum()   # marginals
nj = int(round(JF * N)); nc = N - 2 * nj
A = rng.random(nc) < pa; B = rng.random(nc) < pb; Cc = rng.random(nc) < pc; k = 1 + (~A).astype(int) + 2 * B + 4 * Cc
lab = np.concatenate([k, np.full(nj, -1), np.full(nj, -2)])[rng.permutation(N)]   # -1 blob, -2 broken
def write(name, idx):
    with open('%s/oris_in_%s.txt' % (out, name), 'w') as f:
        for i in idx: f.write('e1=%.4f e2=%.4f e3=%.4f x=0.000 y=0.000 dfx=%.5f dfy=%.5f angast=%.4f\n' % (e1[i], e2[i], e3[i], dfx[i], dfy[i], ang[i]))
for kk in range(1, 9): write('conf%02d' % kk, np.where(lab == kk)[0])
write('blob', np.where(lab == -1)[0]); write('broken', np.where(lab == -2)[0])
print('%d particles: %d non-junk (domain present %.3f, blob %.3f, arm closed %.3f; marginals %.2f / %.2f / %.2f), %d blob + %d broken junk; combination counts %s' %
      (N, nc, A.mean(), B.mean(), Cc.mean(), pa, pb, pc, nj, nj, ' '.join('%d:%d' % (kk, (lab == kk).sum()) for kk in range(1, 9))))
