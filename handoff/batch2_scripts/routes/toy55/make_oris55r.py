#!/usr/bin/env python3
"""Per-conformer simulate_particles oritabs for the RNase P 55-conformer toy: poses (e1 e2 e3) and CTF (dfx dfy angast)
drawn without replacement from the refined real particles (meta_all.npz, flags==1), conformer weights = the paper's
image counts (CSV), plus two junk classes (blob, broken) at JUNK_FRAC each. usage: make_oris55.py META CSV OUT N JUNK_FRAC SEED"""
import sys, csv, numpy as np
meta, csvf, out, N, JF, seed = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6])
rng = np.random.default_rng(seed); M = np.load(meta); ok = np.where(M['flags'] == 1)[0]
pick = rng.choice(ok, N, replace=N > len(ok)); e1, e2, e3, dfx, dfy, ang = (M[k][pick] for k in ('e1', 'e2', 'e3', 'dfx', 'dfy', 'angast'))   # with replacement beyond the real pool (scaling sets)
rows = list(csv.DictReader(open(csvf))); w = np.array([float(r['images']) for r in rows]); w /= w.sum()
nj = int(round(JF * N)); nc = N - 2 * nj; counts = np.floor(w * nc).astype(int); counts[np.argmax(w)] += nc - counts.sum()
lab = np.concatenate([np.full(c, int(r['conformer'])) for c, r in zip(counts, rows)] + [np.full(nj, -1), np.full(nj, -2)])   # -1 blob, -2 broken
lab = lab[rng.permutation(N)]
def write(name, idx):
    with open('%s/oris_in_%s.txt' % (out, name), 'w') as f:
        for i in idx: f.write('e1=%.4f e2=%.4f e3=%.4f x=0.000 y=0.000 dfx=%.5f dfy=%.5f angast=%.4f\n' % (e1[i], e2[i], e3[i], dfx[i], dfy[i], ang[i]))
for r in rows:
    k = int(r['conformer']); write('conf%02d' % k, np.where(lab == k)[0])
write('blob', np.where(lab == -1)[0]); write('broken', np.where(lab == -2)[0])
print('%d particles: %d conformer (weights p10/50/90 %.4f/%.4f/%.4f, counts min %d max %d), %d blob, %d broken; defocus mean %.2f sd %.2f' %
      (N, nc, *np.percentile(w, (10, 50, 90)), counts.min(), counts.max(), nj, nj, dfx.mean(), dfx.std()))
