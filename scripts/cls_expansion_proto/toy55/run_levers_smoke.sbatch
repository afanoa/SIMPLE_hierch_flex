#!/bin/bash
#SBATCH -J toy10_levers_smoke
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Smoke test of anneal_split.py and density_modes.py on a synthetic klem directory (3 cell families, 6 particle states).
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; S=$W/lat/synth_klem; mkdir -p $S
$PY - <<PYEOF
import numpy as np
rng = np.random.default_rng(1); S = '$S'; ncell = 49; i1, rem = np.divmod(np.arange(ncell), 7); i2 = rem
V = 3000; sig = {f: rng.normal(size=V) + 1j * rng.normal(size=V) for f in range(3)}; fam = np.minimum(i1 // 3, 2)
D0 = np.array([sig[fam[k]] * (1 + 0.2 * i2[k]) + 0.8 * (rng.normal(size=V) + 1j * rng.normal(size=V)) for k in range(ncell)]); D1 = np.array([sig[fam[k]] * (1 + 0.2 * i2[k]) + 0.8 * (rng.normal(size=V) + 1j * rng.normal(size=V)) for k in range(ncell)])
n = 30000; centres = np.array([[1, 1], [1, 5], [3, 3], [5, 1], [5, 5], [3, 0.5]], float); truth = rng.integers(0, 6, n); pos = centres[truth] + rng.normal(scale=0.5, size=(n, 2))
C = np.c_[i1, i2].astype(float); d2 = ((pos[:, None, :] - C[None]) ** 2).sum(2); r = np.exp(-0.5 * d2 / 0.8 ** 2); r = np.c_[r, 0.01 * np.ones((n, 2))]; r /= r.sum(1, keepdims=True)
w = r[:, :ncell].sum(0); z = (r[:, :ncell] / r[:, :ncell].sum(1, keepdims=True)) @ C; am = r.argmax(1)
np.savez_compressed(S + '/klem_halves.npz', D0=D0.astype(np.complex64), D1=D1.astype(np.complex64), w=w, i1=i1, i2=i2, n3=31, band=np.arange(V))
np.savez_compressed(S + '/klem_resp.npz', pinds=np.arange(1, n + 1), r=r.astype(np.float32), i1=i1, i2=i2)
with open(S + '/klem_latent.txt', 'w') as f:
    f.write('# pind z1 z2 argmax maxr gridmass layer2mass\n')
    for j in range(n): f.write('%d %.4f %.4f %d %.4f %.4f 0\n' % (j + 1, z[j, 0], z[j, 1], am[j] + 1, r[j].max(), r[j, :ncell].sum()))
with open(S + '/gt_particles.txt', 'w') as f:
    f.write('# pind conformer closed\n')
    for j in range(n): f.write('%d %d 0\n' % (j + 1, truth[j] + 1))
print('synthetic klem dir: 3 cell families (columns 0-2, 3-5, 6), 6 particle states at (1,1) (1,5) (3,3) (5,1) (5,5) (3,0.5)')
PYEOF
echo "== anneal_split (expect 3 families)"; $PY $W/anneal_split.py $S $S/anneal --gt $S/gt_particles.txt
echo "== density_modes (expect 6 modes)"; $PY $W/density_modes.py $S --gt $S/gt_particles.txt --k0 12 --kmeans 6 --nfit 20000 --iters 60
echo SMOKE DONE
