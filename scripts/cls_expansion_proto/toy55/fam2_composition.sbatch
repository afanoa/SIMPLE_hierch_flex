#!/bin/bash
#SBATCH -J rnasep_fam2_composition
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Is real family 2 (197k, matches nothing deposited) a compositional state? Mean cell map of each family's board vs the
# consensus: ordered mass above threshold, the largest connected region of missing density and of extra density.
PY=/mnt/nasapps/production/python/3.12/bin/python3
$PY - <<'PYEOF'
import numpy as np, glob
from scipy.ndimage import gaussian_filter, label
W = '/mnt/beegfs/elmlund/afan/routes/toy55'; E = '/mnt/beegfs/elmlund/afan/rnasep_klem/expand1'
def read_mrc(p):
    h = np.fromfile(p, dtype=np.int32, count=256); nx, ny, nz = h[0], h[1], h[2]; return np.fromfile(p, dtype=np.float32, offset=1024 + h[23]).reshape(nz, ny, nx).astype(float)
fs = {j: sorted(glob.glob(E + '/fam%d/cells/cell_*.mrc' % j)) for j in (1, 2)}; print('cells per family:', {j: len(v) for j, v in fs.items()})
fams = {j: np.mean([read_mrc(f) for f in fs[j]], 0) for j in (1, 2)}; nc = fams[1].shape[0]; smpd_c = 160 * 1.464 / nc
cons = read_mrc(W + '/maps160/toy_consensus.mrc'); n = 160
def crop(v):
    Fv = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v))); c = n // 2; h = nc // 2; sub = Fv[c - h:c - h + nc, c - h:c - h + nc, c - h:c - h + nc]; return np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(sub)))) * (nc / n) ** 3
cc = gaussian_filter(crop(cons), 1.0); g = lambda v: gaussian_filter(v, 1.0)
zz, yy, xx = np.indices((nc, nc, nc)); r = np.hypot(np.hypot(xx - nc / 2, yy - nc / 2), zz - nc / 2) * smpd_c
dens = cc > cc.mean() + 1.5 * cc.std(); thr = cc[dens].mean() * 0.5; mass_c = (cc > thr).sum(); sc = {}
for j in (1, 2):
    v = g(fams[j]); s = np.sum(v[dens] * cc[dens]) / np.sum(v[dens] ** 2); v = v * s; sc[j] = v
    d = v - cc; sd = d[dens].std(); missing = (d < -2.5 * sd) & (cc > thr); extra = (d > 2.5 * sd) & (v > thr)
    for name, m in (('missing', missing), ('extra', extra)):
        lab, nl = label(m); sizes = np.bincount(lab.ravel())[1:] if nl else np.array([0]); k = int(np.argmax(sizes)) + 1 if nl else 0
        print('  family %d: ordered voxels above half-core %.2f of the consensus | %s density: %d voxels in %d regions, largest %d (%.0f A^3) at median radius %.0f A' % (j, (v > thr).sum() / max(mass_c, 1), name, m.sum(), nl, sizes.max() if nl else 0, (sizes.max() if nl else 0) * smpd_c ** 3, np.median(r[lab == k]) if nl else 0))
d12 = sc[1] - sc[2]; sd = d12[dens].std(); big = np.abs(d12) > 3 * sd; lab, nl = label(big); sizes = np.bincount(lab.ravel())[1:] if nl else np.array([0]); k = int(np.argmax(sizes)) + 1 if nl else 0
print('family 1 - family 2: %d voxels beyond 3 sd in %d regions; largest %d voxels (%.0f A^3) at median radius %.0f A, where family %s has more density' % (big.sum(), nl, sizes.max() if nl else 0, (sizes.max() if nl else 0) * smpd_c ** 3, np.median(r[lab == k]) if nl else 0, '1' if nl and d12[lab == k].mean() > 0 else '2'))
PYEOF
