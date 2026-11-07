#!/bin/bash
#SBATCH -J rnasep_sharpen_matched2
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; P=/mnt/beegfs/elmlund/afan/rnasep_klem/expand1/fam1_pipeline/klem_fb/matched_rec
$PY $W/sharpen_matched.py $P/maps $P/for_download_r116 --rmask 116
echo '== deposited conformers: is the arm ordered in the open ones? (density mass above threshold, relative to 18)'
$PY - <<PYEOF
import numpy as np, glob
W = '$W'
def read_mrc(p):
    h = np.fromfile(p, dtype=np.int32, count=256); nx, ny, nz = h[0], h[1], h[2]; return np.fromfile(p, dtype=np.float32, offset=1024 + h[23]).reshape(nz, ny, nx).astype(float)
ids = [1, 10, 18, 27, 28, 37, 46, 48, 49, 51]; C = {i: read_mrc(glob.glob(W + '/maps160/docked/conf%02d_*.mrc' % i)[0]) for i in ids}; thr = {i: C[i].mean() + 2.0 * C[i].std() for i in ids}
m18 = C[18][C[18] > thr[18]].sum(); n = 160; zz, yy, xx = np.indices((n, n, n)); r = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2) * 1.464
for i in ids:
    d = C[i] > thr[i]; print('  conformer %2d: mass above 2 sd %.2f of 18 | voxels %d | max radius %.0f A | mass beyond 90 A %.2f of 18' % (i, C[i][d].sum() / m18, d.sum(), r[d].max(), C[i][d & (r > 90)].sum() / max(C[18][(C[18] > thr[18]) & (r > 90)].sum(), 1e-30)))
PYEOF
