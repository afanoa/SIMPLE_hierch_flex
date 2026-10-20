#!/bin/bash
#SBATCH -J klem_diag_geom
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# slice-geometry check of avgs_states3d.py on an existing direction-bin run (env MODE); also a real-space check:
# the consensus projected along the first bin directions by rotation + sum vs the bin averages (cavgs_bins.mrc)
set -uo pipefail; PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; MODE=${MODE:-tethered}; D=$W/binsroute_$MODE
DIAG_ONLY=1 $PY $W/avgs_states3d.py $D/stitch $D/p $W/maps160/synth_$MODE/consensus_synth.mrc $D/diag --k 6 2>&1 | grep -v Warning
$PY - <<PYEOF
import numpy as np, sys
sys.path.insert(0, '$W'); from model_class import read_mrc
from scipy.ndimage import affine_transform
V = read_mrc('$W/maps160/synth_$MODE/consensus_synth.mrc').astype(np.float32); n = V.shape[0]
h = np.frombuffer(open('$D/p/cavgs_bins.mrc', 'rb').read(1024), np.int32); nsym = int(h[23]); A = np.array(np.memmap('$D/p/cavgs_bins.mrc', np.float32, 'r', offset=1024 + nsym, shape=(int(h[2]), int(h[1]), int(h[0]))))
B = np.load('$D/bins.npz'); C = B['centres']   # unit vectors (x, y, z) of the bin directions, e1/e2 convention of dir_bins
def rot(e1d, e2d):
    a1, a2 = np.radians([e1d, e2d]); c1, s1, c2, s2 = np.cos(a1), np.sin(a1), np.cos(a2), np.sin(a2)
    return np.array([[c1, -s1, 0], [s1, c1, 0], [0, 0, 1]]) @ np.array([[c2, 0, s2], [0, 1, 0], [-s2, 0, c2]])
E2 = np.degrees(np.arccos(np.clip(C[:, 2], -1, 1))); E1 = np.degrees(np.arctan2(C[:, 1], C[:, 0])) % 360.
def project(Rm, T=False):   # rotate the volume so that the projection direction becomes z (array axis 0), then sum
    M = Rm.T if T else Rm; Mz = M[::-1, ::-1]   # (x,y,z) matrix -> array (z,y,x) order
    ctr = np.array([n // 2] * 3); W_ = affine_transform(V, Mz, offset=ctr - Mz @ ctr, order=1); return W_.sum(0)
def cc(a, b): a = a - a.mean(); b = b - b.mean(); return float((a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum() + 1e-30))
x = np.arange(n) - n // 2; m2 = np.hypot(x[:, None], x[None, :]) < 0.45 * n
for T in (False, True):
    for yf in (False, True):
        for xf in (False, True):
            r = []
            for k in range(min(12, len(C))):
                P = project(rot(E1[k], E2[k]), T); P = P[::-1] if yf else P; P = P[:, ::-1] if xf else P; r.append(cc(P[m2], A[k][m2]))
            print('real-space: transpose=%s yflip=%s xflip=%s: mean corr(projection, bin average) %.3f (per bin: %s)' % (T, yf, xf, np.mean(r), ' '.join('%.2f' % v for v in r[:6])))
PYEOF
echo "DIAG GEOM DONE"
