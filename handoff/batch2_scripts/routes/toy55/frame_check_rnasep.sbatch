#!/bin/bash
#SBATCH -J rnasep_frame_check
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Is the re-imported reconstruction the same molecule in another frame? Correlation of the re-imported map (unmoved
# poses) with the family-pose map under identity, axis flips (mirror hands), 180-degree rotations, at 10 A.
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; KD=$A/rnasep_klem/expand1/fam1_pipeline/klem_fb
$PY - <<PYEOF
import numpy as np, sys
sys.path.insert(0, '$W'); from model_class import read_mrc
from scipy.ndimage import gaussian_filter
a = read_mrc('$KD/closed_focal_docked30k/rank/maps/recvol_state01.mrc'); b = read_mrc('$KD/region_align_10k/import_check/plain/maps/recvol_state01.mrc'); n = a.shape[0]
g = lambda v: gaussian_filter(v, 3.0); A_ = g(a); B_ = g(b); A_ = (A_ - A_.mean()) / A_.std()
def cc(x): x = (x - x.mean()) / x.std(); return float(np.mean(A_ * x))
cands = {'identity': B_, 'flip z (mirror)': B_[::-1], 'flip y': B_[:, ::-1], 'flip x': B_[:, :, ::-1], 'rot180 about z': B_[::-1, ::-1], 'rot180 about y': B_[::-1, :, ::-1], 'rot180 about x': B_[:, ::-1, ::-1], 'invert all': B_[::-1, ::-1, ::-1]}
for k, v in cands.items(): print('  %-18s corr %.3f' % (k, cc(v)))
print('  (maps smoothed at ~3 px; a frame/hand relation shows as a corr near 1 for one candidate)')
PYEOF
