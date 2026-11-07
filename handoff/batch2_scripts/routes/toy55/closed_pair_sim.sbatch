#!/bin/bash
#SBATCH -J klem_closed_pair
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# region-local similarities of the two closed zoomed sub-states to the sub-states they were pooled with at cut 0.5
/mnt/nasapps/production/python/3.12/bin/python3 - <<PYEOF
import numpy as np
M = np.load('/mnt/beegfs/elmlund/afan/rnasep_klem/f692c/merge_tree_zoom_reg.npz', allow_pickle=True); S, names, Z = M['S'], list(M['names']), M['Z']
grp = 'f1_motion_1_s01_s01 f1_motion_1_s01_s02 f1_motion_1_s01_s03 f1_motion_1_s01_s04 f1_motion_1_s05_s01 f1_motion_1_s05_s03'.split(); idx = [names.index(g) for g in grp]
print('pairwise region-local similarity of group 7 (cut 0.5):'); print('%22s' % '' + ''.join('%8s' % g[-7:] for g in grp))
for i in idx: print('%22s' % names[i] + ''.join('%8.2f' % S[i, j] for j in idx))
from scipy.cluster.hierarchy import fcluster
for t in (0.5, 0.55, 0.6, 0.65, 0.7):
    lab = fcluster(Z, t=1 - t, criterion='distance'); c = lab[names.index('f1_motion_1_s05_s01')]; mem = [names[k] for k in np.where(lab == c)[0]]
    print('cut %.2f: %d states; the closed pair sits in a group of %d: %s' % (t, lab.max(), len(mem), ' '.join(mem)[:150]))
PYEOF
