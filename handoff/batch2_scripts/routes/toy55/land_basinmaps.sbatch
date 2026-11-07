#!/bin/bash
#SBATCH -J land_basinmaps
#SBATCH -p csbdevel
#SBATCH -c 24
#SBATCH --mem=80G
#SBATCH -t 03:00:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_basinmaps.%j.out
# basin maps of the masked pole-to-pole paths (runpath2_*): particles split at the common sparsest point s* (closed pole 0 -> open
# pole 1), SIMPLE maps per basin via land_rec.sbatch, then the coordinating session's landscape_cmp.py across runs. args: SSTAR
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; L=$K_/landscape; SS=${1:-0.43}
CMP=""
for r in script port1 port2; do D=$L/runpath2_$r
  $PY - <<PYEOF
import numpy as np
arm = [l.split() for l in open('$D/arm.txt') if not l.startswith('#')]; flip = float(arm[0][1]) < float(arm[1][1])
Z = np.load('$D/rw_lmarg_split_0.001/posterior.npz'); s = Z['smean'].astype(float); s = 1. - s if flip else s
lab = np.zeros(932496, int); lab[Z['rows']] = np.where(s < $SS, 1, 2); np.savetxt('$D/basins.txt', lab, fmt='%d')
print('$r: closed-pole basin %d particles (%.3f), open-pole basin %d (%.3f) at s* $SS' % ((lab == 1).sum(), (lab == 1).sum() / len(s), (lab == 2).sum(), (lab == 2).sum() / len(s)))
PYEOF
  OUT=$D/basinrec VIEW=$K_/view_rflex/runpath2_r1_${r}_basins LAB=$D/basins.txt PROJ=$K_/rflex_f692st/node_r1/fstageB/2_refine3D/b.simple NAME=basins_$r bash $T/land_rec.sbatch
  CMP="$CMP $r=$(ls -td $D/basinrec/*_reconstruct3D | head -1):$D/basins.txt"; done
$PY $T/landscape_cmp.py $CMP --lp 15
echo "[$(date '+%F %T')] land_basinmaps DONE"
