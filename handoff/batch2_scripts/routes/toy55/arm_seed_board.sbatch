#!/bin/bash
#SBATCH -J klem_arm_seed
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Focused two-cell EM for the closed arm on family 1 of f692c (335k, refined poses): cell 1 seeded with the particles of the
# leaf state with the highest closed-arm occupancy (SEED list), cell 2 with the rest; focus = the arm's travel mask, no
# balance, no kernel smoothing (1 x 2 lattice); readout = closed-arm occupancy of both cells after the E/M iterations.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/f692c; O=$A/rnasep_klem/arm_seed_f692c; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') armseed] $*"; }
CL=$(ls $W/maps160/docked/conf18_*_m160.mrc | head -1); [ -s $O/arm_mask.mrc ] || $PY $W/make_arm_mask.py $W/maps160/toy_consensus.mrc $CL $O/arm_mask.mrc 2>&1 | grep -v Warning
SEED=${SEED:-$R/fam1/regionmotion_1/level2_s1/states_top.txt}; SUB=$R/expand/family1_subset.txt; PS=$(ls -td $R/fam1/refine/*_refine3D_states | head -1)/p.simple
$PY - <<PYEOF
import numpy as np
seed = np.loadtxt('$SEED', dtype=int); rows = np.where(seed > 0)[0] + 1; sub = np.loadtxt('$SUB', dtype=int).ravel(); s = set(rows.tolist())
with open('$O/init_positions.txt', 'w') as f:
    for p in sub: f.write('%d 0 %d\n' % (p, 0 if p in s else 1))
print('seed: %d closed-ish particles of %d in family 1 -> cell 1, the rest cell 2' % (len(s & set(sub.tolist())), len(sub)))
PYEOF
log "two-cell arm-focused board (10 workers)"; N1=1 N2=2 NFREE=0 BAL=0 NTHR=8 NPARTS=10 CTFABS=1 PSCALE=1 SERVER_LIST= INIT=file INITFILE=$O/init_positions.txt FOCUSMASK=$O/arm_mask.mrc META=$R/fam1/meta_fb.npz bash $W/fboard.sh $O/board $PS $SUB "12 12 10 10 8 8 8 8" "0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3" 2>&1 | grep -v Warning | tail -4
echo "== closed-arm occupancy per cell"; $PY $W/closed_arm_cells.py $O/board --top 2 2>&1 | grep -v Warning | tail -4 | cut -c1-200
$PY - <<PYEOF
import numpy as np
L = np.loadtxt('$O/board/klem_latent.txt', comments='#'); am = L[:, 3].astype(int); print('final assignment: cell 1 %d particles, cell 2 %d' % ((am == 1).sum(), (am == 2).sum()))
PYEOF
log DONE
