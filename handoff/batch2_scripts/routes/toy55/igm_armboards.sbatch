#!/bin/bash
#SBATCH -J igm_armboard
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 03:00:00
#SBATCH --array=1-10
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%A_%a.out
# Factorised test on an IgM-like set: ONE focused board per arm (focus = the arm's sector mask, igm_arm_masks.py), random
# lattice start, all particles at the refined poses, 7x3 lattice + 1 free component; array task t -> arm ((t-1) % NARM) + 1,
# recipe std (oracle schedule) for t <= NARM, coarse (20 -> 16 A, sigmas narrowing to 0.3) after. Readout: armboard_readout.py.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; MODE=${MODE:-armsigm1}; NARM=${NARM:-5}
t=${SLURM_ARRAY_TASK_ID:?}; ARM=$(( (t - 1) % NARM + 1 )); if [ $t -le $NARM ]; then REC=std; LPS="20 14 12 10 10 8 8 8"; SIGS="1.5 1.2 1.0 1.0 0.8 0.8 0.7 0.7"
else REC=coarse; LPS="20 20 20 20 16 16 16 16"; SIGS="1.5 1.2 1.0 0.8 0.6 0.5 0.4 0.3"; fi
O=$W/armboards_$MODE/arm${ARM}_$REC; mkdir -p $O; cd $O; export FB=$A/rnasep_klem/f692e/klem2_bin SIMPLE_QSYS=local
PS=$A/rnasep_klem/synth_$MODE/src.simple; NP=$(grep -vc '^#' $W/sim$MODE/gt_particles.txt); MK=$W/maps160/synth_$MODE/masks/arm_$ARM.mrc
$PY -c "import numpy as np; np.savetxt('$O/subset.txt', np.arange(1, $NP + 1), fmt='%d')"
echo "[$(date '+%F %T') armboard $MODE] arm $ARM, recipe $REC (lps $LPS)"
[ -s $O/board/klem_latent.txt ] || SERVER_LIST= NPARTS= NTHR=16 N1=7 N2=3 NFREE=1 CTFABS=1 PSCALE=1 FOCUSMASK=$MK INIT=random SEED=$t bash $W/fboard.sh $O/board $PS $O/subset.txt "$LPS" "$SIGS" > $O/board.out 2>&1 || { echo "board FAILED: $(tail -2 $O/board.out)"; exit 1; }
$PY $W/armboard_readout.py $O/board $W/sim$MODE/gt_particles.txt $W/truth_$MODE.txt $ARM 2>&1 | grep -v Warning
echo "[$(date '+%F %T') armboard $MODE] DONE"
