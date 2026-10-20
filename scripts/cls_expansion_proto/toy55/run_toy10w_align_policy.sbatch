#!/bin/bash
#SBATCH -J toy10w_align_policy
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 10:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Node alignment policy with the truth known (wide toy, closed arm intact): waits for the toy10w simulation; one-cell
# klem boards of the TRUE particles of conformers 18 (closed), 51 (closed) and 27 (open) at the refined consensus poses
# (angle error ~4 deg), 10 -> 8 A, with klem's strictly local pose step (3x3x3 rotation grid x 2 refinements, Gaussian
# prior, no shifts) against (a) nothing (baseline), (b) the cell reference, everything (lever arm), (c) the cell
# masked to the CORE (density minus the dilated P12 arm), (d) masked to the ARM REGION (deposited 18 minus the
# consensus, dilated). Readouts: pose error vs truth after each pose step (--pose_truth), closed-arm occupancy of the
# conformer-18 cell map (closed_arm_cells.py), pose change. env: STEP (3), PRIOR (4), EVERY (2)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; S=$W/simw; PJ=$W/projw; L=$W/latw; mkdir -p $L; log(){ echo "[$(date '+%F %T') ap] $*"; }; export OMP_NUM_THREADS=1
n=0; while [ $n -lt 300 ]; do f=$(ls -t $W/toy10w_sim.*.out 2>/dev/null | head -1); [ -n "$f" ] && grep -q "DONE toy55" $f && break; sleep 60; n=$((n+1)); done
f=$(ls -t $W/toy10w_sim.*.out | head -1); grep -q "DONE toy55" $f || { log "toy10w simulation not finished"; exit 1; }
grep -a "calibration:\|refined vs true poses" $f | cut -c1-200
cp -f $PJ/ptcl3D_refined.txt $L/ptcl3D_r0.txt; $PY $W/pose_error.py $L/ptcl3D_r0.txt $S/all.txt --tag "toy10w refined consensus poses" 2>&1 | grep -v Warning
RP=$(ls -td $PJ/*_refine3D | head -1); CONS=$RP/recvol_state01.mrc
STATES=$(awk -F, 'NR>1 {printf "%s%s", (NR>2?",":""), $2}' $W/rnasep_conformers10.csv); AXIS=$(echo $STATES | sed 's/[0-9][0-9]*/0/g')
for c in 18 51 27; do awk -v c=$c '!/^#/ && $2 == c {print $1}' $S/gt_particles.txt > $L/conf${c}_subset.txt; log "conformer $c: $(grep -c . $L/conf${c}_subset.txt) true particles"; done
$PY - <<PYEOF
import numpy as np, glob, sys
sys.path.insert(0, '$W'); from model_class import read_mrc, write_mrc
from scipy.ndimage import gaussian_filter, binary_dilation
cons = read_mrc('$W/maps160/toy_consensus.mrc'); c18 = read_mrc(glob.glob('$W/maps160/docked/conf18_*.mrc')[0]); g = lambda v: gaussian_filter(v, 2.5); L18, Lc = g(c18), g(cons)
dens18 = L18 > L18.mean() + 2.0 * L18.std(); diff = L18 - Lc; arm = (diff > diff.mean() + 3 * diff.std()) & dens18
reg = binary_dilation(arm, iterations=3); sr = gaussian_filter(reg.astype(float), 2.0); sr = np.clip(sr / max(sr.max(), 1e-30), 0, 1); write_mrc('$L/arm_region_mask.mrc', sr.astype(np.float32), 1.464)
core = dens18 & ~binary_dilation(arm, iterations=4); sc = gaussian_filter(core.astype(float), 2.0); sc = np.clip(sc / max(sc.max(), 1e-30), 0, 1); write_mrc('$L/core_mask.mrc', sc.astype(np.float32), 1.464)
print('masks: arm %d voxels -> region %d (dilated), core %d (density %d minus the dilated arm)' % (arm.sum(), (sr > 0.5).sum(), (sc > 0.5).sum(), dens18.sum()))
PYEOF
run(){ c=$1; tag=$2; shift 2; D=$L/policy_c${c}_$tag; log "conformer $c, $tag: $*"
  $PY $W/klem.py $S/all.mrcs $L/ptcl3D_r0.txt $S/gt_particles.txt $D --cons $CONS --states $STATES --axis=$AXIS --n1 1 --n2 1 --nfree 0 --nproc 16 --init random --nsub 0 --subset $L/conf${c}_subset.txt \
      --lps "10 10 10 8 8 8" --sigmas "0.7 0.7 0.7 0.7 0.7 0.7" --mskdiam 232 --pose_truth $S/all.txt --tag "c$c $tag" "$@" 2>&1 | grep "subset\|pose step\|Error\|Traceback" | cut -c1-230
  [ $c = 18 ] && $PY $W/closed_arm_cells.py $D --top 1 2>&1 | grep -v Warning | sed "s/^/[c$c $tag] /"; return 0; }
for c in 18 51 27; do
  run $c baseline
  run $c lever  --pose_every ${EVERY:-2} --pose_step ${STEP:-3} --pose_prior ${PRIOR:-4} --pose_ref cell
  run $c core   --pose_every ${EVERY:-2} --pose_step ${STEP:-3} --pose_prior ${PRIOR:-4} --pose_ref cell --pose_mask $L/core_mask.mrc
  run $c region --pose_every ${EVERY:-2} --pose_step ${STEP:-3} --pose_prior ${PRIOR:-4} --pose_ref cell --pose_mask $L/arm_region_mask.mrc
done
log DONE
