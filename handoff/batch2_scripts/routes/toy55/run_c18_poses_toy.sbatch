#!/bin/bash
#SBATCH -J toy10_c18_poses
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Pose blur vs the closed arm, on the toy where the truth is known: the TRUE conformer-18 particles reconstructed (a
# one-cell klem board = plain cross-half reconstruction at 8 A, wide mask) at the TRUE poses, at the realigned poses and
# at the consensus poses; closed-arm occupancy of each on the coarse grid. If the realigned-pose map already loses the
# arm, pose blur is the cap; if not, the real particles' arm is genuinely partial.
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; cd $W; log(){ echo "[$(date '+%F %T') c18] $*"; }
awk '!/^#/ && $2 == 18 {print $1}' $W/sim/gt_particles.txt > $W/lat/conf18_subset.txt; log "true conformer-18 particles: $(grep -c . $W/lat/conf18_subset.txt)"
for v in "true:$W/sim/all.txt" "realigned:$W/lat/realign_realign1/ptcl3D_realigned.txt" "consensus:$W/lat/ptcl3D_r0.txt"; do
  name=${v%%:*}; tab=${v#*:}; log "one-cell board of the true 18 particles at the $name poses"
  env KTAG=_c18_$name INIT=random NSUB=0 N1=1 N2=1 NFREE=0 SUBSET=$W/lat/conf18_subset.txt TABLE=$tab MSKD=232 LPS="8 8" SIGMAS="0.7 0.7" bash $W/run_klem2.sbatch 2>&1 | grep "subset\|Error\|Traceback" | sed "s/^/[$name] /"
  $PY $W/closed_arm_cells.py $W/lat/klem_c18_$name --top 1 2>&1 | grep -v Warning | sed "s/^/[$name] /"
done
log DONE
