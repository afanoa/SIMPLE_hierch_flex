#!/bin/bash
#SBATCH -J toy10_c18_poses6
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Tip ceiling at a band matched to the real maps (6 A): pure true-18 particles, one-cell boards at the true and the realigned poses.
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; cd $W
for v in "true:$W/sim/all.txt" "realigned:$W/lat/realign_realign1/ptcl3D_realigned.txt"; do
  name=${v%%:*}; tab=${v#*:}
  env KTAG=_c18_${name}_lp6 INIT=random NSUB=0 N1=1 N2=1 NFREE=0 SUBSET=$W/lat/conf18_subset.txt TABLE=$tab MSKD=232 LPS="6 6" SIGMAS="0.7 0.7" bash $W/run_klem2.sbatch 2>&1 | grep "Error\|Traceback" | sed "s/^/[$name] /"
  echo "== pure true-18 at the $name poses, 6 A band"; $PY $W/closed_arm_cells.py $W/lat/klem_c18_${name}_lp6 --top 1 2>&1 | grep -v Warning | tail -1
done
echo DONE
