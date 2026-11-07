#!/bin/bash
#SBATCH -J klem_regmerge
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# do a finished run's region boards track the same motions? region_merge.py (canonical correlation of the boards' coordinates,
# one half fitted, the other measured) on each family's region boards; the groups it forms are the boards that would be one
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; R=/mnt/beegfs/elmlund/afan/rnasep_klem/${RTAG:?}; export OMP_NUM_THREADS=8
for F in 1 2; do D=$R/fam$F; C=""
  for RB in $D/region*/; do b=$(basename $RB); reg=${b#region}; case $reg in motion_*) mk=$D/regions_motion/region_${reg#motion_}.mrc;; *) mk=$D/regions/region$reg.mrc;; esac
    [ -s $RB/board/klem_latent.txt ] && [ -s $mk ] && C="$C $mk:${RB%/}/board"; done
  echo "== family $F: $(echo $C | wc -w) region boards"; $PY $W/region_merge.py $R/region_merge_check/fam$F ${THR:-gap:1.5} $C 2>&1 | grep -v Warning
  cat $R/region_merge_check/fam$F/groups.txt 2>/dev/null
done
