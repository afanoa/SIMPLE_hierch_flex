#!/bin/bash
#SBATCH -J klem_soft_maps
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# soft state maps + per-particle weights for every leaf (env RTAG, LEAFGLOB region*|zoomregion*)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/${RTAG:?}
for d in $R/fam?/${LEAFGLOB:-region*}/; do [ -s $d/level1/clusters.txt ] || continue; echo "== $(echo $d | sed "s#$R/##")"; $PY $W/soft_state_maps.py $d/board $d/level1/clusters.txt $W/maps160/toy_consensus.mrc $d/level1/soft 2>&1 | grep -v Warning | head -3; done
echo "SOFT MAPS DONE"
