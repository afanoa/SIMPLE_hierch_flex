#!/bin/bash
#SBATCH -J klem_recover55
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# per-conformer recovery (all 55 apo maps) of the leaf states and of every consolidated set of a run (env RTAG)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/${RTAG:?}
echo "== leaf states"; [ -s $R/readout55/matches.txt ] && $PY $W/recover55.py $R/readout55/matches.txt --names $R/readout55/map_order.txt 2>&1 | grep -v Warning | head -3
for C in $R/consolidated_*; do [ -d $C/fam0/regionall/level1/states_dir ] || continue; echo "== $(basename $C)"; S=$C/readout55; mkdir -p $S; ls $C/fam0/regionall/level1/states_dir/state_*.mrc | xargs -n1 basename > $S/map_order.txt
  $PY $W/deposited_match.py $C/fam0/regionall/level1/states_dir $W/maps160/docked_all $W/maps160/toy_consensus.mrc --save $S/matches.txt > /dev/null 2>&1; $PY $W/recover55.py $S/matches.txt --names $S/map_order.txt 2>&1 | grep -v Warning | head -1; done
echo "RECOVER DONE"
