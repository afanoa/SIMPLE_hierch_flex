#!/bin/bash
#SBATCH -J klem_floor_sweep
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# how many of the level-1 leaves (and zoomed sub-states) are the same state at coarser floors (sphere mask, band FSC)
set -uo pipefail; PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; A=/mnt/beegfs/elmlund/afan/rnasep_klem; export MASK=sphere WHITEN=0
for g in "region*" "zoomregion*"; do for fl in 10 12 15 20; do export LEAFGLOB=$g FLOOR=$fl; echo "== $g floor $fl A"; $PY $W/merge_data_derived.py $A/f692c $W/maps160/toy_consensus.mrc $A/f692c/merge_sweep_$(echo $g | tr -d '*')_$fl 2>&1 | grep -E "^[0-9]+ leaf|pairs at the start|merges ->|remaining"; done; done; echo "SWEEP DONE"
