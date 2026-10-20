#!/bin/bash
#SBATCH -J klem_synth_perarm
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# per-arm (arms*) / per-site (tethered*) truth readout of a synthetic chain's leaf states (env MODE)
set -uo pipefail; PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; MODE=${MODE:?}
case $MODE in arms*) S=arm_angle_readout.py;; *) S=site_occupancy_readout.py;; esac
echo "== $MODE: $S over the leaf states of $A/rnasep_klem/synth_$MODE"; $PY $W/$S $A/rnasep_klem/synth_$MODE $W/sim$MODE/gt_particles.txt $W/truth_$MODE.txt 2>&1 | grep -v Warning | tail -30; echo "PERARM DONE $MODE"
