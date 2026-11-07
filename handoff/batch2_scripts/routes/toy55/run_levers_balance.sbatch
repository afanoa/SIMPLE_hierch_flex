#!/bin/bash
#SBATCH -J toy10_levers_balance
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Lever 3: the OT-balanced E-step (Sinkhorn): the parent (klem_junk recipe, consensus poses) with the grid cells
# constrained to equal shares of the grid mass, fully (lambda 1) and half-way (0.5), against the plain parent; coverage
# per conformer, occupancy, junk capture; then the density-mode reading of each. env: BALS ("1 0.5")
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; cd $W; log(){ echo "[$(date '+%F %T') ba] $*"; }
for B in ${BALS:-1 0.5}; do
  KT=_bal${B/./}; log "parent with the balanced E-step, lambda $B"
  env KTAG=$KT INIT=pca NSUB=0 NFREE=2 BAL=$B bash $W/run_klem.sbatch 2>&1 | grep "balanced E-step\|it [0-9] (lp\|final:\|states: cross-half\|free1\|Error\|Traceback" | sed "s/^/[bal $B] /"
  $PY $W/state_coverage.py $W/lat/klem$KT/klem_latent.txt $W/sim/gt_particles.txt --argmax --tag "balanced lambda $B" | tail -1
done
$PY $W/state_coverage.py $W/lat/klem_junk/klem_latent.txt $W/sim/gt_particles.txt --argmax --tag "plain parent (klem_junk)" | tail -1
for B in ${BALS:-1 0.5}; do echo "== density modes, lambda $B"; $PY $W/density_modes.py $W/lat/klem_bal${B/./} --gt $W/sim/gt_particles.txt 2>&1 | grep -v Warning; done
log DONE
