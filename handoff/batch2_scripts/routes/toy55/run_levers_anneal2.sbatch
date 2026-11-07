#!/bin/bash
#SBATCH -J toy10_levers_anneal2
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Lever 1, second round: the annealing split of the parent (8 families of 10-28k particles) starved 7x7 children (no
# reproducible cells). Same split, child boards sized to the family mass (~1500 particles per cell, 3x3 to 7x7), same
# realigned poses and feedback recipe, merged purity. env: RT (realigned pose table), PERCELL (1500), TAU (1)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; cd $W; log(){ echo "[$(date '+%F %T') an2] $*"; }
RT=${RT:-$W/lat/realign_realign1/ptcl3D_realigned.txt}
log "annealing split of klem_junk_h -> children sized to their families -> merge"
env PARENT=$W/lat/klem_junk_h ETAG=_anneal_g SPLITMODE=anneal TAU=${TAU:-1} AUTOGRID=1 PERCELL=${PERCELL:-1500} NFREE=0 TABLE=$RT PMARG=3 DEVW=0.05 PNOISE=1 LPS="20 20 20 14 14 14 10 10 10" SIGMAS="1.5 1.5 1.2 1.0 1.0 0.8 0.7 0.7 0.7" bash $W/run_expand_toy.sbatch 2>&1 | sed 's/^/[anneal_g] /'
log DONE
