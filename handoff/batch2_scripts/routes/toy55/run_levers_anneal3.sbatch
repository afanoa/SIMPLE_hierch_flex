#!/bin/bash
#SBATCH -J toy10_levers_anneal3
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Lever 1, third round: the annealed families (8, 10-28k particles) with the recipe that let the level-3 grandchildren
# survive on ~800 particles per cell: focus mask on, fine band 12->8 A, pose marginalisation, boards sized to the family
# (~1000 particles per cell), realigned poses; merged purity at 0.5 and 0.3 survival. env: RT, PERCELL (1000), TAU (1)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; cd $W; log(){ echo "[$(date '+%F %T') an3] $*"; }
RT=${RT:-$W/lat/realign_realign1/ptcl3D_realigned.txt}
log "annealing split of klem_junk_h -> focused fine-band children sized to their families -> merge"
env PARENT=$W/lat/klem_junk_h ETAG=_anneal_f SPLITMODE=anneal TAU=${TAU:-1} AUTOGRID=1 PERCELL=${PERCELL:-1000} NFREE=0 TABLE=$RT FOCUS=1 PMARG=3 DEVW=0.05 PNOISE=1 LPS="12 12 10 10 8 8 8 8 8" SIGMAS="1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7" bash $W/run_expand_toy.sbatch 2>&1 | sed 's/^/[anneal_f] /'
CH=$(ls -d $W/lat/klem_expand_anneal_f_fam* | sort -V | tr '\n' ' '); echo "== merge at survival 0.3"; $PY $W/expand_families.py merge $W/lat/expand_anneal_f $CH --gt $W/sim/gt_particles.txt --minrep 0.3 | grep 'MERGED\|best-state' | sed 's/^/[anneal_f 0.3] /'
log DONE
