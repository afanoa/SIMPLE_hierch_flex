#!/bin/bash
#SBATCH -J toy10_levers_anneal
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Lever 1: the deterministic-annealing split criterion for the expansion. The parent (klem_junk recipe) is re-run so
# that its even/odd cell deviations are saved; anneal_split.py then splits it recursively (a family splits while the
# reproducible variance along its principal mode exceeds tau x the noise along it), the children run at the realigned
# poses with the feedback recipe and are merged; the same parent split two ways by the spectral bisection with the same
# child recipe is the control. env: TAU (1), RT (realigned pose table)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; cd $W; log(){ echo "[$(date '+%F %T') an] $*"; }
RT=${RT:-$W/lat/realign_realign1/ptcl3D_realigned.txt}; CH="NFREE=0 TABLE=$RT PMARG=3 DEVW=0.05 PNOISE=1 LPS=20 20 20 14 14 14 10 10 10"
log "parent re-run (klem_junk recipe) with the half deviations saved"
env KTAG=_junk_h INIT=pca NSUB=0 NFREE=2 bash $W/run_klem.sbatch 2>&1 | grep "final:\|states: cross-half\|Error\|Traceback" | sed 's/^/[parent] /'
$PY $W/state_coverage.py $W/lat/klem_junk_h/klem_latent.txt $W/sim/gt_particles.txt --argmax --tag "parent re-run" | tail -1
log "annealing split (tau ${TAU:-1}) -> children at the realigned poses (feedback recipe) -> merge"
env PARENT=$W/lat/klem_junk_h ETAG=_anneal SPLITMODE=anneal TAU=${TAU:-1} NFREE=0 TABLE=$RT PMARG=3 DEVW=0.05 PNOISE=1 LPS="20 20 20 14 14 14 10 10 10" SIGMAS="1.5 1.5 1.2 1.0 1.0 0.8 0.7 0.7 0.7" bash $W/run_expand_toy.sbatch 2>&1 | sed 's/^/[anneal] /'
log "control: spectral 2-way split of the same parent, same child recipe"
env PARENT=$W/lat/klem_junk_h ETAG=_spectral2 NFAM=2 NFREE=0 TABLE=$RT PMARG=3 DEVW=0.05 PNOISE=1 LPS="20 20 20 14 14 14 10 10 10" SIGMAS="1.5 1.5 1.2 1.0 1.0 0.8 0.7 0.7 0.7" bash $W/run_expand_toy.sbatch 2>&1 | sed 's/^/[spectral2] /'
log DONE
