#!/bin/bash
#SBATCH -J rnasep_f692_family
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 20:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# One family of the f692 chain as its own job (so the families run concurrently): child board, k-means,
# refine3D_states, feedback board (f692_family.sh), then the parallel regions stage (f692_regions.sh). Writes
# $R/fam$FAM/DONE. env: R, FAM, SRC, METAALL, FB, NTHR, NP, MAXITS, LPSC, SIGC, MAXREG, NL2, PARMAX
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; log(){ echo "[$(date '+%F %T') famjob$FAM] $*"; }
export R SRC FB NTHR NP MAXITS LPSC SIGC NL2 PARMAX NPARTS
log "family $FAM stage"; FAM=$FAM META=$METAALL bash $W/f692_family.sh || { log "family stage FAILED"; exit 1; }
[ -s $R/fam$FAM/klem_fb/klem_latent.txt ] || { log "no feedback board"; exit 1; }
log "regions stage, family $FAM"; FAM=$FAM META=$R/fam$FAM/meta_fb.npz MAXREG=${MAXREG:-6} NL2=${NL2:-3} bash $W/f692_regions.sh || log "regions stage FAILED"
touch $R/fam$FAM/DONE; log DONE
