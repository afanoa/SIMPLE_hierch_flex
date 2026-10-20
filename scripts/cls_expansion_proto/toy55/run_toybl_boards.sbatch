#!/bin/bash
#SBATCH -J toybl_boards
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Binding with no conformational consequence: on toy bl<KDA>, (A) the plain 7x7 lattice and (B) the two-cell focal
# board at the known blob site; bound/unbound separability (score_ah with bound = conformer >= 201). env: KDA (13),
# NT (16), NPARTS (unset = one process, N = distributed board), FB (build)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; KDA=${KDA:-13}; TAG=bl$KDA; S=$W/sim$TAG; L=$W/lat_$TAG; export TOY_SMPD=1.464; log(){ echo "[$(date '+%F %T') ${TAG}b] $*"; }
n=0; while [ $n -lt 600 ] && [ ! -s $L/ptcl3D_r0.txt ]; do sleep 60; n=$((n+1)); done; [ -s $L/ptcl3D_r0.txt ] || { log "no toy"; exit 1; }
PROJ=$W/proj$TAG/1_fastposes/proj$TAG.simple; GT=$S/gt_particles.txt; CL=$S/closed_ids.txt; SITE=$W/maps160/blob$KDA/blob_site_mask.mrc
NALL=$(grep -vc "^#" $GT); seq 1 $NALL > $L/subset_all.txt
score(){ $PY $W/score_ah.py $1/klem_latent.txt $GT $CL --holo0 201 --tag "$2" 2>&1 | grep -v Warning | sed 's/holo  /bound /'; }
export FB=${FB:-$HOME/SIMPLE_klem2/build} NTHR=${NT:-16}
log "A: plain 7x7 board ($KDA kDa blob, site unknown to the board)"; [ -s $L/fA/klem_latent.txt ] || INIT=random NFREE=2 bash $W/fboard.sh $L/fA $PROJ $L/subset_all.txt "20 20 20 14 14 14 10 10 10" "1.5 1.5 1.2 1.0 1.0 0.8 0.7 0.7 0.7" || exit 1
score $L/fA "A plain board, $KDA kDa"
log "B: focal 1x2 board at the known blob site"; [ -s $L/fB/klem_latent.txt ] || INIT=random NFREE=0 N1=1 N2=2 FOCUSMASK=$SITE bash $W/fboard.sh $L/fB $PROJ $L/subset_all.txt "14 14 12 12 10 10 8 8 8" "0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7 0.7" || exit 1
score $L/fB "B known-site focal 1x2, $KDA kDa"
log DONE
