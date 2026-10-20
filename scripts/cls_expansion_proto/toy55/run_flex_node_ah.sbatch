#!/bin/bash
#SBATCH -J toyah_flex_node
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# flex_pca as the node chart on the apo/holo toy's family 1: at the REFINED poses (refine3D_states output) vs at the
# consensus/fast poses; latent scored on the 18-conformer truth (kNN, d', closed-vs-open). env: NP (10)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; B=$HOME/SIMPLE_kernel/snap/5e1e6aeee
export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel TOY_SMPD=1.464 OMP_NUM_THREADS=8
T=$A/rnasep_klem/toyah_tree; ST=$T/fam1/states.txt; GT=$W/simah/gt_particles.txt; CONS=$W/maps160/toy_consensus.mrc; R=$A/rnasep_klem/toyah_flex_node; mkdir -p $R; log(){ echo "[$(date '+%F %T') fnode] $*"; }
STATES=1,10,18,27,28,37,46,48,49,51,101,102,103,104,105,106,107,108; AXIS=0,0,1,0,0,0,0,0,0,1,0,1,0,0,1,1,0,0   # axis = closed indicator
awk '{print ($1>0)?1:0}' $ST > $R/member.txt   # one population: family membership only (flex_pca refuses several state labels)
for v in "refined:$(ls -td $T/fam1/refine/*_refine3D_states | head -1)/p.simple" "consensus:$T/src.simple"; do tag=${v%%:*}; PJ=${v#*:}; d=$R/$tag; mkdir -p $d && cd $d && cp -f $PJ p.simple
  simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$R/member.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; continue; }
  t0=$(date +%s); log "$tag poses: flex_pca on family 1 (neigs 8, lp 20, nparts ${NP:-10})"
  simple_exec prg=flex_pca projfile=p.simple vol1=$CONS mskdiam=200 lp=20 neigs=8 npreimages=2 rec_states=no nparts=${NP:-10} nthr=8 > flex.log 2>&1 || { log "flex_pca FAILED: $(grep -a -m2 'ERROR\|THROW' flex.log | cut -c1-200)"; continue; }
  FD=$(ls -td $d/*_flex_pca | head -1); log "$tag: flex done in $(( ($(date +%s) - t0) / 60 )) min"
  echo "== node chart at $tag poses"; $PY $W/ptcl_latent_eval.py $FD $GT --states $STATES --axis=$AXIS --subset $R/member.txt --tag node_$tag 2>&1 | grep -v Warning | grep "knn\|closed\|ladder\|states:" | head -6
done; log "ALL DONE"
