#!/bin/bash
#SBATCH -J subcls_rootprep
#SBATCH -p csbdevel
#SBATCH -c 48
#SBATCH --mem=160G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Root node of a synthetic set in the rflex4 node layout, for the sub-class route (subcls_node2.sbatch RUN=rflex_root_<MODE> NODE=r):
# fresh project at the refined poses (focus3 build; stack sim<MODE>/all.mrcs, poses + CTF lat_<MODE>/ptcl3D_r0.txt) -> node mean
# (reconstruct3D) -> flex_pca on the node exactly as rflex4 runs it (lp 15, 8 states, 2 probe iterations), NEIGS covariance components,
# no state maps. Times each stage. env MODE (armsigm1), NEIGS (16), MSKDIAM (220)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=12
A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; MODE=${MODE:-armsigm1}; NEIGS=${NEIGS:-16}; MD=${MSKDIAM:-220}
R=$K_/rflex_root_$MODE; D=$R/node_r; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') rootprep $MODE] $*"; }
NR=$(grep -c 'e1=' $T/lat_$MODE/ptcl3D_r0.txt); t0=$(date +%s)
if [ ! -s $D/n.simple ]; then rm -rf $D/imp && mkdir -p $D/imp && cd $D/imp
  simple_exec prg=new_project projname=n > np.log 2>&1 && cd n && simple_exec prg=import_particles projfile=n.simple stk=$T/sim$MODE/all.mrcs oritab=$T/lat_$MODE/ptcl3D_r0.txt smpd=1.464 kv=300 cs=2.7 fraca=0.1 ctf=yes mkdir=no > import.log 2>&1 \
    && simple_exec prg=update_project projfile=n.simple qsys_name=local nparts=4 nthr=12 mkdir=no > up.log 2>&1 && cp -f n.simple $D/n.simple || { log "import FAILED: $(grep -a -m2 ERROR import.log | cut -c1-150)"; exit 1; }
  cd $D; awk -v n=$NR 'BEGIN {for (i = 1; i <= n; i++) print 1}' > $R/rows_r.txt; echo "r 0 $D/n.simple $R/rows_r.txt" > $R/queue.txt; fi
t1=$(date +%s)
if [ -z "$(ls $D/mean/*_reconstruct3D/recvol_state01.mrc 2>/dev/null)" ]; then mkdir -p $D/mean && cd $D/mean && cp -f $D/n.simple m.simple
  simple_exec prg=reconstruct3D projfile=m.simple pgrp=c1 mskdiam=$MD nstates=1 nparts=4 nthr=12 > rec.log 2>&1 || { log "mean FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-150)"; exit 1; }; cd $D; fi
MEAN=$(ls -td $D/mean/*_reconstruct3D | head -1)/recvol_state01.mrc; t2=$(date +%s)
if [ -z "$(ls $D/*_flex_pca/flex_pca_coordinates.txt 2>/dev/null)" ]; then cp -f $D/n.simple f.simple
  log "flex_pca: $NR particles, neigs $NEIGS, lp 15"
  simple_exec prg=flex_pca projfile=f.simple vol1=$MEAN mskdiam=$MD lp=15 neigs=$NEIGS npreimages=8 rec_states=no min_state_frac=0 umap=no n_probe_iters=2 nparts=4 nthr=12 > flex.log 2>&1 || { log "flex_pca FAILED: $(grep -a -m2 'ERROR\|THROW' flex.log | cut -c1-200)"; exit 1; }; fi
t3=$(date +%s); log "timing: import $((t1 - t0)) s, mean $((t2 - t1)) s, flex_pca $((t3 - t2)) s; coordinates $(ls -td $D/*_flex_pca | head -1)/flex_pca_coordinates.txt"
log DONE
