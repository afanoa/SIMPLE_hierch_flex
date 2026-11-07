#!/bin/bash
#SBATCH -J klem_repro_var
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# search-variant reproducibility test for node f692st r1: rflex13's stage A (same node project n.simple, mask rsmask.mrc, lp 15, box 64,
# nspace 2500, pcg_focus, focus3) from START S (the script's focused flex states + frefs) or P (the port run's focused flex states + maps),
# with the variant's search: REFINE (prob_neigh), MODE (prob_neigh_mode, geom), KEYS (extra, e.g. prob_athres=0.01), MAXITS (3).
# Two starts per variant -> cross-start agreement. Writes rnasep_klem/repro_r1/var_<TAG>_<START>/stageA_states.txt. env TAG START
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=4
K_=/mnt/beegfs/elmlund/afan/rnasep_klem; D=$K_/rflex_f692st/node_r1; P1=$K_/stree_port/tree_f692r1/1_hetero_tree/node_r
O=$K_/repro_r1/var_${TAG:?}_${START:?}; mkdir -p $O/refs $O/stageA; RF=${REFINE:-prob_neigh}
log(){ echo "[$(date '+%F %T') repro_var $TAG $START] $*"; }
st(){ simple_exec prg=print_project_field projfile=$1 oritype=ptcl3D mkdir=no 2>/dev/null | awk '{s = 0; for (i = 1; i <= NF; i++) if ($i ~ /^state=/) {split($i, a, "="); s = a[2] + 0}; print s}'; }
FV=""; for j in $(seq 1 8); do jj=$(printf %02d $j)
  if [ $START = S ]; then src=$D/frefs/ref_state$jj; else src=$P1/fflex/recvol_state$jj; fi
  for x in "" _even _odd; do cp -f $src$x.mrc $O/refs/ref_state$jj$x.mrc || { log "start map $src$x.mrc missing"; exit 1; }; done; FV="$FV vol$j=$O/refs/ref_state$jj.mrc"; done
if [ $START = S ]; then SEEDS=$D/fflex_states.txt; else SEEDS=$K_/repro_r1/labels/P_fflex.txt; fi
XK=""; [ $RF = prob_neigh ] && XK="prob_neigh_mode=${MODE:-geom}"
if [ -z "$(ls $O/stageA/*_refine3D/recvol_state01_even_unfil.mrc 2>/dev/null)" ]; then cd $O/stageA && cp -f $D/n.simple r.simple
  simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$SEEDS mkdir=no > rsel.log 2>&1 || { log "selection FAILED"; exit 1; }
  log "stage A refine3D: refine=$RF $XK ${KEYS:-} maxits ${MAXITS:-3}"
  simple_exec prg=refine3D projfile=r.simple $FV nstates=8 refine=$RF $XK pgrp=c1 mskdiam=232 lp=15 maxits=${MAXITS:-3} center=no rec_backend=pcg \
    pcg_mskfile=$D/rsmask.mrc pcg_focus=yes nspace=2500 box_crop=64 smpd_crop=3.66000000 ${KEYS:-} nparts=4 nthr=6 > ref.log 2>&1 || { log "stage A FAILED: $(grep -a -m2 'ERROR\|THROW' ref.log | cut -c1-200)"; exit 1; }; fi
RDA=$(ls -td $O/stageA/*_refine3D | head -1); st $RDA/r.simple > $O/stageA_states.txt
log "stage A done: $(awk '$1 > 0' $O/stageA_states.txt | sort -n | uniq -c | awk '{printf "%d ", $1}')| state overlap per iteration: $(grep -a 'STATE OVERLAP' $O/stageA/ref.log | awk '{printf "%s ", $NF}')"
log DONE
