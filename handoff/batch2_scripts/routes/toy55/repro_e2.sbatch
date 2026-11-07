#!/bin/bash
#SBATCH -J klem_repro_e2
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# refinement-only reproducibility: rflex13's stage A -> stage-B references -> stage B of node f692st r1 rerun from that node's OWN focused
# flex states and references (fflex_states.txt, frefs/), same keys and build (focus3), BORROW none as in f692st; only refine3D's internal
# random seed differs. Writes rnasep_klem/repro_r1/e2_<REP>/{stageA,stageB}_states.txt for comparison with the original. env REP (a, b, ...)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=4
K_=/mnt/beegfs/elmlund/afan/rnasep_klem; D=$K_/rflex_f692st/node_r1; O=$K_/repro_r1/e2_${REP:?}; MD=232; mkdir -p $O/stageA $O/stageB
log(){ echo "[$(date '+%F %T') repro_e2 $REP] $*"; }
st(){ simple_exec prg=print_project_field projfile=$1 oritype=ptcl3D mkdir=no 2>/dev/null | awk '{s = 0; for (i = 1; i <= NF; i++) if ($i ~ /^state=/) {split($i, a, "="); s = a[2] + 0}; print s}'; }
FKS=$(awk '$1 > 0' $D/fflex_states.txt | sort -n | uniq | wc -l); FVOLS=""; for j in $(seq 1 $FKS); do FVOLS="$FVOLS vol$j=$D/frefs/ref_state$(printf %02d $j).mrc"; done
if [ -z "$(ls $O/stageA/*_refine3D/recvol_state01_even_unfil.mrc 2>/dev/null)" ]; then cd $O/stageA && cp -f $D/n.simple r.simple
  simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$D/fflex_states.txt mkdir=no > rsel.log 2>&1 || { log "selection FAILED"; exit 1; }
  log "stage A refine3D, $FKS states (seeds: the node's focused flex states)"
  simple_exec prg=refine3D projfile=r.simple $FVOLS nstates=$FKS refine=prob_neigh prob_neigh_mode=geom pgrp=c1 mskdiam=$MD lp=15 maxits=3 center=no rec_backend=pcg \
    pcg_mskfile=$D/rsmask.mrc pcg_focus=yes nspace=2500 box_crop=64 smpd_crop=$(awk 'BEGIN {printf "%.8f", 160 * 1.464 / 64}') nparts=8 nthr=8 > ref.log 2>&1 || { log "stage A FAILED: $(grep -a -m2 ERROR ref.log | cut -c1-150)"; exit 1; }; fi
RDA=$(ls -td $O/stageA/*_refine3D | head -1); st $RDA/r.simple > $O/stageA_states.txt; log "stage A done: $(awk '$1 > 0' $O/stageA_states.txt | sort -n | uniq -c | awk '{printf "%d ", $1}')"
cd $O/stageB; [ -s b.simple ] || cp -f $RDA/r.simple b.simple
if [ -z "$(ls $O/stageB/*_reconstruct3D/recvol_state01.mrc 2>/dev/null)" ]; then
  simple_exec prg=reconstruct3D projfile=b.simple pgrp=c1 mskdiam=$MD nstates=$FKS nparts=8 nthr=8 > rec.log 2>&1 || { log "stage B references FAILED"; exit 1; }; fi
RRB=$(ls -td $O/stageB/*_reconstruct3D | head -1); VB=""; for q in $(seq 1 $FKS); do VB="$VB vol$q=$RRB/recvol_state$(printf %02d $q).mrc"; done
if [ -z "$(ls $O/stageB/*_refine3D/recvol_state01_even_unfil.mrc 2>/dev/null)" ]; then log "stage B refine3D"
  simple_exec prg=refine3D projfile=b.simple $VB nstates=$FKS refine=prob_neigh prob_neigh_mode=geom pgrp=c1 mskdiam=$MD lp=10 maxits=3 center=no rec_backend=pcg \
    pcg_mskfile=$D/rsmask.mrc pcg_focus=yes nspace=5000 box_crop=88 smpd_crop=2.66181803 nparts=8 nthr=8 > ref.log 2>&1 || { log "stage B FAILED: $(grep -a -m2 ERROR ref.log | cut -c1-150)"; exit 1; }; fi
RDB=$(ls -td $O/stageB/*_refine3D | head -1); st $RDB/b.simple > $O/stageB_states.txt; log "stage B done: $(awk '$1 > 0' $O/stageB_states.txt | sort -n | uniq -c | awk '{printf "%d ", $1}')"
log DONE
