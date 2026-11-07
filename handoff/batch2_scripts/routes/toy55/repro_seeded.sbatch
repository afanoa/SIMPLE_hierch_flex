#!/bin/bash
#SBATCH -J klem_repro_seeded
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# pole-seeded refinement test for node f692st r1: seeds = a per-row label file (SEEDS: 1..P over the master rows, 0 outside r1; e.g. each
# particle's best pole from one run's reproducible poles), references reconstructed from the seeds, then rflex13's stage A (same node
# project n.simple = the shared input poses, mask rsmask.mrc, lp 15, box 64, nspace 2500, pcg_focus, greedy projections) with
# nstates = P for MAXITS iterations. Seeds from different runs' poles -> cross-run agreement of the refined partitions and pole maps.
# Writes rnasep_klem/repro_r1/seeded_<TAG>/stageA_states.txt. env TAG SEEDS [MAXITS]
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=4
K_=/mnt/beegfs/elmlund/afan/rnasep_klem; D=$K_/rflex_f692st/node_r1; MD=232; O=$K_/repro_r1/seeded_${TAG:?}; mkdir -p $O/refs $O/stageA
log(){ echo "[$(date '+%F %T') repro_seeded $TAG] $*"; }
st(){ simple_exec prg=print_project_field projfile=$1 oritype=ptcl3D mkdir=no 2>/dev/null | awk '{s = 0; for (i = 1; i <= NF; i++) if ($i ~ /^state=/) {split($i, a, "="); s = a[2] + 0}; print s}'; }
cp -f ${SEEDS:?} $O/seeds.txt; NK=$(awk '$1 > 0' $O/seeds.txt | sort -n | uniq | wc -l); log "$NK seed states: $(awk '$1 > 0' $O/seeds.txt | sort -n | uniq -c | awk '{printf "%d ", $1}')"
[ $(awk '$1 > 0' $O/seeds.txt | sort -n | tail -1) -eq $NK ] || { log "seed labels must be 1..P"; exit 1; }
if [ -z "$(ls $O/refs/*_reconstruct3D/recvol_state$(printf %02d $NK).mrc 2>/dev/null)" ]; then cd $O/refs && cp -f $D/n.simple s.simple
  simple_exec prg=selection projfile=s.simple oritype=ptcl3D infile=$O/seeds.txt mkdir=no > sel.log 2>&1 && \
  simple_exec prg=reconstruct3D projfile=s.simple pgrp=c1 mskdiam=$MD nstates=$NK nparts=4 nthr=8 > rec.log 2>&1 || { log "seed references FAILED"; exit 1; }; fi
RR=$(ls -td $O/refs/*_reconstruct3D | head -1); FV=""; for j in $(seq 1 $NK); do FV="$FV vol$j=$RR/recvol_state$(printf %02d $j).mrc"; done
if [ -z "$(ls $O/stageA/*_refine3D/recvol_state01_even_unfil.mrc 2>/dev/null)" ]; then cd $O/stageA && cp -f $D/n.simple r.simple
  simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$O/seeds.txt mkdir=no > rsel.log 2>&1 || { log "selection FAILED"; exit 1; }
  log "stage A refine3D, $NK states, maxits ${MAXITS:-6}, greedy projections"
  simple_exec prg=refine3D projfile=r.simple $FV nstates=$NK refine=prob_neigh prob_neigh_mode=geom prob_athres=0.01 pgrp=c1 mskdiam=$MD lp=15 maxits=${MAXITS:-6} \
    center=no rec_backend=pcg pcg_mskfile=$D/rsmask.mrc pcg_focus=yes nspace=2500 box_crop=64 smpd_crop=3.66000000 nparts=4 nthr=8 > ref.log 2>&1 || { log "stage A FAILED: $(grep -a -m2 'ERROR\|THROW' ref.log | cut -c1-200)"; exit 1; }; fi
RDA=$(ls -td $O/stageA/*_refine3D | head -1); st $RDA/r.simple > $O/stageA_states.txt
log "stage A done: $(awk '$1 > 0' $O/stageA_states.txt | sort -n | uniq -c | awk '{printf "%d ", $1}')| state overlap per iteration: $(grep -a 'STATE OVERLAP' $O/stageA/ref.log | awk '{printf "%s ", $NF}')"
log DONE
