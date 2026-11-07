#!/bin/bash
#SBATCH -J klem_coord
#SBATCH -p csbdevel
#SBATCH -c 64
#SBATCH --mem=120G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# coord_loop.sbatch = lrec_loop.sbatch (realistic-SNR session) with OUT/TAG/CODES parameters, cross-half scoring by default and
# coord_pc.py (even-width sampling of the final coordinate). Called per node by coord_nodes.sbatch.
# continuous-coordinate leaves from labelled particles (a K-state classification): reconstruct the labels' maps -> score every
# labelled particle against them (refine3D, 1 iteration, prob_neigh geom; its pose/state output is discarded) -> per-particle best
# dist to each map -> PC1 of the centred score vectors (lrec_pc.py) -> KL quantile bins = next labels -> ... NIT times; the last
# bins are reconstructed as the leaves (final/). Poses stay those of BASE throughout.
# env SRC (lrec30 run label: BASE = its refined project, LAB0 = its final_states.txt), LABEL, BASE, LAB0, KL (8), NIT (2), LP (8),
# BOXC SMPDC (88 / 2.66181803: the classification's crop; full box needs >300 GB), XHALF (1: score each half against the OTHER half's maps = copies with _even/_odd swapped; no self-inclusion), EXTRA (refine3D keys)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; KD=$A/rnasep_klem
TAG=${TAG:?}; GT=${GT:-$T/sim$TAG/gt_particles.txt}; O=${OUT:?}; LABEL=$(basename $O); KL=${KL:-8}; NIT=${NIT:-2}; LP=${LP:-8}
mkdir -p $O && cd $O; log(){ echo "[$(date '+%F %T') lrec_loop $LABEL] $*"; }
[ -s "${BASE:?}" ] && [ -s "${LAB0:?}" ] || { log "missing BASE $BASE or LAB0 $LAB0"; exit 1; }
cp -f $BASE base.simple; cp -f $LAB0 states_it0.txt; NROWS=$(wc -l < states_it0.txt)
recon(){ # $1 dir, $2 labels -> sets VOLS, K
  mkdir -p $1 && cd $1; cp -f $O/base.simple r.simple; K=$(awk '$1 > 0' $2 | sort -n | uniq | wc -l)
  simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$2 mkdir=no > rsel.log 2>&1 || { log "selection FAILED in $1"; exit 1; }
  simple_exec prg=update_project projfile=r.simple qsys_name=local nparts=4 nthr=16 mkdir=no > rup.log 2>&1
  simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=220 nstates=$K nparts=4 nthr=16 > rec.log 2>&1 || { log "reconstruct3D FAILED in $1"; exit 1; }
  RR=$(ls -td $1/*_reconstruct3D | head -1); rm -f $RR/*_part*_even.mrc $RR/*_part*_odd.mrc $RR/rho_*
  VOLS=""; for s in $(seq 1 $K); do VOLS="$VOLS vol$s=$RR/recvol_state$(printf %02d $s).mrc"; done; }
for it in $(seq 0 $NIT); do
  W=$O/it$it; recon $W $O/states_it$it.txt
  [ $it = 0 ] && [ -s $RR/r.simple ] && cp -f $RR/r.simple $O/base.simple   # keep the e/o partition the first reconstruction made (raw projects have none)
  cd $W; cp -f $O/base.simple s.simple
  simple_exec prg=selection projfile=s.simple oritype=ptcl3D infile=$O/states_it$it.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED it $it"; exit 1; }
  simple_exec prg=update_project projfile=s.simple qsys_name=local nparts=4 nthr=16 mkdir=no > up.log 2>&1
  if [ "${XHALF:-1}" = 1 ]; then XR=$W/xref; mkdir -p $XR; XV=""; for v in $VOLS; do f=${v#*=}; b=$(basename $f .mrc); cp -f $f $XR/$b.mrc; cp -f ${f%.mrc}_odd.mrc $XR/${b}_even.mrc; cp -f ${f%.mrc}_even.mrc $XR/${b}_odd.mrc; XV="$XV ${v%%=*}=$XR/$b.mrc"; done; VOLS=$XV; fi
  log "it $it: scoring all labelled particles against $K maps (lp $LP, cross-half ${XHALF:-1})"
  simple_exec prg=refine3D projfile=s.simple $VOLS nstates=$K refine=prob_neigh prob_neigh_mode=geom pgrp=c1 mskdiam=220 lp=$LP maxits=1 center=no box_crop=${BOXC:-88} smpd_crop=${SMPDC:-2.66181803} nparts=4 nthr=16 ${EXTRA:-} > score.log 2>&1 || { log "scoring refine3D FAILED: $(grep -a -m2 'ERROR' score.log | cut -c1-150)"; exit 1; }
  SD=$(ls -td $W/*_refine3D | head -1)
  $PY $T/coord_pc.py $W/pc --tabs "$SD/dist_part_neigh_*.dat" --nstates $K --nrows $NROWS --k $KL --gt $GT --npc 3 --codes ${CODES:-swing:200:1} --even ${EVEN:-8} --minb ${MINB:-1000} $([ $it -gt 0 ] && echo --prev $O/it$((it - 1))/pc/coord.txt) 2>&1 | sed "s/^/  it $it /"
  [ -s $W/pc/bins.txt ] || { log "lrec_pc.py FAILED it $it"; exit 1; }
  rm -f $SD/*_part*_even.mrc $SD/*_part*_odd.mrc $SD/rho_* $SD/*.mrc; rm -rf $W/xref; cp -f $W/pc/bins.txt $O/states_it$((it + 1)).txt
done
recon $O/final $O/states_it$((NIT + 1)).txt; cp -f $O/states_it$((NIT + 1)).txt $O/final_states.txt; cp -f $O/it$NIT/pc/coord.txt $O/final_coord.txt; cp -f $O/it$NIT/pc/even_bins.txt $O/final_even_bins.txt
log "DONE: leaves = $O/final (maps), labels final_states.txt, coordinate final_coord.txt"
