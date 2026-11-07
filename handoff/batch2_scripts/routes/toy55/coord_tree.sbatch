#!/bin/bash
#SBATCH -J klem_coord
#SBATCH -p csbdevel
#SBATCH -c 64
#SBATCH --mem=160G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# stage 4 of the general pipeline for every branch at once: (1) reconstruct the merge-up's final states; (2) their tree (state_merge.py,
# complete linkage); (3) ONE cross-half scoring pass of all particles against all final-state maps (refine3D 1 iteration, prob_neigh
# geom, the maps' even/odd swapped so each half meets the other half's maps, nspace 5000, cropped box); (4) coord_tree.py: per tree node
# the PC1 coordinate of the node's particles over the node's states' scores and its even-width sampling.
# env TAG, UNITS (final per-row units 1..M), BASE (project with the refined poses), OUT, LP (8), MINP (10000); GT/CODES (environment)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; TAG=${TAG:?}; O=${OUT:?}; U=${UNITS:?}; GT=${GT:-$T/sim$TAG/gt_particles.txt}
mkdir -p $O/rec $O/score && cd $O; log(){ echo "[$(date '+%F %T') coord-tree $TAG] $*"; }; M=$(awk '$1 > 0' $U | sort -n | uniq | wc -l)
# (1) final-state maps
cd $O/rec; if [ -z "$(ls -d $O/rec/*_reconstruct3D 2>/dev/null)" ]; then cp -f ${BASE:?} r.simple
  simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$U mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
  simple_exec prg=update_project projfile=r.simple qsys_name=local nparts=4 nthr=16 mkdir=no > up.log 2>&1
  simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=220 nstates=$M nparts=4 nthr=16 > rec.log 2>&1 || { log "reconstruct3D FAILED"; exit 1; }; fi
RR=$(ls -td $O/rec/*_reconstruct3D | head -1); log "$M final-state maps in $RR"
# (2) their tree
[ -s $O/tree/merge_tree.txt ] || $PY $T/state_merge.py $O/tree --rdir $RR --labels $U --link complete --ncut 2 2>&1 | grep -v Warning | grep "states,\|leaf pair" | cut -c1-160
# (3) one cross-half scoring pass of all particles against all final-state maps
cd $O/score; if [ -z "$(ls $O/score/*_refine3D/dist_part_neigh_*.dat 2>/dev/null)" ]; then cp -f $RR/r.simple s.simple 2>/dev/null || cp -f $O/rec/r.simple s.simple
  XR=$O/score/xref; mkdir -p $XR; VOLS=""
  for s in $(seq 1 $M); do b=recvol_state$(printf %02d $s); cp -f $RR/$b.mrc $XR/$b.mrc; cp -f $RR/${b}_odd.mrc $XR/${b}_even.mrc; cp -f $RR/${b}_even.mrc $XR/${b}_odd.mrc; VOLS="$VOLS vol$s=$XR/$b.mrc"; done
  log "scoring all particles against $M maps (cross-half, lp ${LP:-8}, nspace 5000, box 88)"
  simple_exec prg=refine3D projfile=s.simple $VOLS nstates=$M refine=prob_neigh prob_neigh_mode=geom pgrp=c1 mskdiam=220 lp=${LP:-8} maxits=1 center=no nspace=5000 box_crop=88 smpd_crop=2.66181803 nparts=4 nthr=16 > score.log 2>&1 || { log "scoring FAILED: $(grep -a -m2 ERROR score.log | cut -c1-150)"; exit 1; }
  SD=$(ls -td $O/score/*_refine3D | head -1); rm -f $SD/*.mrc; rm -rf $XR; fi
SD=$(ls -td $O/score/*_refine3D | head -1)
# (4) per tree node: coordinate and even sampling
$PY $T/coord_tree.py $O/nodes --tabs "$SD/dist_part_neigh_*.dat" --nstates $M --units $U --order $(seq -s, 1 $M) --tree $O/tree/merge_tree.txt --minp ${MINP:-10000} --gt $GT --codes ${CODES:-swing:200:1} 2>&1 | grep -v Warning
log DONE
