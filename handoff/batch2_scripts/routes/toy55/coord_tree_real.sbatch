#!/bin/bash
#SBATCH -J klem_coord
#SBATCH -p csbdevel
#SBATCH -c 64
#SBATCH --mem=240G
#SBATCH -t 08:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# stage 4 on real data, one family of a run: the family's merge-up states relabelled 1..M -> (1) their maps at the family-refined
# poses; (2) their tree; (3) ONE cross-half scoring pass of the family's particles against all M maps (nspace 5000, cropped box); (4) per
# tree node the PC1 coordinate, continuity test and even-width sampling (coord_tree.py, closed-arm of the member states as the check);
# (5) the sampled nodes' even bins reconstructed and read out (closed-arm / open-arm occupancy per bin, the trajectory along the branch).
# env RTAG (f692r), FAM (1|2), CDIR (mergeup/consolidated_family: states.txt over master rows + arm_states.txt), OUT, SNAP (build with
# more than vol20: focus3), MSKDIAM (232), LP (8), MINP (10000)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=${SNAP:-$HOME/SIMPLE_focus/snap/focus3}; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem; R=$K_/${RTAG:?}; F=${FAM:?}
C=$R/${CDIR:-mergeup/consolidated_family}; O=${OUT:-$R/coordtree_fam$F}; MD=${MSKDIAM:-232}; mkdir -p $O/rec $O/score && cd $O; log(){ echo "[$(date '+%F %T') coord-tree $RTAG fam$F] $*"; }
BASE=$(ls -td $R/fam$F/refine/*_refine3D_states | head -1)/p.simple; U=$O/units.txt
[ -s $U ] || $PY - <<PYEOF
import numpy as np
st = np.loadtxt('$C/states.txt', dtype=int); sub = np.loadtxt('$R/expand/family${F}_subset.txt', dtype=int).ravel() - 1; m = np.zeros(len(st), bool); m[sub] = True
ids = sorted(int(x) for x in np.unique(st[m & (st > 0)])); u = np.zeros(len(st), int)
for i, s in enumerate(ids, 1): u[m & (st == s)] = i
np.savetxt('$U', u, fmt='%d'); open('$O/names.txt', 'w').write(''.join('%d %d\n' % (i, s) for i, s in enumerate(ids, 1)))
print('family $F: %d merge-up states, %d particles' % (len(ids), (u > 0).sum()))
PYEOF
M=$(awk '$1 > 0' $U | sort -n | uniq | wc -l); log "$M states, poses from $BASE"
cd $O/rec; if [ -z "$(ls -d $O/rec/*_reconstruct3D 2>/dev/null)" ]; then cp -f $BASE r.simple
  simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$U mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
  simple_exec prg=update_project projfile=r.simple qsys_name=local nparts=4 nthr=16 mkdir=no > up.log 2>&1
  simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=$MD nstates=$M nparts=4 nthr=16 > rec.log 2>&1 || { log "reconstruct3D FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-150)"; exit 1; }; fi
RR=$(ls -td $O/rec/*_reconstruct3D | head -1); log "maps in $RR"
[ -s $O/tree/merge_tree.txt ] || $PY $T/state_merge.py $O/tree --rdir $RR --labels $U --link complete --ncut 2 2>&1 | grep -v Warning | grep "states,\|leaf pair" | cut -c1-160
cd $O/score; if [ -z "$(ls $O/score/*_refine3D/dist_part_neigh_*.dat 2>/dev/null)" ]; then cp -f $O/rec/r.simple s.simple
  XR=$O/score/xref; mkdir -p $XR; VOLS=""
  for s in $(seq 1 $M); do b=recvol_state$(printf %02d $s); cp -f $RR/$b.mrc $XR/$b.mrc; cp -f $RR/${b}_odd.mrc $XR/${b}_even.mrc; cp -f $RR/${b}_even.mrc $XR/${b}_odd.mrc; VOLS="$VOLS vol$s=$XR/$b.mrc"; done
  log "scoring the family's particles against $M maps (cross-half, lp ${LP:-8}, nspace 5000, box 88)"
  simple_exec prg=refine3D projfile=s.simple $VOLS nstates=$M refine=prob_neigh prob_neigh_mode=geom pgrp=c1 mskdiam=$MD lp=${LP:-8} maxits=1 center=no nspace=5000 box_crop=88 smpd_crop=2.66181803 nparts=4 nthr=16 > score.log 2>&1 || { log "scoring FAILED: $(grep -a -m2 ERROR score.log | cut -c1-150)"; exit 1; }
  SD=$(ls -td $O/score/*_refine3D | head -1); rm -f $SD/*.mrc; rm -rf $XR; fi
SD=$(ls -td $O/score/*_refine3D | head -1)
$PY $T/coord_tree.py $O/nodes --tabs "$SD/dist_part_neigh_*.dat" --nstates $M --units $U --order $(seq -s, 1 $M) --tree $O/tree/merge_tree.txt --minp ${MINP:-10000} --names $O/names.txt --arm $C/arm_states.txt 2>&1 | grep -v Warning
# (5) the sampled branches: their even bins reconstructed and read out
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
for nd in $(cat $O/nodes/sampled_nodes.txt 2>/dev/null); do D=$O/nodes/node$(printf %03d $nd); [ -s $D/arm_states.txt ] && continue; mkdir -p $D/rec && cd $D/rec
  NB=$(awk '$1 > 0' $D/even_bins.txt | sort -n | uniq | wc -l); cp -f $O/rec/r.simple b.simple
  simple_exec prg=selection projfile=b.simple oritype=ptcl3D infile=$D/even_bins.txt mkdir=no > sel.log 2>&1 && simple_exec prg=update_project projfile=b.simple qsys_name=local nparts=4 nthr=16 mkdir=no > up.log 2>&1 \
    && simple_exec prg=reconstruct3D projfile=b.simple pgrp=c1 mskdiam=$MD nstates=$NB nparts=4 nthr=16 > rec.log 2>&1 || { log "node $nd bins reconstruct FAILED"; continue; }
  BR=$(ls -td $D/rec/*_reconstruct3D | head -1); mkdir -p $D/maps/halves; for s in $(seq 1 $NB); do b=recvol_state$(printf %02d $s); cp -f $BR/$b.mrc $D/maps/state_$(printf %02d $s).mrc; cp -f $BR/${b}_even.mrc $D/maps/halves/state_$(printf %02d $s)_even.mrc; cp -f $BR/${b}_odd.mrc $D/maps/halves/state_$(printf %02d $s)_odd.mrc; done
  $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $D/arm_states.txt $D/maps 2>&1 | grep -v Warning | tail -1
  echo "node $nd sampled trajectory (bin: closed-arm open-arm, particles): $(awk '!/^#/ {n = split($1, p, "_"); b = p[n]; sub(/\.mrc/, "", b); printf "%s:%.2f/%.2f ", b, $2, $3}' $D/arm_states.txt) | particles $(awk '$1 > 0 {c[$1]++} END {for (b = 1; b <= 20; b++) if (c[b]) printf "%d ", c[b]}' $D/even_bins.txt)"
done
log DONE
