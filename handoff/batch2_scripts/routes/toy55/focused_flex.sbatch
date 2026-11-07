#!/bin/bash
#SBATCH -J klem_focused_flex
#SBATCH -p norm
#SBATCH -c 64
#SBATCH --mem=240G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# mask-focused flex of one real rflex node: flex_pca on the node's particles (parent's poses) with pcg_mskfile = the node's real-space variation
# mask (focused_arm_<RT>_<NODE>/arm_mask.mrc) -> the latent describes the arm only. Readouts: separability of an earlier run's closed states
# (closed_latent.py: LDA AUC, PCs, ICs) and the cheap route on the focused latent (view-centred PC1: 10 quantile + 12 equal-width bins)
# reconstructed at the parent's poses with the closed-arm readout. env RT (f692rf), NODE (r22), LPF (15), MSKDIAM (232)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; RT=${RT:-f692rf}; NODE=${NODE:-r22}; MD=${MSKDIAM:-232}
R=$K_/rflex_$RT; ND=$R/node_$NODE; CK=$K_/cheapkids_${RT}_$NODE; MASK=$K_/focused_arm_${RT}_$NODE/arm_mask.mrc; O=$K_/focused_flex_${RT}_$NODE; mkdir -p $O && cd $O
log(){ echo "[$(date '+%F %T') focused-flex $RT $NODE] $*"; }; CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1); set -- $(awk -v u=$NODE '$1 == u' $R/queue.txt); PP=${3:-}; ROWS=${4:-}
MEAN=$(ls -td $ND/mean/*_reconstruct3D | head -1)/recvol_state01.mrc
[ -s $ND/n.simple ] && [ -s $MASK ] && [ -s $MEAN ] || { log "inputs missing"; exit 1; }
if [ -z "$(ls $O/*_flex_pca/flex_pca_coordinates.txt 2>/dev/null)" ]; then cp -f $ND/n.simple f.simple
  log "flex_pca on the node with the arm mask as the basis support (lp ${LPF:-15}, 8 states)"
  simple_exec prg=flex_pca projfile=f.simple vol1=$MEAN mskdiam=$MD lp=${LPF:-15} npreimages=8 rec_states=yes min_state_frac=0 umap=no n_probe_iters=${NPROBE:-2} pcg_mskfile=$MASK nparts=4 nthr=16 > flex.log 2>&1 || { log "flex_pca FAILED: $(grep -a -m2 'ERROR\|THROW' flex.log | cut -c1-200)"; exit 1; }; fi
FD=$(ls -td $O/*_flex_pca | head -1)
echo "== separability of the earlier closed states (f692r 16+15) in the ARM-FOCUSED latent (whole-molecule latent: LDA AUC 0.80, PC1 0.76)"
$PY $T/closed_latent.py --flex $FD/flex_pca_coordinates.txt --oris $CK/ptcl3D_parent.txt --rows $ROWS --ref $K_/f692r/coordtree_fam1/units.txt --ids 16,15 2>&1 | grep -v Warning | grep -v "  IC"
for ax in 1 w1; do kk=10; [ $ax = w1 ] && kk=12
  [ -s $O/children_$ax.txt ] || $PY $T/cheapkids_rf2.py --oris $CK/ptcl3D_parent.txt --rows $ROWS --flex $FD/flex_pca_coordinates.txt --out $O --k $kk --axis $ax --p 0.5 2>&1 | grep -v Warning | sed "s/^/focused latent /"
  NS=$(awk '$1 > 0' $O/children_$ax.txt | sort -n | uniq | wc -l)
  if [ -z "$(ls $O/rec_$ax/*_reconstruct3D/recvol_state$(printf %02d $NS).mrc 2>/dev/null)" ]; then mkdir -p $O/rec_$ax && cd $O/rec_$ax
    cp -f $PP r.simple && simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$O/children_$ax.txt mkdir=no > rsel.log 2>&1 \
      && simple_exec prg=update_project projfile=r.simple qsys_name=local nparts=4 nthr=16 mkdir=no > rup.log 2>&1 \
      && simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=$MD nstates=$NS nparts=4 nthr=16 > rec.log 2>&1 || { log "axis $ax reconstruct FAILED"; cd $O; continue; }
    RR=$(ls -td $O/rec_$ax/*_reconstruct3D | head -1); rm -f $RR/*_part*_even.mrc $RR/*_part*_odd.mrc $RR/rho_*; cd $O; fi
  RR=$(ls -td $O/rec_$ax/*_reconstruct3D | head -1); V=$K_/view_rflex/focused_flex_${RT}_${NODE}_$ax
  : > $O/spec_$ax.txt; for s in $(seq 1 $NS); do b=$RR/recvol_state$(printf %02d $s); echo "child_$(printf %02d $s) ${b}_even_unfil.mrc ${b}_odd_unfil.mrc" >> $O/spec_$ax.txt; done
  $PY $T/make_view.py $V $O/spec_$ax.txt > /dev/null 2>&1
  $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/arm_$ax.txt $(ls $V/*.mrc) 2>&1 | grep -v Warning | tail -1 | sed "s/^/axis $ax: /"
  awk 'NR == FNR {n[$1]++; next} !/^#/ {k = $1; sub(/.*child_/, "", k); sub(/\.mrc/, "", k); printf "   child %s (%d ptcls): closed-arm %.2f open-arm %.2f\n", k, n[k + 0], $2, $3}' $O/children_$ax.txt $O/arm_$ax.txt
done
log DONE
