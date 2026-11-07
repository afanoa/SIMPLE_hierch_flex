#!/bin/bash
#SBATCH -J klem_leaf_plainrec
#SBATCH -p csbdevel
#SBATCH -c 24
#SBATCH --mem=96G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Is a tree leaf's short mobile tip the particles or the focused refinement? The leaf's particles (state STATE of node NODE of rflex_<RT>,
# read-only) reconstructed with plain reconstruct3D (no pcg focus) at the node's refined poses -> lp 10 view -> closed-/open-arm readout and
# radial extent of the closed arm, next to the node's focused stage-B state map and the node's focus mask. Outputs in OUT only.
# env RT (f692st), NODE (r1), STATE (5), OUT (rnasep_klem/factorial/<RT>_<NODE>_s<STATE>_plain)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; RT=${RT:-f692st}; NODE=${NODE:-r1}; ST=${STATE:-5}
D=$K_/rflex_$RT/node_$NODE; O=${OUT:-$K_/factorial/${RT}_${NODE}_s${ST}_plain}; V=$K_/view_rflex/factorial_${RT}_${NODE}_s${ST}; mkdir -p $O $V; cd $O
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1); CONS=$T/maps160/toy_consensus.mrc; RD=$(ls -td $D/fstageB/*_refine3D | head -1)
awk -v s=$ST '{print ($1 == s) ? 1 : 0}' $D/refined_f.txt > $O/rows.txt; echo "leaf ${NODE}_s$ST: $(awk '$1 > 0' $O/rows.txt | wc -l) particles, poses $RD/b.simple"
if [ -z "$(ls $O/*_reconstruct3D/recvol_state01_even_unfil.mrc 2>/dev/null)" ]; then cp -f $RD/b.simple $O/r.simple
  simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$O/rows.txt mkdir=no > sel.log 2>&1 && simple_exec prg=update_project projfile=r.simple qsys_name=local nparts=3 nthr=8 mkdir=no > up.log 2>&1 \
    && simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=${MSKDIAM:-232} nstates=1 nparts=3 nthr=8 > rec.log 2>&1 || { echo "reconstruct FAILED: $(grep -a -m2 ERROR rec.log)"; exit 1; }; fi
RR=$(ls -td $O/*_reconstruct3D | head -1); b=$RD/recvol_state$(printf %02d $ST)
printf "plain %s %s\nfocused %s %s\n" $RR/recvol_state01_even_unfil.mrc $RR/recvol_state01_odd_unfil.mrc ${b}_even_unfil.mrc ${b}_odd_unfil.mrc > $O/spec.txt
$PY $T/make_view.py $V $O/spec.txt --box 160 > /dev/null 2>&1; $PY $T/arm_readout.py $CONS $CL $O/arm.txt $V/plain.mrc $V/focused.mrc 2>&1 | grep -v Warning
$PY $T/radial_extent.py $CONS $CL --maps plain_reconstruct3D:$V/plain.mrc focused_stageB:$V/focused.mrc --masks node_mask:$D/rsmask.mrc 2>&1 | grep -v Warning
