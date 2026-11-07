#!/bin/bash
#SBATCH -J land_leafext
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 01:00:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_leafext.%j.out
# extremes of four independent r1 runs: pooled stable leaf directions, each run's two extreme leaves per axis, compared map against
# map against half-map noise (node mask lp 15 / lp 10, focus region lp 15); lp-10 views + arm readouts of the extreme leaves; and the
# shared input poses of node_r1/n.simple dumped for the path side readout
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 SIMPLE_EMAIL=nil
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; L=$K_/landscape
O=$L/leafext; V=$K_/view_rflex/leafext_r1; mkdir -p $O $V/halves && cd $O; log(){ echo "[$(date '+%F %T') leafext] $*"; }
RUNS="script=$K_/rflex_f692st/node_r1/fstageB/2_refine3D:$K_/repro_r1/labels/S_B.txt port1=$K_/stree_port/tree_f692r1/1_hetero_tree/node_r/stageB:$K_/repro_r1/labels/P_B.txt"
RUNS="$RUNS port2=$K_/stree_port/tree_f692r1_seed2/1_hetero_tree/node_r/stageB:$K_/repro_r1/dscore_port2/labels.txt port3=$K_/stree_port/tree_f692r1_rep3/1_hetero_tree/node_r/stageB:$K_/repro_r1/dscore_port3/labels.txt"
[ -s $O/poses_n.txt ] || simple_exec prg=print_project_field projfile=$K_/rflex_f692st/node_r1/n.simple oritype=ptcl3D mkdir=no > $O/poses_n.txt 2> /dev/null
log "shared input poses (node_r1/n.simple): $(grep -c . $O/poses_n.txt) rows"
log "node mask, lp 15"; $PY -u $T/land_leafext.py $RUNS --lp 15 --naxes 2 --spec $O/spec.txt
log "node mask, lp 10"; $PY -u $T/land_leafext.py $RUNS --lp 10 --naxes 2
log "focus region (rsmask_variance clipped to 0-1), lp 15"; $PY -u $T/land_leafext.py $RUNS --lp 15 --naxes 2 --mask $K_/rflex_f692st/node_r1/rsmask_variance.mrc
$PY $T/make_view.py $V $O/spec.txt --lp 10 --box 160 > /dev/null 2>&1
while read nm e o; do cp -f $e $V/halves/${nm}_even.mrc; cp -f $o $V/halves/${nm}_odd.mrc; done < $O/spec.txt
log "views: $(ls $V/state_*.mrc | wc -l) extreme-leaf maps lp 10 unsharpened in $V, halves in $V/halves"
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
$PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/arm.txt $V 2>&1 | grep --line-buffered -v "Warning\|lp = lambda" | sed 's#^#  #'
log "DONE: $O"
