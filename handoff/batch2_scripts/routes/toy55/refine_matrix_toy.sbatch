#!/bin/bash
#SBATCH -J toyah_refine_matrix
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 08:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# refine3D_states settings matrix on the apo/holo toy's family 1 (87k particles, 10 k-means states, consensus/fast
# poses in): for each setting the wall time, the pose error vs the TRUE poses (pose_error.py), and the identification
# of the state maps against the 18 toy conformers. Baseline = the chain's existing run (local 20->8, maxits 6, 65 min).
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; B=$HOME/SIMPLE_kernel/snap/5e1e6aeee
export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel TOY_SMPD=1.464 OMP_NUM_THREADS=8
T=$A/rnasep_klem/toyah_tree; SRC=$T/src.simple; ST=$T/fam1/states.txt; TRUTH=$W/simah/all.txt; R=$A/rnasep_klem/toyah_refmat; mkdir -p $R; log(){ echo "[$(date '+%F %T') refmat] $*"; }
mkdir -p $W/maps160/docked_ah; for m in $W/maps160/docked/conf??_*_m160.mrc $W/maps160/holo_docked/conf1??_*_m160.mrc; do ln -sf $m $W/maps160/docked_ah/; done
score(){ d=$1; tag=$2; RD=$(ls -td $d/*_refine3D_states 2>/dev/null | head -1); [ -n "$RD" ] || { log "$tag: no output"; return; }
  simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $d/refined_ptcl3D.txt
  echo "== $tag: pose error vs truth"; $PY $W/pose_error.py $d/refined_ptcl3D.txt $TRUTH --states $ST --tag $tag 2>&1 | grep -v Warning | tail -3
  mkdir -p $d/states_dir; for w in $(seq -f %02g 1 10); do cp -f $RD/recvol_state$w.mrc $d/states_dir/state_$w.mrc 2>/dev/null; done
  echo "== $tag: identification vs the 18 toy conformers"; $PY $W/deposited_match.py $d/states_dir $W/maps160/docked_ah $W/maps160/toy_consensus.mrc 2>&1 | grep -v Warning | grep "identified"
  grep -a "POSTPROCESS: FSC" $RD/../refine.log 2>/dev/null | awk '{print $NF, $(NF-1)}' | sort -n | awk 'NR==1{b=$0} {l=$0} END{print "   FSC 0.143 range (best/worst state): " b " / " l}'; }
echo "== starting poses (fast/consensus) vs truth"; simple_exec prg=print_project_field projfile=$SRC oritype=ptcl3D mkdir=no 2>/dev/null > $R/start_ptcl3D.txt; $PY $W/pose_error.py $R/start_ptcl3D.txt $TRUTH --states $ST --tag start 2>&1 | grep -v Warning | tail -2
score $T/fam1/refine "baseline local 20->8 maxits6 (65 min)"
run(){ tag=$1; shift; d=$R/$tag; mkdir -p $d && cd $d && cp -f $SRC p.simple
  simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$ST mkdir=no > sel.log 2>&1 || { log "$tag selection FAILED"; return; }
  t0=$(date +%s); log "$tag: refine3D_states $*"; simple_exec prg=refine3D_states projfile=p.simple nstates=10 mskdiam=220 pgrp=c1 nparts=10 nthr=8 "$@" > refine.log 2>&1 || { log "$tag FAILED: $(grep -a -m2 ERROR refine.log | cut -c1-200)"; return; }
  log "$tag: $(( ($(date +%s) - t0) / 60 )) min, $(grep -a -c ITERATION refine.log) iterations"; score $d "$tag"; cd $R; }
run nomarch  pose_policy=local lpstart=10 lpstop=8 maxits=5
run band8    pose_policy=local lpstart=8  lpstop=8 maxits=5
run fixed    pose_policy=fixed lpstart=20 lpstop=8 maxits=5
run coarse   pose_policy=local lpstart=10 lpstop=8 maxits=5 nspace=1000 nspace_sub=100
run small    pose_policy=local lpstart=10 lpstop=8 maxits=5 nsample=20000
log "ALL DONE"
