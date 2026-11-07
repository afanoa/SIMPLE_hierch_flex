#!/bin/bash
#SBATCH -J klem_phaseB
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 00:22:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Phase B: 5-bin map series at the CAHRA M s1 root (open) and the cdisc M root (A, closed) along the best population-wide score
# (pick_score.py), with the tree's own reconstruct3D line (nstates=5, postprocess=no); then the step-profile and mixture tests
# (phaseB_maptests.py) on these and the existing IgM series
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_stree/snap/m8c52t; NP=2; NT=4
export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local SIMPLE_EMAIL=nil OMP_NUM_THREADS=$NT SIMPLE_SEED=1
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; S=$K/stree_port
OUT=$K/factorial/phaseB; mkdir -p $OUT; log(){ echo "[$(date "+%F %T") phaseB] $*"; }
CA=$S/arm_M_c11_s1/1_hetero_tree; CD=$S/fulltoy_m52M_cd1/1_hetero_tree; GT=/mnt/beegfs/elmlund/EMPIAR_downloads/13991/data/source.star; ST=$K/factorial/junk_test/simsites_truth_jt.txt
for e in "cahra_open:$CA:open:cahra:$GT" "cdisc_A:$CD:A:cdisc:$ST" "cdisc_closed:$CD:C:cdisc:$ST"; do IFS=: read nm D tgt kind tf <<< "$e"
  W=$OUT/$nm; [ -e $W ] && mv $W ${W}_old_$(date +%s); mkdir -p $W/halves $W/work; cd $W/work
  $PY $T/pick_score.py score.txt $tgt $kind $tf $D > ../pick.txt 2>&1 || { log "$nm pick FAILED: $(tail -2 ../pick.txt)"; continue; }
  $PY $T/motion_bins3.py score.txt states.txt $kind $tf > ../bins.txt 2>&1 || { log "$nm bins FAILED: $(tail -1 ../bins.txt)"; continue; }
  cp $D/node_r/stageB/b.simple p.simple; simple_exec prg=update_project projfile=p.simple qsys_name=local mkdir=no > up.log 2>&1
  simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=states.txt mkdir=no > sel.log 2>&1 || { log "$nm selection FAILED"; continue; }
  REC=$(grep -a -o "prg=reconstruct3D[^>]*" $D/node_r/recB0/hetero_step.sh | sed -e 's/projfile=[^ ]*/projfile=p.simple/' -e 's/nstates=[^ ]*/nstates=5/' -e "s/nparts=[^ ]*/nparts=$NP/" -e "s/nthr=[^ ]*/nthr=$NT/")
  REC="$REC postprocess=no"; t0=$(date +%s); simple_exec $REC > rec.log 2>&1 || { log "$nm reconstruct3D FAILED: $(grep -a -m2 ERROR rec.log)"; continue; }
  for k in 1 2 3 4 5; do s=$(printf %02d $k); cp recvol_state$s.mrc ../bin$k.mrc; cp recvol_state${s}_even.mrc ../halves/bin${k}_even.mrc; cp recvol_state${s}_odd.mrc ../halves/bin${k}_odd.mrc; done
  log "$nm reconstruct3D $(( $(date +%s) - t0 )) s: $REC"; { echo "######## $nm"; cat ../pick.txt ../bins.txt; } | tee -a $OUT/summary.txt
done
I=$K/factorial/motion_series_igm; MI=$S/fulltoy_m52M_igm/1_hetero_tree/node_r/mask.mrc; MB=$S/fulltoy_igm/1_hetero_tree/node_r/mask.mrc
$PY $T/phaseB_maptests.py $OUT/cahra_open:$CA/node_r/mask.mrc:160 $OUT/cdisc_A:$CD/node_r/mask.mrc:225.456 $OUT/cdisc_closed:$CD/node_r/mask.mrc:225.456 \
  $I/m52M_r_s1_a4th_q5:$MI:220 $I/m52M_r2_s4_a2th_q5:$MI:220 $I/m8c17_r2_s1_a2th_q5:$MB:220 2>&1 | grep -v Warning | tee $OUT/maptests.txt
log DONE
