#!/bin/bash
#SBATCH -J lsamp_ftest
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 01:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# the Fortran landscape_sample (build ~/SIMPLE_lsamp) on stree_port poles_r1_p1 / p2 node r, then lsamp_fcmp.py against the Python pole-pair
# coordinates (repro_r1/pq) and arm readouts of the held-out and final maps (padded to box 160).
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=16
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_lsamp/build; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; W=$K_/lsamp/repro_r1
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
for p in ${RUNS:-p1 p2}; do N=$K_/stree_port/poles_r1_$p/1_hetero_tree/node_r; O=$W/fortran/$p; mkdir -p $O; cd $O
  echo "[$(date '+%F %T')] $p"
  /usr/bin/time -v simple_exec prg=landscape_sample projfile=$N/stageB/b.simple dir=$N/stageB infile=$N/site/site_poles.txt \
      infile2=$N/split/split_children.txt mskdiam=232 nthr=16 mkdir=no > ls.log 2> time.log
  echo "exit $?"; grep "LANDSCAPE_SAMPLE\|resolvable\|CHOSEN\|steps CV\|check 3\|sample \|ERROR\|error" ls.log | cut -c1-220
  grep "Elapsed\|Maximum resident" time.log
  $PY $T/lsamp_fcmp.py $O $W/pq $p 160 1.464 2>&1 | grep -v Warning
  for d in $O/landscape_*/heldout_A_pad $O/landscape_*/heldout_B_pad $O/landscape_*/final_pad; do [ -d $d ] || continue
    $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $d/arm.txt $d 2>&1 | grep -v Warning > /dev/null
    echo "  ${d#$O/}: $(awk '!/^#/ {k = $1; sub(/.*state_/, "", k); sub(/\.mrc/, "", k); printf "%s %.2f/%.2f  ", k + 0, $2, $3}' $d/arm.txt)"; done
done
echo "[$(date '+%F %T')] DONE"
