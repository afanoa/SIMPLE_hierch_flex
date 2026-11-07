#!/bin/bash
#SBATCH -J klem_view_check
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 00:40:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Do the bins of each 5-bin series differ in views or defocus? print_project_field on each series' work/p.simple (its states are the
# bins) piped into view_check.py, 8 at a time; VIEW lines in factorial/phaseE/view_check.txt
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_stree/snap/m8c52t; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local SIMPLE_EMAIL=nil OMP_NUM_THREADS=1
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; E=$K/factorial/phaseE; L=()
for d in $E/*/*/ $E/rnasep_m33/node_*/*/ $K/factorial/phaseB/*/ $K/factorial/motion_series_igm/*_q5/ /mnt/beegfs/elmlund/afan/ole_12707/motion_series/*/; do
  b=$(basename $d); [ -s $d/work/p.simple ] && [ -s $d/bin5.mrc ] || continue; case $d in *_old_*) continue;; esac
  case $d in */phaseB/*|*/motion_series_igm/*) b=ex_$b;; */ole_12707/*) b=ole_$b;; esac
  case $b in *cahra*) L+=("$b:$d:11");; *) L+=("$b:$d:1");; esac
done
echo "[$(date "+%F %T") view_check] ${#L[@]} series"; mkdir -p $E/view_parts
for p in 0 1 2 3 4 5 6 7; do
  ( for ((i = p; i < ${#L[@]}; i += 8)); do IFS=: read nm d ns <<< "${L[$i]}"; cd $d/work
      simple_exec prg=print_project_field projfile=p.simple oritype=ptcl3D 2>/dev/null | $PY $T/view_check.py $nm $ns; done > $E/view_parts/part$p.txt 2>&1 ) &
done; wait
cat $E/view_parts/part*.txt | grep "^VIEW" | sort > $E/view_check.txt; echo "[$(date "+%F %T") view_check] DONE: $(wc -l < $E/view_check.txt) VIEW lines"
