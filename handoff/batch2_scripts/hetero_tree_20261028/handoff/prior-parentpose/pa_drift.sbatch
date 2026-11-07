#!/bin/bash
#SBATCH -J pa_drift
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 00:40:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/anchor/logs/%x_%j.out
# production pose drift of the depth-1 gate nodes (incoming -> stage A -> stage B) against the parent's last-iteration spread
set -uo pipefail
A=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/anchor
G=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/gates
B=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/merge_master2/install
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$B/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 SIMPLE_PATH=$B
export OMP_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1
D=$A/drift; mkdir -p $D
dump() { [ -s $2 ] || simple_exec prg=print_project_field projfile=$1 oritype=ptcl3D mkdir=no > $2 2> $2.err; }
for t in mm2_r131 mm2_qtoy2 mm2_qtoyflex; do
  T=$G/$t/1_hetero_tree
  dump $T/node_r/stageB/b.simple $D/${t}_r_b.txt
  for n in $(ls -d $T/node_r? | xargs -n1 basename); do
    k=${t}_${n#node_}
    dump $T/$n/n.simple $D/${k}_n.txt
    dump $T/$n/stageA/r.simple $D/${k}_a.txt
    dump $T/$n/stageB/b.simple $D/${k}_b.txt
    $PY $A/pa_drift.py $k $D/${t}_r_b.txt $D/${k}_n.txt $D/${k}_a.txt $D/${k}_b.txt 1.464 $T/$n/rows.txt > $D/${k}_drift.txt 2>&1
    cat $D/${k}_drift.txt
  done
done
echo done $(date)
