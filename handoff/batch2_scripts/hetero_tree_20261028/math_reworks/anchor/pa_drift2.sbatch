#!/bin/bash
#SBATCH -J pa_drift2
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/anchor/logs/%x_%j.out
# pose drift by depth in two older trees that kept their node projects (qtoy2_off to depth 4, the RNase P subtree r1211 to depth 3)
set -uo pipefail
A=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/anchor
AD=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/ancestral_dsnr/gates
B=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/merge_master2/install
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$B/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 SIMPLE_PATH=$B
export OMP_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1 PPINFO=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo
D=$A/drift2; mkdir -p $D
dump() { [ -s $2 ] || simple_exec prg=print_project_field projfile=$1 oritype=ptcl3D mkdir=no > $2 2> $2.err; }
declare -A TRUTH=([qtoy2_off]=qtoy:/mnt/beegfs/elmlund/afan/routes/toy55/simqtoy2/truth_qtoy.txt \
  [r1211_s1_off]=closed:$AD/r1211_s1_off/fixture_rows.txt:/mnt/beegfs/elmlund/afan/hetero_tree_regression_20261026/runs/real_root_split_contrast_m8c52_s1/closed_reference_rows.txt)
for t in qtoy2_off r1211_s1_off; do
  T=$AD/$t/1_hetero_tree
  for nd in $(ls -d $T/node_r? $T/node_r?? $T/node_r??? 2>/dev/null); do
    n=$(basename $nd); id=${n#node_}; pid=${id%?}; P=$T/node_$pid
    for f in n.simple stageA/r.simple stageB/b.simple; do [ -s $nd/$f ] || continue 2; done
    [ -s $P/stageB/b.simple ] || continue
    k=${t}_$id
    dump $P/stageB/b.simple $D/${t}_${pid}_b.txt
    dump $nd/n.simple $D/${k}_n.txt; dump $nd/stageA/r.simple $D/${k}_a.txt; dump $nd/stageB/b.simple $D/${k}_b.txt
    $PY - $P/rows.txt $nd/rows.txt $D/${k}_prow.txt <<'PYE'
import sys, numpy as np
p = np.loadtxt(sys.argv[1], dtype=int, ndmin=1); c = np.loadtxt(sys.argv[2], dtype=int, ndmin=1)
pos = np.searchsorted(p, c)
assert (p[np.minimum(pos, len(p) - 1)] == c).all(), 'child rows not in parent rows'
np.savetxt(sys.argv[3], pos + 1, fmt='%d')
PYE
    echo "######## $k (depth ${#id} - 1)"
    $PY $A/pa_drift.py $k $D/${t}_${pid}_b.txt $D/${k}_n.txt $D/${k}_a.txt $D/${k}_b.txt 1.464 $D/${k}_prow.txt 2>&1 | grep -vE "^in->A|^A->B" | grep -vE "^       (geodesic|direction|in-plane|shift)" | head -0
    $PY $A/pa_drift.py $k $D/${t}_${pid}_b.txt $D/${k}_n.txt $D/${k}_a.txt $D/${k}_b.txt 1.464 $D/${k}_prow.txt 2>&1 | sed -n '1,4p;/^in->B/,/shift/p'
    $PY $A/pa_tail.py $k $nd ${TRUTH[$t]} $D/${k}_n.txt $D/${k}_b.txt 2>&1
  done
done
echo done $(date)
