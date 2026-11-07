#!/bin/bash
#SBATCH -J mrw_readout
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/logs/%x_%j.out
# gate readouts of one arm (ARM=base|sub|r2): null largest-leaf shares, qtoy2 leaves vs truth, r131 closed-corner readout
set -uo pipefail
: ${ARM:?}
M=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks; S=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/fixtures/scripts
PY=/mnt/nasapps/production/python/3.12/bin/python3; export PYTHONDONTWRITEBYTECODE=1
OUT=$M/gates/${ARM}_readout.txt; : > $OUT
for ds in qtoynull olenull micnull qtoy2 r131 ole; do
  T=$M/gates/${ARM}_$ds/1_hetero_tree
  test -e $T/hetero_tree_leaves.txt || { echo "## $ds: no leaves yet" >> $OUT; continue; }
  echo "## $ds ($(grep -c . $T/hetero_tree_leaves.txt) leaves; $(grep -h finished $M/gates/${ARM}_$ds/run_manifest.txt | head -1))" >> $OUT
  awk '$2 != "junk" {n[NR]=$3; s+=$3; if ($3>m) m=$3} END {printf "   largest leaf %d of %d = %.3f\n", m, s, (s>0?m/s:0)}' $T/hetero_tree_leaves.txt >> $OUT
  grep -E "HETERO_TREE node [a-z0-9]*: (split group|states)" $M/gates/${ARM}_$ds/tree.log | cut -c1-140 >> $OUT
done
echo "## qtoy2 leaves vs truth" >> $OUT
test -e $M/gates/${ARM}_qtoy2/1_hetero_tree/hetero_tree_labels.txt && $PY -I $S/qtoy_leaves.py $ARM=$M/gates/${ARM}_qtoy2/1_hetero_tree >> $OUT 2>&1
echo "## r131 closed-corner readout (stage-B maps, closed corner >= 0.95 / <= 0.6)" >> $OUT
test -e $M/gates/${ARM}_r131/1_hetero_tree/hetero_tree_nodes.txt && $PY -I $S/closed_nodes2.py $M/gates/${ARM}_r131_readout $ARM=$M/gates/${ARM}_r131/1_hetero_tree >> $OUT 2>&1
echo "## core-hours (sacct, this arm's controller jobs and their children by WorkDir)" >> $OUT
for ds in qtoynull olenull micnull qtoy2 r131 ole; do
  j=$(grep -h '^job=' $M/gates/${ARM}_$ds/run_manifest.txt 2>/dev/null | cut -d= -f2); test -n "$j" || continue
  sacct -j $j -X -n -o Elapsed,CPUTimeRAW 2>/dev/null | awk -v d=$ds '{printf "   %s controller elapsed %s\n", d, $1}' >> $OUT
done
sacct -u atifao -S $(date -d '-12 hours' +%Y-%m-%dT%H:%M) -X -n -o JobID,JobName%20,WorkDir%120,CPUTimeRAW 2>/dev/null | awk -v a=$ARM '$3 ~ "math_reworks/gates/"a"_" {split($3,p,"/"); for(i=1;i<=length(p);i++) if (p[i] ~ "^"a"_") ds=p[i]; c[ds]+=$4} END {for (k in c) printf "   %s %.1f core-h\n", k, c[k]/3600}' >> $OUT
echo done >> $OUT
