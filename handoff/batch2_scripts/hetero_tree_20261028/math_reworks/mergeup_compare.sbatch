#!/bin/bash
#SBATCH -J mrw_mucmp
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/logs/%x_%j.out
# the merge-up of arm B against arm A on the gate trees: leaf tables, merged units (leaf sets), the equivalence decisions and
# the merge-up step's wall (A=bc B=mu by default); writes gates/<B>_mergeup_compare.txt
set -uo pipefail
A=${A:-bc}; B=${B:-mu}; M=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks
OUT=$M/gates/${B}_mergeup_compare.txt; : > $OUT
wall() { local d=$1; local s e; s=$(stat -c %Y $d/hetero_step.sh 2>/dev/null); e=$(stat -c %Y $d/hetero_step.exit 2>/dev/null); [ -n "$s" ] && [ -n "$e" ] && echo "$(( (e - s) / 60 )) min" || echo "n/a"; }
for ds in qtoynull olenull micnull qtoy2 r131 r131_s1 ole; do
  TA=$M/gates/${A}_$ds/1_hetero_tree; TB=$M/gates/${B}_$ds/1_hetero_tree
  echo "## $ds" >> $OUT
  test -e $TB/hetero_tree_leaves.txt || { echo "   $B: no leaves yet" >> $OUT; continue; }
  na=$(grep -c . $TA/hetero_tree_leaves.txt); nb=$(grep -c . $TB/hetero_tree_leaves.txt)
  if diff -q <(awk '{print $1, $2, $3}' $TA/hetero_tree_leaves.txt | sort) <(awk '{print $1, $2, $3}' $TB/hetero_tree_leaves.txt | sort) > /dev/null; then
    echo "   leaves: $A $na = $B $nb (identical names and sizes)" >> $OUT
  else
    echo "   leaves: $A $na vs $B $nb (differ):" >> $OUT
    diff <(awk '{print $1, $2, $3}' $TA/hetero_tree_leaves.txt | sort) <(awk '{print $1, $2, $3}' $TB/hetero_tree_leaves.txt | sort) | sed 's/^/      /' >> $OUT
  fi
  for arm in $A $B; do
    T=$M/gates/${arm}_$ds/1_hetero_tree
    if [ -e $T/mergeup/merged_leaves.txt ]; then
      echo "   $arm merge-up: $(grep -vc '^#' $T/mergeup/merged_leaves.txt) units from $(grep -c . $T/hetero_tree_leaves.txt) leaves, wall $(wall $T/mergeup); merged sets: $(grep -v '^#' $T/mergeup/merged_leaves.txt | awk '$2 ~ /,/ {printf "%s ", $2}')" >> $OUT
      grep -E "merge rule inputs" $T/mergeup/hetero_step.log | cut -c1-200 | sed 's/^/      /' >> $OUT
      grep -E "^equivalence:|weakest kept pair|re-test: weakest|^>>>   level [0-9]+ node [0-9]+: [0-9]+ units ->" $T/mergeup/hetero_step.log | cut -c1-200 | sed 's/^/      /' >> $OUT
    else
      echo "   $arm merge-up: none (one leaf or not run)" >> $OUT
    fi
  done
done
echo "## core-hours of the $B arm (sacct by WorkDir, last 24 h)" >> $OUT
sacct -u atifao -S $(date -d '-24 hours' +%Y-%m-%dT%H:%M) -X -n -o WorkDir%120,CPUTimeRAW 2>/dev/null | awk -v a=$B '$1 ~ "math_reworks/gates/"a"_" {split($1,p,"/"); for(i=1;i<=length(p);i++) if (p[i] ~ "^"a"_") ds=p[i]; c[ds]+=$2} END {for (k in c) printf "   %s %.1f core-h\n", k, c[k]/3600}' >> $OUT
echo done >> $OUT
