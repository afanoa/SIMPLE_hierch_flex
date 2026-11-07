#!/bin/bash
# rescoring lines for finished runs: the run's own configuration (its list line minus --labels) with --labels file:<its posteriors>,
# so the hold-out curve and the pole rule's seed test are recomputed for them. usage: make_rescore.sh RUN [RUN ...] > list_rescore.txt
P=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo
for run in "$@"; do
  line=$(grep -h "^$run " $P/list_phase*.txt | head -1); [ -n "$line" ] || { echo "no list line for $run" >&2; continue; }
  posts=$(ls $P/runs/$run/site_post_*.txt 2>/dev/null | tr '\n' ',' | sed 's/,$//'); [ -n "$posts" ] || { echo "no posteriors in $run" >&2; continue; }
  extra=$(echo "$line" | awk '{$1=$2=$3=$4=$5=$6=""; print $0}' | sed -E 's/--labels [^ ]+//; s/--maxits [0-9]+//; s/  +/ /g')
  echo "rs_$run $(echo "$line" | awk '{print $2, $3, $4, $5, $6}')$extra --labels file:$posts"
done
