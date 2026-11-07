#!/bin/bash
#SBATCH -J s30split2_truth
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Truth scorecards of the split session's s30split2 nodes r1-r5 (read-only): their flex states, focused-refined states and my cheap children.
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; GT=$T/simswing30/gt_particles.txt
for u in r1 r2 r3 r4 r5; do [ $u = r1 ] && { D=$K_/rflex_s30split2/node_r1; C=$K_/cheapkids_s30split2_r1; } || { D=$K_/rflex_s30split2_$u/node_$u; C=$K_/cheapkids_s30split2_${u}_$u; }
  echo "===== s30split2 $u"; for f in $D/flex_states.txt:flex $D/refined_states.txt:refined $C/children.txt:cheap_children; do [ -s ${f%%:*} ] && $PY $T/scorecard.py ${f%%:*} $GT --codes swing:200:1 --tag "${f#*:}" 2>&1 | grep -v Warning; done; done
