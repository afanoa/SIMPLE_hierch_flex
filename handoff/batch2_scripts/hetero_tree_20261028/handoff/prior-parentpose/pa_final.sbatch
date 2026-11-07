#!/bin/bash
#SBATCH -J pa_final
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:40:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/anchor/final_readout.txt
# every readout of the anchor test in one file: labels and agreement (pa_agree), qtoy2 purity (pa_qtoy), the qtoyflex arm cut
# (pa_qflex), drift and tail per replay (pa_post body)
A=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/anchor
G=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/gates
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PYTHONDONTWRITEBYTECODE=1 PPINFO=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo
echo "==================== labels, truth readouts, pairwise agreement"
awk '!s[$1]++ {print $1, $2, $5}' $A/pa_list1.txt | while read -r KEY NODE TRUTH; do $PY $A/pa_agree.py $KEY $NODE $TRUTH; done
echo "==================== qtoy2 children purity"
for n in r1 r2 r3 r4 r5; do $PY $A/pa_qtoy.py mm2_qtoy2_$n $G/mm2_qtoy2/1_hetero_tree/node_$n; done
echo "==================== qtoyflex r1 arm cut"
$PY $A/pa_qflex.py mm2_qtoyflex_r1 $G/mm2_qtoyflex/1_hetero_tree/node_r1
echo "==================== r131 stage-B map readout (closed-arm / open-arm)"
for r in $A/replays/mm2_r131_*; do [ -s $r/readout/arm.txt ] && grep -v '^#' $r/readout/arm.txt | awk -v n=$(basename $r) '{k=split($1,p,"/"); print "  " n, p[k], $2, $3}'; done
echo "==================== drift (in -> B) and tail per replay"
while read -r KEY NODE ARM SEED TRUTH XA XB; do
  R=$A/replays/${KEY}_${ARM}_s$SEED
  [ -s $R/stageB/b_ptcl3D.txt ] || continue
  echo "######## ${KEY}_${ARM}_s$SEED"
  $PY $A/pa_drift.py ${KEY}_${ARM}_s$SEED $A/drift/${KEY%_r?}_r_b.txt $A/drift/${KEY}_n.txt $R/stageA/a_ptcl3D.txt $R/stageB/b_ptcl3D.txt 1.464 $NODE/rows.txt 2>&1 | sed -n '/^in->B/,/shift/p'
  $PY $A/pa_tail.py ${KEY}_${ARM}_s$SEED $NODE $TRUTH $A/drift/${KEY}_n.txt $R/stageB/b_ptcl3D.txt 2>&1 | sed -n '1p'
done < $A/pa_list1.txt
echo "done $(date)"
