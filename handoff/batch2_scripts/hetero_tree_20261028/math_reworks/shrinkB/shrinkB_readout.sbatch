#!/bin/bash
#SBATCH -J shrinkBr
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/shrinkB/logs/%x_%j.out
# All finished arms of every fixture side by side: stage-B labels and split children against the truth and each other
# (f2r_ab_readout.py), the per-call populations, and on r131 the closed-arm / open-arm occupancy of every stage-B map.
set -uo pipefail
M=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/shrinkB
L=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent
S=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/fixtures/scripts
R160=/mnt/beegfs/elmlund/afan/routes/toy55/maps160
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PPINFO=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2
while read -r KEY NODE TRUTH <&3; do
  [[ -z "$KEY" || "$KEY" == \#* ]] && continue
  TAG=$(basename $NODE | sed 's/^node_//'); arms=""; first=""
  for R in $M/runs/${KEY}_${TAG}_*; do
    m=${R##*_${TAG}_}; [ -s $R/stageB/b_ptcl3D.txt ] || continue
    [ -z "$first" ] && first=$R && arms="--arm production=$R/prod"
    arms="$arms --arm $m=$R/stageB"
  done
  [ -z "$arms" ] && continue
  echo "######## $KEY $TAG"
  for R in $M/runs/${KEY}_${TAG}_*; do echo "  ${R##*/} populations per call: $(for it in 1 2 3; do awk -v K=9 '/PCG DISTRIBUTED \| STATE=/ && /KIND=base/ { split($0, a, "|"); s = a[2]; gsub(/[^0-9]/, "", s); n = a[5]; gsub(/[^0-9]/, "", n); if (s + 0 > 0) t[s + 0] += n } END { for (k = 1; k <= 9; k++) if (t[k] > 0) printf "%d:%d ", k, t[k] }' $R/stageB/step$it.log 2> /dev/null; echo -n "| "; done)"; done
  $PY $L/f2r_ab_readout.py --node $NODE --truth $TRUTH $arms 2>&1
  if [[ $KEY == *r131* ]]; then
    O=$M/readout_${KEY}_${TAG}; mkdir -p $O; : > $O/spec.txt
    for e in $NODE/stageB/recvol_state??_even_unfil.mrc $M/runs/${KEY}_${TAG}_*/stageB/recvol_state??_even_unfil.mrc; do
      o=${e/_even_/_odd_}; k=$(basename $e | cut -c13-14)
      if [[ $e == $M/* ]]; then r=${e%/stageB/*}; tag=${r##*_${TAG}_}__s$k; else tag=production__s$k; fi
      [ -e $o ] && echo "$tag $e $o" >> $O/spec.txt
    done
    CL=$(ls $R160/docked/conf18_*_m160.mrc | head -1)
    $PY $S/make_view.py $O/views $O/spec.txt --box 160 > /dev/null 2>&1
    $PY $S/arm_readout.py $R160/toy_consensus.mrc $CL $O/arm.txt $O/views/*.mrc > /dev/null 2>&1
    echo "closed-arm / open-arm occupancy per stage-B map (closed corner >= 0.95 / <= 0.6):"; grep -v '^#' $O/arm.txt | awk '{n=split($1,p,"/"); print "   " p[n], $2, $3, $4}'
  fi
done 3< ${LIST:-$M/shrinkB_list.txt}
echo "done $(date)"
