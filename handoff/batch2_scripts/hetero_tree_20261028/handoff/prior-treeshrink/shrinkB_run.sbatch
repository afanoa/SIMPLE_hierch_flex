#!/bin/bash
#SBATCH -J shrinkB
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/shrinkB/logs/%x_%j.out
# Stage B as three single-iteration refine3D calls in one directory (the node's own stageB call with maxits=1, startit=1..3),
# the state references between calls passed through shrinkB.py MODE (copy = the control, dev1 / dev2 = deviation shrinkage),
# then the node's split step under the gate's split environment; readouts against production (and the control when it exists).
# LIST lines: "KEY NODE TRUTH"; MODE and LIST from the environment.
set -uo pipefail
: ${MODE:?} ${LIST:?}
M=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/shrinkB
L=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent
B=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/merge_master2/install
S=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/fixtures/scripts
R160=/mnt/beegfs/elmlund/afan/routes/toy55/maps160
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$B/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 SIMPLE_PATH=$B
export OMP_NUM_THREADS=8 PPINFO=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo PYTHONDONTWRITEBYTECODE=1
export SIMPLE_SPLIT_CONTRAST=maps SIMPLE_TREE_KEEPMAPS=yes
counts() { # LOG K: even + odd base-pair populations per state of the iteration's PCG assembly
  awk -v K=$2 '/PCG DISTRIBUTED \| STATE=/ && /KIND=base/ { split($0, a, "|"); s = a[2]; gsub(/[^0-9]/, "", s); n = a[5]; gsub(/[^0-9]/, "", n)
    if (s + 0 > 0) { if (a[3] ~ /even/) e[s + 0] = n; else o[s + 0] = n } }
    END { for (k = 1; k <= K; k++) printf "%s%d", (k > 1 ? "," : ""), e[k] + o[k] }' $1
}
while read -r KEY NODE TRUTH <&3; do
  [[ -z "$KEY" || "$KEY" == \#* ]] && continue
  TAG=$(basename $NODE | sed 's/^node_//')
  R=$M/runs/${KEY}_${TAG}_$MODE
  [ -e $R ] && { echo "EXISTS $R"; continue; }
  mkdir -p $R/stageB $R/split $R/prod $R/readout
  SEED=$(grep -o 'SIMPLE_SEED=[0-9]*' $NODE/stageB/hetero_step.sh | cut -d= -f2); export SIMPLE_SEED=$SEED
  base=$(grep -E 'simple_exec prg=refine3D' $NODE/stageB/hetero_step.sh | sed -e "s#^[^ ]*/simple_exec#$B/bin/simple_exec#" -e 's# > hetero_step.log.*##' -e 's#maxits=[0-9.]*#maxits=1#')
  K=$(echo "$base" | grep -o 'nstates=[0-9]*' | cut -d= -f2); MSK=$(echo "$base" | grep -o 'mskdiam=[0-9.]*' | cut -d= -f2)
  echo "== $KEY $TAG $MODE: K $K seed $SEED $(date +%T)"
  cp $NODE/recB/b.simple $R/stageB/b.simple
  cd $R/stageB; ok=1
  for it in 1 2 3; do
    cmd="$base startit=$it"
    if [ $it -gt 1 ]; then
      dst=$R/stageB/refs_it$((it - 1)); cnt=$(counts step$((it - 1)).log $K)
      $PY $M/shrinkB.py $R/stageB $dst $MODE --counts $cnt --mskdiam $MSK > shrink_it$((it - 1)).log 2>&1 || { echo "   shrink failed"; ok=0; break; }
      for s in $(seq 1 $K); do cmd=$(echo "$cmd" | sed "s#vol$s=[^ ]*#vol$s=$dst/recvol_state$(printf %02d $s).mrc#"); done
    fi
    echo "$cmd" > cmd$it.txt
    eval "$cmd" > step$it.log 2>&1; rc=$?
    echo "   call $it exit $rc $(date +%T) counts $(counts step$it.log $K)"
    [ $rc = 0 ] && grep -q "SIMPLE_REFINE3D NORMAL STOP" step$it.log || { echo "   call $it failed"; ok=0; break; }
  done
  [ $ok = 1 ] || continue
  simple_exec prg=print_project_field projfile=$R/stageB/b.simple oritype=ptcl3D mkdir=no > $R/stageB/b_ptcl3D.txt 2> /dev/null
  cp $R/stageB/b.simple $R/split/s.simple; cd $R/split
  scmd=$(grep -E 'simple_exec prg=state_split' $NODE/split/hetero_step.sh | sed -e "s#^[^ ]*/simple_exec#$B/bin/simple_exec#" -e 's# > hetero_step.log.*##')
  echo "$scmd" > cmd.txt; eval "$scmd" > step.log 2>&1; echo "   split exit $? $(date +%T)"
  ln -sfn $R/split/split_children.txt $R/stageB/split_children.txt
  simple_exec prg=print_project_field projfile=$NODE/stageB/b.simple oritype=ptcl3D mkdir=no > $R/prod/b_ptcl3D.txt 2> /dev/null
  ln -sfn $NODE/split/split_children.txt $R/prod/split_children.txt
  arms="--arm production=$R/prod"
  C=$M/runs/${KEY}_${TAG}_copy/stageB
  [ $MODE != copy ] && [ -s $C/b_ptcl3D.txt ] && arms="$arms --arm copy=$C"
  $PY $L/f2r_ab_readout.py --node $NODE --truth $TRUTH $arms --arm $MODE=$R/stageB 2>&1 | tee $R/readout/labels.txt
  if [[ $KEY == *r131* ]]; then
    : > $R/readout/spec.txt
    for e in $R/stageB/recvol_state??_even_unfil.mrc $NODE/stageB/recvol_state??_even_unfil.mrc; do
      o=${e/_even_/_odd_}; k=$(basename $e | cut -c13-14); tag=$([[ $e == $R/* ]] && echo $MODE || echo production)__s$k
      [ -e $o ] && echo "$tag $e $o" >> $R/readout/spec.txt
    done
    CL=$(ls $R160/docked/conf18_*_m160.mrc | head -1)
    $PY $S/make_view.py $R/readout/views $R/readout/spec.txt --box 160 > /dev/null 2>&1
    $PY $S/arm_readout.py $R160/toy_consensus.mrc $CL $R/readout/arm.txt $R/readout/views/*.mrc > /dev/null 2>&1
    echo "closed-arm / open-arm occupancy per stage-B map (closed corner >= 0.95 / <= 0.6):"; grep -v '^#' $R/readout/arm.txt | awk '{n=split($1,p,"/"); print "   " p[n], $2, $3, $4}'
  fi
done 3< $LIST
echo "done $(date)"
