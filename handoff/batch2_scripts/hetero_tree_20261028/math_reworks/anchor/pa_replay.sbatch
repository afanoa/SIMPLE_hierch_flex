#!/bin/bash
#SBATCH -J pa_replay
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/anchor/logs/%x_%A_%a.out
# pose-anchor test: a depth-1 gate node replayed from its production seeds with its own saved stage A -> recB -> stage B ->
# split calls (the f2r replay harness without the label swap), the arm's keys appended to the stage A and stage B refine3D
# calls. Readouts: stage-B labels and split children against the truth and production (f2r_ab_readout.py), the pose drift
# from the incoming poses (pa_drift.py), on r131 the closed-arm / open-arm readout of the stage-B maps.
# LIST lines: "KEY NODE ARM SEED TRUTH XA XB" (XA/XB: stage A/B extra keys joined by '+', or '-')
set -uo pipefail
: ${LIST:?}
A=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/anchor
L=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent
B=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/merge_master2/install
BSM=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/merge_sm/install
S=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/fixtures/scripts
R160=/mnt/beegfs/elmlund/afan/routes/toy55/maps160
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$B/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 SIMPLE_PATH=$B
export OMP_NUM_THREADS=8 PPINFO=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo PYTHONDONTWRITEBYTECODE=1
export SIMPLE_SPLIT_CONTRAST=maps SIMPLE_TREE_KEEPMAPS=yes
read -r KEY NODE ARM SEED TRUTH XA XB < <(sed -n "${SLURM_ARRAY_TASK_ID:-1}p" $LIST)
export SIMPLE_SEED=$SEED
R=$A/replays/${KEY}_${ARM}_s$SEED
[ -e $R ] && { echo "EXISTS $R"; exit 3; }
XA=$([ "$XA" = "-" ] && echo "" || echo " ${XA//+/ }"); XB=$([ "$XB" = "-" ] && echo "" || echo " ${XB//+/ }")
LABELS=$NODE/seeds/seeds.txt
K=$(grep -v '^0$' $LABELS | sort -n | uniq | wc -l)
mkdir -p $R/stageA $R/recB $R/stageB $R/split $R/prod $R/readout
echo "replay $KEY $ARM seed $SEED: K $K; stage A +[$XA] stage B +[$XB]; $(date)"
inject() { ( cd $1 && $BSM/bin/simple_exec prg=write_states projfile=$2 infile=$3 mkdir=no > $2.sel.log 2>&1 ) || echo "write_states failed in $1"; }
step() { # DIR SAVED SEDX: the saved step's simple_exec call with the binary and SEDX swaps, run in $R/DIR
  local d=$1 src=$2 sedx=$3 cmd
  cd $R/$d
  cmd=$(grep -E 'simple_exec prg=' $NODE/$src/hetero_step.sh | sed -e "s#^[^ ]*/simple_exec#$B/bin/simple_exec#" -e 's# > hetero_step.log.*##' -e "$sedx")
  echo "$cmd" > cmd.txt
  echo "== $d $(date +%T)"; eval "$cmd" > step.log 2>&1; echo "   exit $? $(date +%T)"
}
cp $NODE/n.simple $R/stageA/r.simple; inject $R/stageA r.simple $LABELS
step stageA stageA "s#\$#$XA#"
simple_exec prg=print_project_field projfile=$R/stageA/r.simple oritype=ptcl3D mkdir=no > $R/stageA/a_ptcl3D.txt 2> /dev/null
KB=$($PY $L/f2r_fold_states.py $R/stageA/a_ptcl3D.txt $R/recB/labels_B.txt)
echo "stage A states fold to $KB"
drop=""; for n in $(seq $((KB + 1)) 9); do drop="$drop; s# vol$n=[^ ]*##"; done
cp $R/stageA/r.simple $R/recB/b.simple; inject $R/recB b.simple $R/recB/labels_B.txt
step recB recB "s#nstates=[0-9.]*#nstates=$KB#"
cp $R/recB/b.simple $R/stageB/b.simple
step stageB stageB "s#$NODE/recB/#$R/recB/#g; s#nstates=[0-9.]*#nstates=$KB#$drop; s#\$#$XB#"
cp $R/stageB/b.simple $R/split/s.simple
step split split "s#x#x#"
grep -E 'contrast-aware|membership by the choose rule|K\* = ' $R/split/step.log | cut -c1-160
simple_exec prg=print_project_field projfile=$R/stageB/b.simple oritype=ptcl3D mkdir=no > $R/stageB/b_ptcl3D.txt 2> /dev/null
ln -sfn $R/split/split_children.txt $R/stageB/split_children.txt
ln -sfn $A/drift/${KEY}_b.txt $R/prod/b_ptcl3D.txt; ln -sfn $NODE/split/split_children.txt $R/prod/split_children.txt
echo "== readout"
$PY $L/f2r_ab_readout.py --node $NODE --truth $TRUTH --arm production=$R/prod --arm ${ARM}_s$SEED=$R/stageB 2>&1 | tee $R/readout/labels.txt
$PY $A/pa_drift.py ${KEY}_${ARM}_s$SEED $A/drift/${KEY%_r?}_r_b.txt $A/drift/${KEY}_n.txt $R/stageA/a_ptcl3D.txt $R/stageB/b_ptcl3D.txt 1.464 $NODE/rows.txt 2>&1 | grep -vE "^parent|^  |incoming" | tee $R/readout/drift.txt
if [[ $KEY == *r131* ]]; then
  : > $R/readout/spec.txt
  for e in $R/stageB/recvol_state??_even_unfil.mrc $NODE/stageB/recvol_state??_even_unfil.mrc; do
    o=${e/_even_/_odd_}; k=$(basename $e | cut -c13-14); tag=$([[ $e == $R/* ]] && echo ${ARM}_s${SEED} || echo production)__s$k
    [ -e $o ] && echo "$tag $e $o" >> $R/readout/spec.txt
  done
  CL=$(ls $R160/docked/conf18_*_m160.mrc | head -1)
  $PY $S/make_view.py $R/readout/views $R/readout/spec.txt --box 160 > /dev/null 2>&1
  $PY $S/arm_readout.py $R160/toy_consensus.mrc $CL $R/readout/arm.txt $R/readout/views/*.mrc > /dev/null 2>&1
  echo "closed-arm / open-arm occupancy per stage-B map (closed corner >= 0.95 / <= 0.6):"; grep -v '^#' $R/readout/arm.txt | awk '{n=split($1,p,"/"); print "   " p[n], $2, $3, $4}'
fi
echo "done $(date)"
