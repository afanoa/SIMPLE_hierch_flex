#!/bin/bash
#SBATCH -J replay3
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent/logs/%x_%A_%a.out
# latent-seeding test: a gate node replayed from seed labels with its own saved step calls (seeds reconstruct3D, stage A
# refine3D, recB reconstruct3D, stage B refine3D, state_split), swapping only the project (the node's n.simple with the arm's
# labels injected by prg=selection infile=) and the starting maps; the stage-A labels fold and compact as the tree's do
# (f2r_fold_states.py). The masters submit their parts as the tree's did. Readouts: stage-B labels and split children against
# the truth and production's stage B (f2r_ab_readout.py); on r131 the gates' closed-arm / open-arm map readout of the stage-B
# maps beside production's. LIST lines: "KEY NODE ARM SEED LABELS TRUTH" (LABELS one integer per node row; SEED = SIMPLE_SEED).
set -uo pipefail
: ${LIST:?}
L=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent
B=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/merge_master2/install
S=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/fixtures/scripts
R160=/mnt/beegfs/elmlund/afan/routes/toy55/maps160
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$B/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 SIMPLE_PATH=$B
export OMP_NUM_THREADS=8 PPINFO=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo PYTHONDONTWRITEBYTECODE=1
read -r KEY NODE ARM SEED LABELS TRUTH < <(sed -n "${SLURM_ARRAY_TASK_ID:-1}p" $LIST)
export SIMPLE_SEED=$SEED
R=$L/ab_seed/${KEY}_${ARM}_s$SEED
[ -e $R ] && { echo "EXISTS $R"; exit 3; }
test -s $LABELS || { echo "no labels $LABELS"; exit 3; }
test -s $L/dumps/$KEY/b_ptcl3D.txt || { echo "no production stage-B dump for $KEY"; exit 3; }
K=$(grep -v '^0$' $LABELS | sort -n | uniq | wc -l)
KP=$(grep -v '^0$' $NODE/seeds/seeds.txt | sort -n | uniq | wc -l)
[ "$K" = "$KP" ] || { echo "labels have $K states, production's seeds $KP"; exit 3; }
mkdir -p $R/seeds $R/stageA $R/recB $R/stageB $R/split $R/prod $R/readout
echo "replay $KEY $ARM seed $SEED: K $K, $(wc -l < $LABELS) labels; $(date)"
BSM=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/merge_sm/install
inject() { ( cd $1 && $BSM/bin/simple_exec prg=write_states projfile=$2 infile=$3 mkdir=no > $2.sel.log 2>&1 ) || echo "write_states failed in $1"; }
step() { # DIR SAVED SEDX: the saved step's simple_exec call with the binary and SEDX swaps, run in $R/DIR
  local d=$1 src=$2 sedx=$3 cmd
  cd $R/$d
  cmd=$(grep -E 'simple_exec prg=' $NODE/$src/hetero_step.sh | sed -e "s#^[^ ]*/simple_exec#$B/bin/simple_exec#" -e 's# > hetero_step.log.*##' -e "$sedx")
  echo "$cmd" > cmd.txt
  echo "== $d $(date +%T)"; eval "$cmd" > step.log 2>&1; echo "   exit $? $(date +%T)"
}
cp $NODE/n.simple $R/seeds/s.simple;  inject $R/seeds s.simple $LABELS
step seeds seeds "s#nstates=[0-9.]*#nstates=$K#"
cp $NODE/n.simple $R/stageA/r.simple; inject $R/stageA r.simple $LABELS
step stageA stageA "s#$NODE/seeds/#$R/seeds/#g"
simple_exec prg=print_project_field projfile=$R/stageA/r.simple oritype=ptcl3D mkdir=no > $R/stageA/a_ptcl3D.txt 2> /dev/null
KB=$($PY $L/f2r_fold_states.py $R/stageA/a_ptcl3D.txt $R/recB/labels_B.txt)
echo "stage A states fold to $KB"
drop=""; for n in $(seq $((KB + 1)) 9); do drop="$drop; s# vol$n=[^ ]*##"; done
cp $R/stageA/r.simple $R/recB/b.simple; inject $R/recB b.simple $R/recB/labels_B.txt
step recB recB "s#nstates=[0-9.]*#nstates=$KB#"
cp $R/recB/b.simple $R/stageB/b.simple
step stageB stageB "s#$NODE/recB/#$R/recB/#g; s#nstates=[0-9.]*#nstates=$KB#$drop"
cp $R/stageB/b.simple $R/split/s.simple
step split split "s#x#x#"
simple_exec prg=print_project_field projfile=$R/stageB/b.simple oritype=ptcl3D mkdir=no > $R/stageB/b_ptcl3D.txt 2> /dev/null
ln -sfn $R/split/split_children.txt $R/stageB/split_children.txt
ln -sfn $L/dumps/$KEY/b_ptcl3D.txt $R/prod/b_ptcl3D.txt; ln -sfn $NODE/split/split_children.txt $R/prod/split_children.txt
echo "== readout"
$PY $L/f2r_ab_readout.py --node $NODE --truth $TRUTH --arm production=$R/prod --arm ${ARM}_s$SEED=$R/stageB 2>&1 | tee $R/readout/labels.txt
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
