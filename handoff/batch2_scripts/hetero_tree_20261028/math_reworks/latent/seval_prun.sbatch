#!/bin/bash
#SBATCH -J sevprun
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 01:30:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent/logs/%x_%A_%a.out
# refine=state_eval with the node's flex prior (seval_prior.py; SIMPLE_STATE_PRIOR at weight ALPHA, none = the objective alone): the node's
# stage-B call on a copy of its final project and maps, MAXITS its; seval_score.py vs stage B and the truth. LIST: "KEY NODE MAXITS ALPHA TRUTHSPEC".
set -uo pipefail
: ${LIST:?}
L=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent
SE=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/state_eval/install
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$SE/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 SIMPLE_PATH=$SE
read -r KEY NODE MAXITS ALPHA TRUTH < <(sed -n "${SLURM_ARRAY_TASK_ID:-1}p" $LIST)
D=$L/seval_prior_runs/${KEY}_a${ALPHA}_m$MAXITS; mkdir -p $D/stderrout; cd $D
if [ "$ALPHA" != none ]; then export SIMPLE_STATE_PRIOR=$L/seval_prior/$KEY.txt SIMPLE_STATE_PRIOR_ALPHA=$ALPHA; test -s $SIMPLE_STATE_PRIOR || { echo "no prior file"; exit 2; }; fi
cp $NODE/stageB/b.simple s.simple
K=$(ls $NODE/stageB/recvol_state[0-9][0-9].mrc | wc -l)
VOLS=$(for s in $(seq 1 $K); do printf 'vol%d=%s/stageB/recvol_state%02d.mrc ' $s $NODE $s; done)
CMD=$(grep -o 'simple_exec prg=refine3D .*' $NODE/stageB/hetero_step.sh | sed 's/ > hetero_step.log.*//' | sed -E "s/ vol[0-9]+=[^ ]+//g; s/ refine=[^ ]+/ refine=state_eval/; s/ prob_neigh_mode=[^ ]+//; s/ prob_beta=[^ ]+//; s/ maxits=[^ ]+/ maxits=$MAXITS/; s/ projfile=[^ ]+/ projfile=s.simple/")
CMD="$SE/bin/${CMD} $VOLS"
echo "$CMD" > seval_cmd.txt
dump(){ simple_exec prg=print_project_field projfile=s.simple oritype=ptcl3D mkdir=no 2>/dev/null | sed 's/\x1b\[[0-9;]*m//g' | grep 'e1=' > $1; }
dump before.txt
OMP_NUM_THREADS=8 /usr/bin/time -f "refine3D wall %e s" $CMD > seval.log 2> seval_time.txt; rc=$?
echo "state_eval $KEY alpha $ALPHA ($K states, $MAXITS its) exit $rc: $(tail -1 seval_time.txt)"
[ $rc = 0 ] || { grep -a -m3 'ERROR' seval.log; exit $rc; }
dump after.txt
export PPINFO=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo OPENBLAS_NUM_THREADS=1
$PY $L/seval_score.py $NODE before.txt after.txt $TRUTH
