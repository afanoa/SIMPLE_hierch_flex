#!/bin/bash
#SBATCH -J md_proto
#SBATCH -p csbdevel
#SBATCH --constraint="x6342|x6740P"
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 01:30:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/mask_disc/logs/%x_%A_%a.out
# one node of the mask / discovery / stop prototype: line $SLURM_ARRAY_TASK_ID of SPEC = "TAG NODE MSKDIAM [cheap=F|dump] [parent=DIR] [arm=r131|arms]"
set -uo pipefail
: ${SPEC:?}
W=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/mask_disc; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55
PY=/mnt/nasapps/production/python/3.12/bin/python3; export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
BIN=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/build_dump/install/bin
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$BIN:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SPEC); set -- $line; TAG=$1; NODE=$2; MD=$3; shift 3
O=$W/proto/$TAG; mkdir -p $O/dump; cd $O/dump
if [ "$MD" = auto ]; then MD=$(grep -o 'mskdiam=[0-9.]*' $NODE/../../run_manifest.txt | head -1 | cut -d= -f2); fi
X=""
dump() { # dump STATES_FILE PROJECT: one state per row
  $BIN/simple_exec prg=print_project_field projfile=$2 oritype=ptcl3D mkdir=no 2>/dev/null | grep -oE "(^| )state=[0-9]+" | cut -d= -f2 > $1
  echo "dumped $(grep -c . $1) states from $2"; }
for kv in "$@"; do
  case $kv in
    cheap=dump) dump $O/dump/cheap_states.txt $NODE/cheap/c.simple; X="$X --labels-cheap $O/dump/cheap_states.txt" ;;
    cheap=*) X="$X --labels-cheap ${kv#cheap=}" ;;
    parent=*) X="$X --parent ${kv#parent=}" ;;
    arm=r131) X="$X --cons $T/maps160/toy_consensus.mrc --closed $(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)" ;;
    arm=arms) X="$X --truth $A/rnasep_klem/synth_arms/arm1_travel.mrc,$A/rnasep_klem/synth_arms/arm2_travel.mrc,$A/rnasep_klem/synth_arms/arm3_travel.mrc,$A/rnasep_klem/synth_arms/arm4_travel.mrc" ;;
    nodisc) X="$X --no-disc" ;;
  esac
done
if [ -e $NODE/stageB/b.simple ]; then dump $O/dump/stageb_states.txt $NODE/stageB/b.simple; X="$X --labels-stageb $O/dump/stageb_states.txt"; fi
echo "node $NODE mskdiam $MD args $X"
$PY -u $W/scripts/vmask_proto.py --node $NODE --tag $TAG --out $O --mskdiam $MD $X 2>&1 | grep -v Warning
