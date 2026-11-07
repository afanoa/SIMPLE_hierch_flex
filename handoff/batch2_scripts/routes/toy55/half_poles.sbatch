#!/bin/bash
#SBATCH -J klem_half_poles
#SBATCH -p csbdevel
#SBATCH -c 32
#SBATCH --mem=200G
#SBATCH -t 08:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# discovery + poles independently on two disjoint halves of a node (half_poles.py); env RT NODE OUT H CP CTFABS MSKDIAM NSN NPOLE MAPN PARS
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=${QPART:-csbdevel} OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=4
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; R=$K_/rflex_${RT:?}
$PY -u $T/half_poles.py --tree $R --node ${NODE:-r1} --out ${OUT:?} --h ${H:?} --cache ${CP:?} --ctfabs ${CTFABS:-1} --mskdiam ${MSKDIAM:-232} --nsn ${NSN:-3} --npole ${NPOLE:-7500} --mapn ${MAPN:-1250} \
  --pars ${PARS:-4} --cons $T/maps160/toy_consensus.mrc --closed $(ls $T/maps160/docked/conf18_*_m160.mrc | head -1) 2>&1 | grep -v Warning
