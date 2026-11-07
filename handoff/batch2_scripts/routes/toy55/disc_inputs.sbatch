#!/bin/bash
#SBATCH -J klem_disc_inputs
#SBATCH -p csbdevel
#SBATCH -c 32
#SBATCH --mem=200G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# discovery on pre-stage-A inputs (fflex states, cheap children) vs recB per run, and the resulting closed pole (disc_inputs.py); env RT NODE BASE BASEPOST SEED PORTS OUT H CP CTFABS MSKDIAM PARS
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=4
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; R=$K_/rflex_${RT:?}
$PY -u $T/disc_inputs.py --base ${BASE:?} --basepost ${BASEPOST:?} --seed ${SEED:?} --tree $R --node ${NODE:-r1} --out ${OUT:?} --h ${H:?} --cache ${CP:?} --ctfabs ${CTFABS:-1} --mskdiam ${MSKDIAM:-232} \
  --pars ${PARS:-4} --ports ${PORTS:-} 2>&1 | grep -v Warning
