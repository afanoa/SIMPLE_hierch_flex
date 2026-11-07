#!/bin/bash
#SBATCH -J klem_mask_series
#SBATCH -p csbdevel
#SBATCH -c 32
#SBATCH --mem=200G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# closed-pole mask sensitivity (mask_series.py): fixed start, shared input poses, only the site mask varies; env RT NODE BASE SEED PORTS OUT H CP CTFABS MSKDIAM PARS
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=4
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; R=$K_/rflex_${RT:?}
$PY -u $T/mask_series.py --base ${BASE:?} --seed ${SEED:?} --tree $R --node ${NODE:-r1} --out ${OUT:?} --h ${H:?} --cache ${CP:?} --ctfabs ${CTFABS:-1} --mskdiam ${MSKDIAM:-232} --pars ${PARS:-4} \
  --cons $T/maps160/toy_consensus.mrc --closed $(ls $T/maps160/docked/conf18_*_m160.mrc | head -1) --ports ${PORTS:-} 2>&1 | grep -v Warning
