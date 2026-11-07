#!/bin/bash
#SBATCH -J routes_realprops
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
set -uo pipefail
export LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64:${LD_LIBRARY_PATH:-} OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; cd $W
DATA=$(ls -d $A/rnasep_coupled/basis_span_all/data $A/rnasep_coupled/latent2d3d/cls_full 2>/dev/null | head -1); echo "data dir: $DATA"; ls $DATA | head -3
$PY $W/real_props.py $DATA $A/rnasep_coupled/avgflex/cons_box80.mrc $W/props --nptcl 4000
