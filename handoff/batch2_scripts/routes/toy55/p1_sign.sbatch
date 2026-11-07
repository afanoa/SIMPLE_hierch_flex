#!/bin/bash
#SBATCH -J klem_p1sign
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# signedness of the end-to-end difference of the rflex_s30a coordinate nodes (problem 1, p1_sign.py) with truth
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; R=/mnt/beegfs/elmlund/afan/rnasep_klem/rflex_s30a
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
NODES=""; for nd in r22 r11 r12 r21; do NODES="$NODES --node $nd:$(ls -td $R/coord_$nd/leaves/*_reconstruct3D | head -1):$R/coord_$nd/leaf_states.txt"; done
$PY $T/p1_sign.py /mnt/beegfs/elmlund/afan/rnasep_klem/p1sign $NODES --gt $T/simswing30/gt_particles.txt --truth $T/simswing30/maps 2>&1 | grep -v Warning
echo DONE
