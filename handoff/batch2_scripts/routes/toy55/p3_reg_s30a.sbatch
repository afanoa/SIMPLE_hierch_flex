#!/bin/bash
#SBATCH -J klem_p3reg
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# rigid re-registration of the rflex_s30a flex states (true support known) before the member mask (problem 3, p3_masks2.py --reg 1)
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem; O=$K_/p3masks/reg; mkdir -p $O
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
for nd in r22 r2; do $PY $T/p3_masks2.py $O --name s30a_$nd --list $K_/rflex_s30a/node_$nd/members.txt --labels $K_/rflex_s30a/node_$nd/flex_states.txt --gt $T/simswing30/gt_particles.txt --truth $T/simswing30/maps --mskdiam 220 --reg 1 > $O/log_$nd.txt 2>&1 & done
wait; cat $O/*_masks2.txt; echo DONE
