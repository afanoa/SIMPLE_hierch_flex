#!/bin/bash
#SBATCH -J klem_arm_int
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:15:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# closed-arm / open-arm readouts of the integrated (pole-seeded) hetero_tree r1 final maps, runs poles_r1_p1/p2/p3 (raw and emp2, 10 A)
export OMP_NUM_THREADS=2; PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; P=/mnt/beegfs/elmlund/afan/rnasep_klem/stree_port
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
for r in p1 p2 p3; do F=$P/poles_r1_$r/1_hetero_tree/final_maps; for k in raw emp2; do echo "######## $r $k: $(grep -v '^#' $F/units_index.txt | awk '{printf "%s %s %d | ", $1, $2, $3}')"
  $PY $T/arm_readout_any.py $T/maps160/toy_consensus.mrc $CL $F/$k/unit_0?.mrc 2>&1 | grep -v Warning; done; done
