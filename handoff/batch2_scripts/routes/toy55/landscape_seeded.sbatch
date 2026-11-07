#!/bin/bash
#SBATCH -J klem_landscape_seeded
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# map-level reproducibility of the pole-seeded K=2 and K=3 stage-A states across script / port1 / port2 poles (landscape_cmp.py)
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8; K=/mnt/beegfs/elmlund/afan/rnasep_klem; PY=/mnt/nasapps/production/python/3.12/bin/python3
for p in p2 p3; do A=""; for r in script port1 port2; do d=$K/repro_r1/seeded_${p}_$r; A="$A $r=$(ls -td $d/stageA/*_refine3D | head -1):$d/stageA_states.txt"; done
  echo "######## pole-seeded $p"; $PY /mnt/beegfs/elmlund/afan/routes/toy55/landscape_cmp.py $A --lp 15 2>&1; done > $K/repro_r1/landscape_seeded.txt; cat $K/repro_r1/landscape_seeded.txt
