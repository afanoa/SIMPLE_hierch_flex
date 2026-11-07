#!/bin/bash
#SBATCH -J klem_viz_prep
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# data for the f692rf tree and r22 UMAP figures (viz_prep.py) -> rnasep_klem/viz_f692rf/
export OMP_NUM_THREADS=16 NUMBA_NUM_THREADS=16
/mnt/nasapps/production/python/3.12/bin/python3 -u /mnt/beegfs/elmlund/afan/routes/toy55/viz_prep.py /mnt/beegfs/elmlund/afan/rnasep_klem/viz_f692rf 2>&1 | grep -v Warning
