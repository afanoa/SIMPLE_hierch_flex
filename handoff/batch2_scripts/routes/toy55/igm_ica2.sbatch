#!/bin/bash
#SBATCH -J igm_ica2
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# ICA variants (sub-Gaussian pursuit, contrasts, whitening dims) of the IgM-like root flex latent (igm_ica2.py), scored per arm angle. Inputs from the rootprep + subcls_igm1 runs.
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
S=$K_/subcls_igm1_r_b150s8d
$PY $T/igm_ica2.py --flex $K_/rflex_root_armsigm1/node_r/1_flex_pca/flex_pca_coordinates.txt --bins $S/bins.npz --embed $S/embed3/subcls_embed.npz \
  --gt $T/simarmsigm1/gt_particles.txt --ttable $T/truth_armsigm1.txt 2>&1 | grep -v Warning
echo "[$(date '+%F %T') igm_ica] DONE"
