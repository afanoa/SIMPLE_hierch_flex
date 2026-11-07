#!/bin/bash
#SBATCH -J klem_latent_cca
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# flex latent stability across runs of node f692st r1 (latent_cca.py): whole flex script vs port; focused flex script vs port run 1 vs port
# rerun (seed2) vs the script's earlier cheap-mask run
K_=/mnt/beegfs/elmlund/afan/rnasep_klem; S=$K_/rflex_f692st/node_r1; P1=$K_/stree_port/tree_f692r1/1_hetero_tree/node_r; P2=$K_/stree_port/tree_f692r1_seed2/1_hetero_tree/node_r
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4; PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55
echo "== whole flex"; $PY $T/latent_cca.py script=$S/1_flex_pca/flex_pca_coordinates.txt port1=$P1/flex/flex_pca_coordinates.txt 2>&1 | grep -v Warning
echo "== focused flex"; $PY $T/latent_cca.py script=$S/focus/1_flex_pca/flex_pca_coordinates.txt port1=$P1/fflex/flex_pca_coordinates.txt \
  port2=$P2/fflex/flex_pca_coordinates.txt oldmask=$S/old_cheap20/focus/1_flex_pca/flex_pca_coordinates.txt 2>&1 | grep -v Warning
