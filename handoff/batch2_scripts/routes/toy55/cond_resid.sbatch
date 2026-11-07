#!/bin/bash
#SBATCH -J klem_condres
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 00:25:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# weak IgM motions: conditional (within-branch) and residualised (population-wide) score tests (cond_resid.py) on the M IgM tree; env OUT
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem
D=$K/stree_port/fulltoy_m52M_igm/1_hetero_tree; P=$K/stree_port/promote_igm
$PY $T/cond_resid.py $K/factorial/armsigm1_truth/arms_truth.txt $D $D/node_r/fflex/flex_pca_embedding.bin \
  a5_th:node_r211122/fsite_r211122_s3/site_post_z1hi.txt:r22_s6:$P/a5_th/site_post_z1hi.txt \
  a3_th:node_r211212/fsite_r211212_s2/site_post_z3lo.txt:r1_s3:$P/a3_th/site_post_z3hi.txt \
  a5_ph:node_r22112/fsite_r22112_s5/site_post_z3hi.txt:r111_s2:$P/a5_ph/site_post_z3hi.txt 2>&1 | grep -v Warning | tee ${OUT:?}
