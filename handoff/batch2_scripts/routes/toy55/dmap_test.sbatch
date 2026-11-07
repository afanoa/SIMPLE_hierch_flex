#!/bin/bash
#SBATCH -J klem_dmapA
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Phase A, discrete vs continuous by diffusion maps + sweep conductance (dmap_test.py) on the IgM M, CAHRA M s1 and cdisc M roots; then flex's GMM AUTO lines
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; S=$K/stree_port
O=$K/factorial/validation/dmap_phaseA; mkdir -p $O; cd $O
e(){ D=$S/$1/1_hetero_tree/node_r; echo "$D/flex/flex_pca_embedding.bin:$D/fflex/flex_pca_embedding.bin"; }
$PY $T/dmap_test.py igm_M:igm:$K/factorial/armsigm1_truth/arms_truth.txt:$(e fulltoy_m52M_igm) cahra_Ms1:cahra:/mnt/beegfs/elmlund/EMPIAR_downloads/13991/data/source.star:$(e arm_M_c11_s1) \
  cdisc_M:cdisc:$K/factorial/junk_test/simsites_truth_jt.txt:$(e fulltoy_m52M_cd1) 2>&1 | grep -v Warning | tee dmap_test.txt
for t in fulltoy_m52M_igm arm_M_c11_s1 fulltoy_m52M_cd1; do for l in flex fflex; do echo "######## $t $l"; grep -a "GMM AUTO" $S/$t/1_hetero_tree/node_r/$l/hetero_step.log | cut -c1-200; done; done | tee gmm_auto.txt
