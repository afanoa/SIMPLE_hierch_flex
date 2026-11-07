#!/bin/bash
#SBATCH -J klem_paths
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# paths between motion extremes (paths.py) on the IgM trees: state paths (latent kNN / map MST), continuous-score bins, site coverage; env OUT
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; S=$K/stree_port
AT=$K/factorial/armsigm1_truth/arms_truth.txt; OUT=${OUT:-$K/factorial/validation/paths_igm}
spec(){ D=$S/$2/1_hetero_tree; echo "$1:$D:factorial/table.txt:$D/node_r/flex/flex_pca_embedding.bin:$D/node_r/fflex/flex_pca_embedding.bin"; }
$PY $T/paths.py $AT $OUT $(spec m52M_igm fulltoy_m52M_igm) $(spec igm_m8c17 fulltoy_igm) 2>&1 | grep -v Warning | tee $OUT.txt
