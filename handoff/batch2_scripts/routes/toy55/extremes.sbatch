#!/bin/bash
#SBATCH -J klem_extremes
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# motion-extreme capture (extremes.py) on the IgM trees: leaves, units, site poles / score tails, root flex / fflex latent tails vs the oracle; env OUT
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; S=$K/stree_port
AT=$K/factorial/armsigm1_truth/arms_truth.txt; OUT=${OUT:-$K/factorial/validation/extremes_igm}
spec(){ D=$S/$2/1_hetero_tree; echo "$1:$D:factorial/table.txt:$D/node_r/flex/flex_pca_embedding.bin:$D/node_r/fflex/flex_pca_embedding.bin"; }
$PY $T/extremes.py $AT $OUT $(spec m52M_igm fulltoy_m52M_igm) $(spec igm_m8c17 fulltoy_igm) 2>&1 | grep -v Warning | tee $OUT.txt
