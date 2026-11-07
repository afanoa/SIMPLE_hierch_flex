#!/bin/bash
#SBATCH -J subcls_oracle
#SBATCH -p csbdevel
#SBATCH -c 32
#SBATCH --mem=120G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# ORACLE node: a node's members restricted by truth to one conformer range (default the swing 200..230), as a pseudo rflex run
# (rflex_oracle_<RUN>_<NODE>: queue line, node mean and flex latent linked from the real node) -> subcls_node3.sbatch (bins,
# cls_expansion with the 2D-latent dump) -> subcls_regjob.sbatch (2D vs 3D latent per bin, windows, truth), inline in this job.
# env RUN (rflex_s30a), NODE (r22), TAG (swing30), LO/HI (200/230), NBINS (80), NSUB (8)
set -uo pipefail
A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; RUN=${RUN:-rflex_s30a}; NODE=${NODE:-r22}; TAG=${TAG:-swing30}; LO=${LO:-200}; HI=${HI:-230}
RD=$K_/$RUN; ND=$RD/node_$NODE; set -- $(awk -v u=$NODE '$1 == u' $RD/queue.txt); PP=$3; ROWS=$4
OR=$K_/rflex_oracle_${RUN#rflex_}_${NODE}; mkdir -p $OR/node_${NODE}o
awk 'NR == FNR {if ($0 !~ /^#/) c[$1] = $2; next} {r++; print ($1 > 0 && c[r] >= '$LO' && c[r] <= '$HI') ? 1 : 0}' $T/sim$TAG/gt_particles.txt $ROWS > $OR/rows_${NODE}o.txt
echo "${NODE}o 2 $PP $OR/rows_${NODE}o.txt" > $OR/queue.txt
ln -sfn $ND/mean $OR/node_${NODE}o/mean; ln -sfn $(ls -td $ND/*_flex_pca | head -1) $OR/node_${NODE}o/1_flex_pca
echo "[$(date '+%F %T') oracle] $(awk '$1 > 0' $OR/rows_${NODE}o.txt | wc -l) of $(awk '$1 > 0' $ROWS | wc -l) members of $NODE have conformer $LO..$HI"
O=$K_/subcls_oracle_${RUN#rflex_}_${NODE}o_b${NBINS:-80}s${NSUB:-8}d
RUN=$(basename $OR) NODE=${NODE}o TAG=$TAG CODES=swing:$LO:$HI NBINS=${NBINS:-80} NSUB=${NSUB:-8} DUMP=1 SKIPEMBED=1 OUT=$O bash $T/subcls_node3.sbatch
RUN=$(basename $OR) NODE=${NODE}o O=$O TAG=$TAG CODES=swing:$LO:$HI bash $T/subcls_regjob.sbatch
