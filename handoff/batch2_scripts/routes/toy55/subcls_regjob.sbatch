#!/bin/bash
#SBATCH -J subcls_reg
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# subcls_reg.py on a sub-class node run made with DUMP=1 (O/p holds cls_expansion_latent_class*.txt): per-bin 2D latent -> 3D latent
# regression (cross-fitted), k-means windows, truth. env RUN, NODE, O, TAG (swing30), CODES (swing:200:230), TTABLE, EXTRA
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
RUN=${RUN:?}; NODE=${NODE:?}; O=${O:?}; TAG=${TAG-swing30}; CODES=${CODES-swing:200:230}
RD=$K_/$RUN; ND=$RD/node_$NODE; set -- $(awk -v u=$NODE '$1 == u' $RD/queue.txt); ROWS=${4:-}; FLEX=$(ls -td $ND/*_flex_pca | head -1)/flex_pca_coordinates.txt
TRU=""; [ -n "$TAG" ] && TRU="--gt $T/sim$TAG/gt_particles.txt $([ -n "$CODES" ] && echo --codes $CODES) $([ -n "${TTABLE:-}" ] && echo --ttable $TTABLE)"
echo "[$(date '+%F %T') subcls_reg $RUN $NODE] $O"
$PY $T/subcls_reg.py --dumps $O/p --labels $O/p/cls_expansion_labels.txt --flex $FLEX --rows $ROWS --out $O/reg $TRU ${EXTRA:-} 2>&1 | grep -v Warning
echo "[$(date '+%F %T') subcls_reg $RUN $NODE] DONE"
