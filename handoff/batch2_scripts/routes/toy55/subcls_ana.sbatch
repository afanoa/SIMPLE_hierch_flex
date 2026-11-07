#!/bin/bash
#SBATCH -J subcls_ana
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Analysis-only pass over an existing sub-class node run (subcls_node*.sbatch output O): SCRIPT (subcls_embed3.py) with the node's
# inputs (rows and flex coordinates from the run's queue / node dir, labels and oris from O), truth by TAG + CODES (colon form) or TTABLE.
# env RUN, NODE, O, TAG (swing30), CODES (swing:200:230), TTABLE, SCRIPT (subcls_embed3.py), EOUT (embed3), EXTRA (space-separated extra args)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
RUN=${RUN:?}; NODE=${NODE:?}; O=${O:?}; TAG=${TAG-swing30}; CODES=${CODES-swing:200:230}; SCRIPT=${SCRIPT:-subcls_embed3.py}; EOUT=${EOUT:-embed3}
RD=$K_/$RUN; ND=$RD/node_$NODE; set -- $(awk -v u=$NODE '$1 == u' $RD/queue.txt); ROWS=${4:-}; FLEX=$(ls -td $ND/*_flex_pca | head -1)/flex_pca_coordinates.txt
REFS=""; [ -s $ND/refined_states.txt ] && REFS="rflex_refined_states=$ND/refined_states.txt,rflex_flex_states=$ND/flex_states.txt"; [ -s $RD/coord_$NODE/leaf_states.txt ] && REFS="$REFS,rflex_coord_leaves=$RD/coord_$NODE/leaf_states.txt"
TRU=""; [ -n "$TAG" ] && TRU="--gt $T/sim$TAG/gt_particles.txt $([ -n "$CODES" ] && echo --codes $CODES) $([ -n "${TTABLE:-}" ] && echo --ttable $TTABLE) $([ -n "$REFS" ] && echo --refs ${REFS#,})"
echo "[$(date '+%F %T') subcls_ana $RUN $NODE] $SCRIPT on $O -> $O/$EOUT ${EXTRA:-}"
$PY $T/$SCRIPT --labels $O/p/cls_expansion_labels.txt --flex $FLEX --rows $ROWS --oris $O/oris_node.txt $TRU --out $O/$EOUT ${EXTRA:-} 2>&1 | grep -v Warning
echo "[$(date '+%F %T') subcls_ana $RUN $NODE] DONE"
