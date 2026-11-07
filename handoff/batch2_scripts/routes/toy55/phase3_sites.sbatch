#!/bin/bash
#SBATCH -J klem_phase3
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Phase 3 head-to-head on a synthetic set: the same scorecard (score_tree_sites.py) for (1) the tree's leaves, (2) the site layer's
# combination classes, (3) the hybrid = tree leaves x site classes (hybrid_labels.py, pieces >= MINN), then the side-by-side table
# (compare_tree_sites.py). env RT (tree), SITES (site-layer OUT dir), TRUTH, TMAPS / TMASKS (optional), ORACLE (oracle .out, optional), OUT
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; R=$K_/rflex_${RT:?}; S=${SITES:?}; O=${OUT:-$S/phase3}; mkdir -p $O
MX=""; [ -n "${TMAPS:-}" ] && MX="--maps $TMAPS --masks $TMASKS"
$PY $T/score_tree_sites.py --tree $R --truth $TRUTH $MX --minn 500 --json $O/tree.json > $O/tree.txt 2>&1
$PY $T/score_tree_sites.py --tree $R --truth $TRUTH --labels $S/combine/combo_labels.txt --leafmaps $S/leafmaps.txt $MX --minn 500 --json $O/sites.json > $O/sites.txt 2>&1
SL=$(ls $S/combine/site_s*_labels.txt | sort -V | tr '\n' ',' | sed 's/,$//')
$PY $T/hybrid_labels.py $R/final_labels.txt $SL $O/hybrid_labels.txt --minn ${MINN:-500}
$PY $T/score_tree_sites.py --tree $R --truth $TRUTH --labels $O/hybrid_labels.txt --minn 500 --json $O/hybrid.json > $O/hybrid.txt 2>&1
$PY $T/compare_tree_sites.py tree:$O/tree.json sitelayer:$O/sites.json hybrid:$O/hybrid.json ${ORACLE:+--oracle $ORACLE} 2>&1 | grep -v Warning
