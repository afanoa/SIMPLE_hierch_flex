#!/bin/bash
#SBATCH -p csbdevel
#SBATCH --mem=96G
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# RNase P readout of the switch / swing call on the full692_m33_s1 tree. env NODE (r, r1, r2) and SCORES: 'pole' (the kept pole's
# endpoint posteriors, null_score.py) and / or endpoint names (z1hi ..: that endpoint's site posterior). Per score: 5 at-node bins
# over all node rows (the junk screen was off, so the node has no junk rows; motion_bins5.py), the node's own project and recB0
# reconstruct3D with nstates=5, postprocess=no; phaseE_maptests.py in the node's mask. Scores' Spearman matrix in scores.txt
set -uo pipefail
NP=${NP:-4}; NT=${NT:-8}; TAG=m33_${NODE:?}; source /mnt/beegfs/elmlund/afan/routes/toy55/seriesE_lib.sh
H=$S/full692_m33_s1/1_hetero_tree; N=$H/node_$NODE; OUT=$K/factorial/phaseE/rnasep_m33/node_$NODE; mkdir -p $OUT/scores; CS=""
for sc in ${SCORES:?}; do
  if [ $sc = pole ]; then $PY $T/null_score.py $N $OUT/scores/pole.txt > $OUT/scores/pole_rule.txt 2>&1 || { log "pole score FAILED: $(tail -1 $OUT/scores/pole_rule.txt)"; continue; }
  else cp $N/site/site_post_$sc.txt $OUT/scores/$sc.txt; fi; CS="$CS $sc=$OUT/scores/$sc.txt"
done
E=""; for e in z1lo z1hi z3lo z3hi; do E="$E $e=$N/site/site_post_$e.txt"; done
$PY $T/score_corr.py $CS $E > $OUT/scores.txt 2>&1
for sc in $SCORES; do [ -s $OUT/scores/$sc.txt ] || continue; nm=m33_${NODE}_$sc; newdir $nm
  { $PY $T/motion_bins5.py $OUT/scores/$sc.txt states.txt none - > ../bins.txt 2>&1 && recon $nm $N/n.simple $N/recB0/hetero_step.sh $N/mask.mrc; } || log "$nm FAILED: $(tail -1 ../bins.txt)"
done
finish
