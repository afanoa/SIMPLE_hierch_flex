#!/bin/bash
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# OLE reference calls for the port's node_motion. env TREE (6_hetero_tree | 7_hetero_tree), SEED (a tag), CANDS ('node:site ..').
# Per candidate site: site_score.py's pole score, 5 at-node bins over the node's rows (motion_bins5.py), the node's own project, recB0
# reconstruct3D with nstates=5, postprocess=no, box_crop=200 and mask (seriesF_lib.sh); then phaseE_maptests.py, phaseF_rel.py and
# view_check.py (C2) on each, in the background
set -uo pipefail
NP=2; NT=8; TAG=ole_${SEED:?}; REC_EXTRA="box_crop=200"; source /mnt/beegfs/elmlund/afan/routes/toy55/seriesF_lib.sh
O=/mnt/beegfs/elmlund/afan/ole_12707; H=$O/${TREE:?}; OUT=$O/motion_calls/$SEED; mkdir -p $OUT
for c in ${CANDS:?}; do nd=${c%%:*}; k=${c#*:}; N=$H/node_$nd
  if [ $k = 1 ]; then SD=$N/site; else SD=$N/nsite_${nd}_s$k; [ -d $SD ] || SD=$N/fsite_${nd}_s$k; fi
  nm=ole_${SEED}_${nd}_s$k; newdir $nm
  { $PY $T/site_score.py $SD $N score.txt > ../score_rule.txt 2>&1 && $PY $T/motion_bins5.py score.txt states.txt none - > ../bins.txt 2>&1 \
    && recon $nm $N/n.simple $N/recB0/hetero_step.sh $N/mask.mrc; } || { log "$nm FAILED: $(tail -qn1 ../score_rule.txt ../bins.txt 2>&1)"; continue; }
  ( OMP_NUM_THREADS=1 $PY $T/phaseF_rel.py $nm:$W:$N/mask.mrc > $W/relF.txt 2>&1
    cd $W/work && simple_exec prg=print_project_field projfile=p.simple oritype=ptcl3D 2>/dev/null | $PY $T/view_check.py $nm 2 > $W/view.txt 2>&1 ) &
done
finish
cat $OUT/*/relF.txt 2>/dev/null | grep "^RELF" > $OUT/relF.txt; cat $OUT/*/view.txt 2>/dev/null | grep "^VIEW" > $OUT/view_check.txt
log "phaseF $(wc -l < $OUT/relF.txt), views $(wc -l < $OUT/view_check.txt)"
