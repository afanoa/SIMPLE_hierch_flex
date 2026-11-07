#!/bin/bash
#SBATCH -J klem_umaps_now
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Embeddings for the status report: (1) toy, truth-coloured UMAPs of the parent at consensus poses vs the realigned
# balanced parent; (2) real family-1 feedback board, UMAP of the responsibilities coloured by the matched deposited
# conformer (and the usual panels), the consensus-pose child as the reproducibility partner; (3) closed-state montage:
# consensus, our conformer-18/46/10 maps and the deposited 18/46/10 at 8 A.
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; P=$A/rnasep_klem/expand1/fam1_pipeline; F=$W/figs_now; mkdir -p $F; cd $W
echo "== toy UMAPs"; $PY $W/toy_umap_compare.py $W/sim/gt_particles.txt "$W/lat/klem_junk:parent, consensus poses (0.48)" "$W/lat/klem_realign1_bal05:parent, realigned poses + balanced E-step (0.61)" $F/umap_toy_parent_vs_realigned_bal.png --independent 2>&1 | grep -v Warning
echo "== real labels"; $PY - <<PYEOF
import numpy as np
ids = [1, 10, 18, 27, 28, 37, 46, 48, 49, 51]; st = np.loadtxt('$P/klem_fb/matched_rec/states.txt', dtype=int); lab = np.where(st > 0, np.array([0] + ids)[st], 0); np.savetxt('$F/matched_conformer_labels.txt', lab, fmt='%d'); print('labels: ' + ' '.join('%d:%d' % (i, (lab == i).sum()) for i in ids))
PYEOF
echo "== real UMAP (family-1 feedback board, matched conformer colours)"; $PY $W/klem_umap.py $P/klem_fb $P/meta_fb.npz --labels $F/matched_conformer_labels.txt --other $A/rnasep_klem/expand1/fam1 --out $F/umap_real_fam1_fb_matched.png --out2 $F/umap_real_fam1_flex_vs_klem.png 2>&1 | grep -v Warning
echo "== closed-state montage"; M=$P/klem_fb/matched_rec/maps; D=$W/maps160/docked
$PY $W/montage_maps.py $F/closed_state_maps.png --lp 8 "$W/maps160/toy_consensus.mrc:consensus" "$M/recvol_state03.mrc:ours 18 (arm 0.76)" "$(ls $D/conf18_*.mrc | head -1):deposited 18 (closed)" "$M/recvol_state07.mrc:ours 46" "$(ls $D/conf46_*.mrc | head -1):deposited 46" "$M/recvol_state02.mrc:ours 10" "$(ls $D/conf10_*.mrc | head -1):deposited 10" "$M/recvol_state04.mrc:ours 27 (open)" "$(ls $D/conf27_*.mrc | head -1):deposited 27" 2>&1 | grep -v Warning
echo "[$(date '+%F %T')] umaps DONE"
