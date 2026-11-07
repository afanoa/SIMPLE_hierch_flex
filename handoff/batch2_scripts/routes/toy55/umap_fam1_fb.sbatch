#!/bin/bash
#SBATCH -J umap_fam1_fb
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 01:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# UMAPs of a real-data board: responsibility vectors coloured by the deposited conformer matched by the argmax cell
# (identification readout, not truth), posterior means, junk mass, defocus, and a second board's z1 (reproducibility);
# plus the flex latent of the same particles coloured by the board. env: RTAG (f692), FAM (1), BOARD (klem_fb), OTHER (klem_fb_distr6/board)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/${RTAG:-f692}; D=$R/fam${FAM:-1}; B=$D/${BOARD:-klem_fb}; O=$D/${OTHER:-klem_fb_distr6/board}
$PY - <<PYEOF
import numpy as np
M = np.loadtxt('$B/cell_matches.txt', comments='#'); best = {int(r[0]): int(r[1]) for r in M}
L = np.loadtxt('$B/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int); am = L[:, 3].astype(int); gm = L[:, 5]
N = len(np.load('$D/meta_fb.npz')['flags']); lab = np.zeros(N, int)
for p, c, g in zip(pind, am, gm):
    if g > 0.5 and c in best: lab[p - 1] = best[c]
np.savetxt('$B/labels_conformer.txt', lab, fmt='%d'); u, n = np.unique(lab[lab > 0], return_counts=True); print('conformer labels (argmax cell -> best deposited match):', dict(zip(u.tolist(), n.tolist())))
PYEOF
OTH=""; [ -s $O/klem_latent.txt ] && OTH="--other $O"
$PY $W/klem_umap.py $B $D/meta_fb.npz --labels $B/labels_conformer.txt --label_title "deposited conformer matched by the argmax cell (readout, not truth)" --label_name conformer $OTH --out $D/fig_umap_fam${FAM:-1}_${BOARD:-klem_fb}.png --out2 $D/fig_umap_fam${FAM:-1}_${BOARD:-klem_fb}_flexlatent.png 2>&1 | grep -v Warning | tail -6
echo "UMAP DONE"
