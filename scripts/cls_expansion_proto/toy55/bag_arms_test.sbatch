#!/bin/bash
#SBATCH -J klem_bag_arms
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Pooling by 2D sub-class inside the board: the Python board on the particles of the 5 expanded classes of the light
# arms set (cls_expansion kernel sub-classes as bags), focus = arm 1's travel, 3x3 lattice; bag prior strength 0 (per
# particle) vs 2 (sub-class pooled); readout = arm angle std ratio per cell.
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/synth_arms; O=$A/rnasep_klem/bag_arms; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') bagarms] $*"; }
$PY - <<PYEOF
import numpy as np
M = dict(np.load('$R/meta_toy.npz')); bag = np.loadtxt('$W/route2d_arms/subclass_labels.txt', dtype=int); par = np.loadtxt('$W/route2d_arms/expansion/class_cluster.txt', dtype=int)[:, 1]
n = len(M['flags']); bag = bag[:n]; par = par[:n]; M['flags'] = (bag > 0).astype(np.int32); np.savez('$O/meta_bags.npz', **M)
np.savetxt('$O/bags.txt', bag, fmt='%d'); np.savetxt('$O/parents.txt', par, fmt='%d'); print('bagged particles: %d in %d sub-classes of %d parents' % ((bag > 0).sum(), len(set(bag[bag > 0])), len(set(par[bag > 0]))))
PYEOF
for BETA in 0 2; do
  log "board, bag prior beta $BETA"; $PY $W/klem.py $W/simarms/all.mrcs none none $O/board_beta$BETA --cons $W/maps160/synth_arms/consensus_synth.mrc --states 1,2,3,4 --axis=0,0,0,0 --init random --meta $O/meta_bags.npz --bags $O/bags.txt --bag_parents $O/parents.txt --bag_beta $BETA --n1 3 --n2 3 --nfree 0 --focus_mask $R/arm1_travel.mrc --lps "12 10 10 8 8 8" --sigmas "1.0 0.8 0.8 0.6 0.5 0.5" --nproc 16 --mskdiam 230 > $O/board_beta$BETA.log 2>&1 || { log "board FAILED: $(grep -a -m2 -i "error" $O/board_beta$BETA.log | cut -c1-160)"; tail -3 $O/board_beta$BETA.log; continue; }
  $PY - <<PYEOF
import numpy as np
L = np.loadtxt('$O/board_beta$BETA/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int); am = L[:, 3].astype(int); N = int(open('$W/simarms/gt_particles.txt').read().count('\n')) + 2
lab = np.zeros(N, int); lab[pind - 1] = am; np.savetxt('$O/cells_beta$BETA.txt', lab, fmt='%d')
PYEOF
  echo "== arm readout, beta $BETA (per cell std ratio of each arm; 1 = no information)"; $PY $W/label_truth_readout.py $O/cells_beta$BETA.txt $W/simarms/gt_particles.txt $W/truth_arms.txt arms 2>&1 | grep -v Warning | tail -11
done
log DONE
