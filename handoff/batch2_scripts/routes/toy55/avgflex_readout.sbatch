#!/bin/bash
#SBATCH -J synth_avgflex_readout
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
set -uo pipefail; PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; MODE=${MODE:?}; D=$W/binsroute_$MODE; O=$D/avgflex_k6; GT=$W/sim$MODE/gt_particles.txt
echo "== $MODE: view-centred flex latent of the bin sub-class averages"; $PY $W/avgflex_readout.py $O $D/p $GT $W/truth_$MODE.txt $MODE --k 6 --q ${Q:-6} 2>&1 | grep -v Warning
NP=$(grep -vc '^#' $GT); $PY - <<PYEOF
import numpy as np; pl = np.loadtxt('$O/particle_states_kmeans.txt', dtype=int); out = np.zeros($NP, int); out[:min(len(pl), $NP)] = pl[:$NP]; np.savetxt('$O/particle_states_kmeans.txt', out, fmt='%d')
PYEOF
echo "== truth readout of the k-means states"; $PY $W/label_truth_readout.py $O/particle_states_kmeans.txt $GT $W/truth_$MODE.txt $MODE 2>&1 | grep -v Warning | tail -8; echo "READOUT DONE $MODE"
