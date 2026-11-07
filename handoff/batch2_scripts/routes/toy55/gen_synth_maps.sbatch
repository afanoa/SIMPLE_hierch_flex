#!/bin/bash
#SBATCH -J synth_maps
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55
for m in arms tethered; do $PY $W/make_synth_maps.py $m $W/maps160/synth_$m $W/synth_$m.csv $W/truth_$m.txt $W/maps160/toy_consensus.mrc --nconf 100; head -3 $W/truth_$m.txt; done
$PY - <<PY
import sys, numpy as np; sys.path.insert(0, '$W'); from model_class import read_mrc
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax = plt.subplots(2, 5, figsize=(20, 8))
for i, m in enumerate(('arms', 'tethered')):
    for j, k in enumerate((1, 2, 3, 4)):
        V = read_mrc('$W/maps160/synth_%s/conf%02d_synth.mrc' % (m, k)); ax[i, j].imshow(V.sum(0), cmap='gray'); ax[i, j].set_title('%s conf %d (z-projection)' % (m, k)); ax[i, j].axis('off')
    V = read_mrc('$W/maps160/synth_%s/consensus_synth.mrc' % m); ax[i, 4].imshow(V.sum(0), cmap='gray'); ax[i, 4].set_title('%s consensus (40 confs)' % m); ax[i, 4].axis('off')
fig.savefig('$W/figs/synth_maps_preview.png', dpi=90, facecolor='white'); print('preview written')
PY
echo "ALL DONE"
