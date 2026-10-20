#!/bin/bash
#SBATCH -J klem_synth_truth
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# truth readout of a synthetic benchmark chain (env MODE arms|tethered): per region board the majority configuration,
# purity and size of every leaf state, and the coverage of the 100 configurations across all leaves
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/synth_${MODE:?}; GT=$W/sim$MODE/gt_particles.txt
echo "== $MODE: ground truth $(grep -vc '^#' $GT) rows, $(awk '!/^#/{print $2}' $GT | sort -u | wc -l) configurations (+junk 0)"
echo "== parent / feedback boards (argmax cell vs configuration): NMI"; for b in $R/parent $R/fam?/klem_fb; do [ -s $b/klem_latent.txt ] || continue; $PY - <<PYEOF
import numpy as np
L = np.loadtxt('$b/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int); am = L[:, 3].astype(int)
gt = {}; [gt.__setitem__(int(w[0]), int(w[1])) for w in (l.split() for l in open('$GT')) if w and not w[0].startswith('#')]
g = np.array([gt.get(p, 0) for p in pind]); m = g > 0
from collections import Counter
def H(x): c = np.array(list(Counter(x).values()), float); p = c / c.sum(); return -(p * np.log(p)).sum()
hx, hy = H(am[m]), H(g[m]); hxy = H(list(zip(am[m], g[m]))); print('  %-50s NMI(cell, configuration) %.3f over %d particles' % ('$b'.split('synth_$MODE/')[-1], (hx + hy - hxy) / max(min(hx, hy), 1e-9), m.sum()))
PYEOF
done
cov=""; for f in $R/fam?/region*/level1/states_clusters.txt; do t=$(echo $f | sed "s#$R/##; s#/level1/states_clusters.txt##"); $PY $W/score_states.py $f $GT --tag $t 2>&1 | grep -v Warning | head -1; done
$PY - <<PYEOF
import numpy as np, glob
gt = {}; [gt.__setitem__(int(w[0]), int(w[1])) for w in (l.split() for l in open('$GT')) if w and not w[0].startswith('#')]
cov = set(); pur = []
for f in glob.glob('$R/fam?/region*/level1/states_clusters.txt'):
    st = np.loadtxt(f, dtype=int)
    for s in range(1, st.max() + 1):
        rows = np.where(st == s)[0] + 1; g = np.array([gt.get(int(p), 0) for p in rows]); g = g[g > 0]
        if len(g): v, c = np.unique(g, return_counts=True); cov.add(int(v[c.argmax()])); pur.append(c.max() / len(g))
print('== coverage over all leaf states: %d of 100 configurations are the majority of some leaf; mean leaf purity %.2f (chance ~0.01)' % (len(cov), np.mean(pur) if pur else 0))
PYEOF
echo "SYNTH TRUTH DONE"
