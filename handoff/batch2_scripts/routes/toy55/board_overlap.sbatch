#!/bin/bash
#SBATCH -J klem_overlap
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# do the region boards of a family find the same groups? for every pair of level-1 states from DIFFERENT region boards of family FAM
# of run RTAG: shared particles / smaller state (before consolidation hands each particle to one board); the closed-leaning
# leaves (s01 of each board, CLOSED list "reg:state ...") reported first
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; R=/mnt/beegfs/elmlund/afan/rnasep_klem/${RTAG:?}
$PY - <<PYEOF
import glob, os, numpy as np
D = '$R/fam${FAM:-1}'; L = {}
for f in sorted(glob.glob(D + '/region*/level1/states_clusters.txt')):
    b = f.split('/')[-3]; L[b] = np.loadtxt(f, dtype=int)
names = list(L); print('family ${FAM:-1}: %d region boards: %s' % (len(names), ' '.join('%s(%d states)' % (b, L[b].max()) for b in names)))
sets = {(b, s): np.where(L[b] == s)[0] for b in names for s in range(1, L[b].max() + 1) if (L[b] == s).sum() >= 1000}
def ov(p, q): return len(np.intersect1d(sets[p], sets[q], assume_unique=True)) / min(len(sets[p]), len(sets[q]))
cl = [tuple([w.split(':')[0], int(w.split(':')[1])]) for w in '${CLOSED:-}'.split()]
if cl:
    print('closed-leaning leaves (shared particles / smaller state):')
    for i, p in enumerate(cl):
        for q in cl[i + 1:]:
            if p in sets and q in sets: print('  %s s%d  vs  %s s%d : %.2f' % (p[0], p[1], q[0], q[1], ov(p, q)))
best = []
for p in sets:
    o = [(ov(p, q), q) for q in sets if q[0] != p[0]]
    if o: m = max(o); best.append(m[0])
best = np.array(best); print('every state vs its best match on another board: median overlap %.2f, %d of %d states >= 0.5, %d >= 0.8' % (np.median(best), (best >= 0.5).sum(), len(best), (best >= 0.8).sum()))
PYEOF
