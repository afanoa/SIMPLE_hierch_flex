#!/bin/bash
#SBATCH -J closed_seeds
#SBATCH -p csbdevel
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 00:10:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# polish seeds = a final-pass labelling (env SEED dir: children_{A,B}.txt) + the f692e closed group(s) (env CL labels, default 64)
# in f692e cand_{A,B}.txt as extra states (K+1, ...) overriding their seed label; env OUT dir -> states_{A,B}.txt
PY=/mnt/nasapps/production/python/3.12/bin/python3; K=/mnt/beegfs/elmlund/afan/rnasep_klem; C=$K/dtree_rnasep692_f692e; mkdir -p ${OUT:?}
$PY - <<PYEOF
import numpy as np
for h in 'AB':
    s = np.loadtxt('${SEED:?}/children_%s.txt' % h).astype(int); c = np.loadtxt('$C/cand_%s.txt' % h).astype(int); k = s.max()
    for i, l in enumerate([int(v) for v in '${CL:-64}'.split(':')]):
        m = c == l; s[m] = k + i + 1; print('half %s: f692e label %d -> state %d (%d particles)' % (h, l, k + i + 1, m.sum()))
    np.savetxt('${OUT}/states_%s.txt' % h, s, fmt='%d'); print('half %s: sizes %s' % (h, np.bincount(s)[1:].tolist()))
PYEOF
