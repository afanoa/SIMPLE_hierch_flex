#!/bin/bash
#SBATCH -J truth_seed
#SBATCH -p csbdevel
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 00:10:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# synthetic seed = a final-pass labelling (env SEED: children_{A,B}.txt) + the TRUE particles of conformer CONF (18) as state K+1 (oracle seed); env OUT
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; mkdir -p ${OUT:?}
$PY - <<PYEOF
import numpy as np
for h in 'AB':
    s = np.loadtxt('${SEED:?}/children_%s.txt' % h).astype(int); g = np.loadtxt('$W/sim55_700k_%s/gt_particles.txt' % h, comments='#')[:, 1].astype(int)[:len(s)]
    m = (g == ${CONF:-18}) & (s > 0); k = s.max(); s[m] = k + 1; np.savetxt('${OUT}/children_%s.txt' % h, s, fmt='%d')
    print('half %s: conformer ${CONF:-18} -> seed state %d (%d particles); sizes %s' % (h, k + 1, m.sum(), np.bincount(s)[1:].tolist()))
PYEOF
