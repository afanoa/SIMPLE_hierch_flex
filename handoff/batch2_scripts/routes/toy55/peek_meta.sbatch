#!/bin/bash
#SBATCH -J peek_meta
#SBATCH -p csbdevel
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 00:10:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan
$PY - <<PYEOF
import numpy as np, glob
M = np.load('$A/rnasep_coupled/basis_span_all/meta_all.npz'); print('keys', list(M.keys()))
f = M['flags'].astype(int); eo = M['eo'].astype(int); print('rows', len(f), 'refined', (f == 1).sum(), 'eo among refined', np.bincount(eo[f == 1]).tolist())
print('stkind range', M['stkind'].min(), M['stkind'].max(), 'indstk range', M['indstk'].min(), M['indstk'].max())
for k in ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast'): v = M[k][f == 1]; print(k, '%.3f %.3f %.3f' % (v.min(), np.median(v), v.max()))
for s in sorted(glob.glob('/mnt/beegfs/elmlund/rnasep/22_denoise_project/denoise_raw_particles_part*.mrcs')):
    h = np.fromfile(s, np.int32, 256); print(s.split('/')[-1], 'nx ny nz mode', h[0], h[1], h[2], h[3], 'nsym', h[23])
PYEOF
