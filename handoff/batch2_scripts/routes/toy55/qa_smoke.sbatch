#!/bin/bash
#SBATCH -p norm
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:10:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/scorecard/%x.%j.out
# qa_ratio.py on a synthetic one-node tree: bin 5's corr is lower, so L = bin 5 and the fake QS3 line becomes a ratio
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; W=/mnt/beegfs/elmlund/afan/rnasep_klem/scorecard/quality_a/smoke_$SLURM_JOB_ID
mkdir -p $W; cd $W
$PY - <<'PYEOF'
import os, numpy as np
rng = np.random.default_rng(0); n = 5000; os.makedirs('tree/node_r/motion/c1', exist_ok=True); os.makedirs('ref_r', exist_ok=True)
np.savetxt('tree/node_r/rows.txt', np.arange(1, n + 1), fmt='%d'); b = rng.integers(1, 6, n); np.savetxt('tree/node_r/motion/c1/bins.txt', b, fmt='%d')
for x in ['recvol_state%02d%s.mrc' % (k, s) for k in (1, 5) for s in ('', '_even', '_odd', '_even_unfil', '_odd_unfil')] + ['hetero_step.sh']: open('tree/node_r/motion/c1/' + x, 'w').close()
open('tree/node_r/mask.mrc', 'w').close(); open('tree/node_r/motion/c1/motion_call.txt', 'w').write('switch dipole\n')
c = 0.4 + 0.01 * rng.normal(size=n) - 0.004 * (b == 5); f = rng.uniform(1, 3, n)
np.savez('ref_r/dump.npz', state=np.ones(n), corr=c, dfx=f, e1=np.zeros(n), e2=np.zeros(n))
st = np.zeros(n, int); r = np.argsort(np.argsort(c)) / (n - 1.); st[r >= 0.8] = 1; st[r <= 0.2] = 2; np.savetxt('ref_r/states.txt', st, fmt='%d')
open('ref_r/dBref.txt', 'w').write('QS3 x ref | 24 shells | dBx +200.0 A^2 (SE 10.0) | dqx -1 (SE 0.1, z -6)\n')
PYEOF
$PY $T/qa_ratio.py cases $W/tree other r $W $W/ref_ && cat qscreen.txt cases_info.tsv
echo "QS3 ref_r_c1 cand | 24 shells | dBx +30.0 A^2 (SE 3.0) | dqx -0.1 (SE 0.01, z -10)" > qs3.txt
$PY $T/qa_ratio.py table $W/tree other r $W $W/ref_
