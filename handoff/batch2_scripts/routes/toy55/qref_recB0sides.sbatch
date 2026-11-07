#!/bin/bash
#SBATCH -p norm
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:15:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/scorecard/%x.%j.out
# The root replay's recB0/b.simple (node_motion's frame): same poses / corr as the site candidates' motion/cK/m.simple? Each site's
# lower-corr end bin (1 vs 5) by recB0's corr.
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_stree/snap/m8c52t; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local SIMPLE_EMAIL=nil OMP_NUM_THREADS=1
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; G=/mnt/beegfs/elmlund/afan/rnasep_klem/stree_port/gate/m33_root_m79/1_hetero_tree/node_r
W=/mnt/beegfs/elmlund/afan/rnasep_klem/scorecard/quality_ref/m79_recB0_sides; mkdir -p $W; cd $W
simple_exec prg=print_project_field projfile=$G/recB0/b.simple oritype=ptcl3D 2>/dev/null | $PY $T/dump_ptcl.py recB0.npz > recB0.log 2>&1; cat recB0.log
$PY - $G <<'PYEOF'
import sys, numpy as np
G = sys.argv[1]; R = np.load('recB0.npz'); M = np.load('/mnt/beegfs/elmlund/afan/rnasep_klem/scorecard/quality_ref/m79_root_mframe/dump.npz')
c = np.isfinite(R['corr']) & np.isfinite(M['corr'])
print('recB0 vs c6 m.simple: dfx max |diff| %.1e; corr max |diff| %.1e (r %.4f); e1 max |diff| %.2f, e2 max |diff| %.2f deg' % (np.nanmax(np.abs(R['dfx'] - M['dfx'])),
      np.abs(R['corr'][c] - M['corr'][c]).max(), np.corrcoef(R['corr'][c], M['corr'][c])[0, 1], np.nanmax(np.abs(((R['e1'] - M['e1']) + 180) % 360 - 180)), np.nanmax(np.abs(R['e2'] - M['e2']))))
for s in (1, 2, 3, 5, 6, 7):
    b = np.loadtxt(G + '/motion/c%d/bins.txt' % s); m1, m5 = np.nanmedian(R['corr'][b == 1]), np.nanmedian(R['corr'][b == 5])
    print('site %d: recB0 corr median bin1 %.4f, bin5 %.4f -> lower-corr end bin%d (%d / %d particles)' % (s, m1, m5, 1 if m1 < m5 else 5, (b == 1).sum(), (b == 5).sum()))
PYEOF
