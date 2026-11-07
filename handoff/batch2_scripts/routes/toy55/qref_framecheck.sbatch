#!/bin/bash
#SBATCH -p norm
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/scorecard/%x.%j.out
# base80's root: does n.simple carry a corr with spread (quintile medians inside dfx deciles), and how far are its poses / corr from the
# node-rows-only stage-B b.simple (mapped through rows.txt)?
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_stree/snap/m8c52t; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local SIMPLE_EMAIL=nil OMP_NUM_THREADS=1
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; H=/mnt/beegfs/elmlund/afan/rnasep_klem/stree_port/full692_base80_s1/1_hetero_tree/node_r
W=/mnt/beegfs/elmlund/afan/rnasep_klem/scorecard/quality_ref/base80_framecheck; mkdir -p $W; cd $W
simple_exec prg=print_project_field projfile=$H/n.simple oritype=ptcl3D 2>/dev/null | $PY $T/dump_ptcl.py n.npz > n.log 2>&1 &
simple_exec prg=print_project_field projfile=$H/stageB/b.simple oritype=ptcl3D 2>/dev/null | $PY $T/dump_ptcl.py b.npz > b.log 2>&1 &
wait; cat n.log b.log
$PY - $H/rows.txt <<'PYEOF'
import sys, numpy as np
N, Bz = np.load('n.npz'), np.load('b.npz'); r = np.loadtxt(sys.argv[1], dtype=int) - 1
def qm(D, rows):
    st, c, f = np.nan_to_num(D['state'][rows]), D['corr'][rows], D['dfx'][rows]; node = (st > 0) & np.isfinite(c) & np.isfinite(f); out = np.zeros(len(st), int)
    dd = np.searchsorted(np.quantile(f[node], np.linspace(0.1, 0.9, 9)), f)
    for j in range(10):
        m = np.nonzero(node & (dd == j))[0]; k = np.argsort(np.argsort(c[m])) / max(len(m) - 1, 1); out[m[k >= 0.8]] = 1; out[m[k <= 0.2]] = 2
    return node.sum(), np.percentile(c[node], [5, 50, 95]), np.median(c[out == 1]), np.median(c[out == 2])
rn = r if len(N['dfx']) > len(r) else np.arange(len(N['dfx']))
for nm, D, rows in (('n.simple', N, rn), ('stageB b.simple', Bz, np.arange(len(Bz['dfx'])))):
    n, p, hi, lo = qm(D, rows); print('%-16s %d rows (%d node): corr 5/50/95 %% %.4f / %.4f / %.4f; quintile medians %.4f / %.4f' % (nm, len(D['dfx']), n, *p, hi, lo))
x, y = N['corr'][rn], Bz['corr']; c = np.isfinite(x) & np.isfinite(y)
print('n vs b at mapped rows: dfx max |diff| %.2e; corr r %.3f; median |de1| %.1f, |de2| %.1f deg' % (np.nanmax(np.abs(N['dfx'][rn] - Bz['dfx'])), np.corrcoef(x[c], y[c])[0, 1],
      np.nanmedian(np.abs(((N['e1'][rn] - Bz['e1']) + 180) % 360 - 180)), np.nanmedian(np.abs(N['e2'][rn] - Bz['e2']))))
PYEOF
