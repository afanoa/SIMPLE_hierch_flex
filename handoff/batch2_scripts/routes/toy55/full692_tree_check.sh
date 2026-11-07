#!/bin/bash
# full692 (port) tree leaves: arm readout of final_maps/raw/unit_NN (leaf NN of hetero_tree_leaves.txt) + each leaf's overlap with the f692st
# closed leaf (rflex_f692st node r1 state 5); run via run_py.sbatch. Outputs rnasep_klem/factorial/full692_sites_check/.
K=/mnt/beegfs/elmlund/afan/rnasep_klem; T=/mnt/beegfs/elmlund/afan/routes/toy55; PY=/mnt/nasapps/production/python/3.12/bin/python3
P=$K/stree_port/full692/1_hetero_tree; O=$K/factorial/full692_sites_check; V=$K/view_rflex/full692_tree_units; mkdir -p $O $V
for m in $P/final_maps/raw/unit_*.mrc; do n=$(basename $m .mrc); echo "$n $m $m"; done > $O/spec_units.txt
[ -s $O/arm_units.txt ] || { $PY $T/make_view.py $V $O/spec_units.txt --box 160 > /dev/null 2>&1
  CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1); $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/arm_units.txt $(awk -v v=$V '{print v "/" $1 ".mrc"}' $O/spec_units.txt) 2>&1 | grep -v Warning | tail -1; }
$PY - <<PYEOF
import numpy as np, os
L = np.loadtxt('$P/hetero_tree_labels.txt', dtype=int); cl = np.loadtxt('$K/rflex_f692st/node_r1/refined_f.txt', dtype=int) == 5
lv = [l.split()[0] for l in open('$P/hetero_tree_leaves.txt') if l.strip()]
A = {os.path.basename(l.split()[0])[:-4]: (float(l.split()[1]), float(l.split()[2])) for l in open('$O/arm_units.txt') if not l.startswith('#')}
R = []
for k in range(1, len(lv) + 1):
    m = L == k; c, o = A.get('unit_%02d' % k, (np.nan, np.nan)); R.append((c, o, k, lv[k - 1], m.sum(), (m & cl).sum() / max(m.sum(), 1), (m & cl).sum() / cl.sum()))
print('full692 tree: %d leaves, labels %d; closed leaf of f692st: %d particles' % (len(lv), (L > 0).sum(), cl.sum()))
for r in sorted(R, reverse=True)[:8]: print('  unit %2d %-9s closed-arm %.2f open-arm %.2f | %6d particles | share in the closed leaf %.3f | closed-leaf recall %.3f' % (r[2], r[3], r[0], r[1], r[4], r[5], r[6]))
PYEOF
