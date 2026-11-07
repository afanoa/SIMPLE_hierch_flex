#!/bin/bash
# full692 factorial_v3 (port): arm readout (closed / open arm occupancy) of every site pole map + each pole's overlap with the f692st closed
# leaf (rflex_f692st node r1 state 5, 1.07 / 0.51); run via run_py.sbatch. Outputs rnasep_klem/factorial/full692_sites_check/.
K=/mnt/beegfs/elmlund/afan/rnasep_klem; T=/mnt/beegfs/elmlund/afan/routes/toy55; PY=/mnt/nasapps/production/python/3.12/bin/python3
F=$K/stree_port/full692/1_hetero_tree/factorial_v3; O=$K/factorial/full692_sites_check; V=$K/view_rflex/full692_factorial_v3; mkdir -p $O $V
for m in $F/maps/*_pole*.mrc; do n=$(basename $m .mrc); echo "$n $m $m"; done > $O/spec.txt
$PY $T/make_view.py $V $O/spec.txt --box 160 > /dev/null 2>&1
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1); CONS=$T/maps160/toy_consensus.mrc
$PY $T/arm_readout.py $CONS $CL $O/arm.txt $(awk -v v=$V '{print v "/" $1 ".mrc"}' $O/spec.txt) 2>&1 | grep -v Warning | tail -1
$PY - <<PYEOF
import numpy as np, os
F, O = '$F', '$O'; cl = np.loadtxt('$K/rflex_f692st/node_r1/refined_f.txt', dtype=int) == 5; print('reference closed leaf (f692st r1_s5): %d particles' % cl.sum())
arm = {os.path.basename(l.split()[0])[:-4]: (float(l.split()[1]), float(l.split()[2])) for l in open(O + '/arm.txt') if not l.startswith('#')}
rows = []; lab = {}
for n, (c, o) in arm.items():
    site, pole = n.rsplit('_pole', 1); j = int(pole)
    if site not in lab: lab[site] = np.loadtxt(F + '/labels/%s.txt' % site, dtype=int)
    m = lab[site] == j; k = m.sum(); rows.append((c, o, n, k, (m & cl).sum() / max(k, 1), (m & cl).sum() / cl.sum()))
print('%d pole maps read; %d of them read closed-arm >= 0.90' % (len(rows), sum(r[0] >= 0.9 for r in rows)))
print('TOP 12 by closed-arm occupancy: pole | closed-arm open-arm | labelled | share in the closed leaf | closed-leaf recall')
for r in sorted(rows, reverse=True)[:12]: print('  %-18s %.2f / %.2f | %7d | %.3f | %.3f' % (r[2], r[0], r[1], r[3], r[4], r[5]))
print('TOP 8 by share of the closed leaf among the pole\'s labels:')
for r in sorted(rows, key=lambda r: -r[4])[:8]: print('  %-18s %.2f / %.2f | %7d | %.3f | %.3f' % (r[2], r[0], r[1], r[3], r[4], r[5]))
PYEOF
