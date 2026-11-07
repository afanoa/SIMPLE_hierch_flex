#!/bin/bash
# collapse-guard evidence on fulltoy_cdisc's root site-EM rounds: site_size.py hold-out curve per round (OUT <run>/size), class-2 occupancy,
# Spearman with the round it refines; truth AUC from the earlier logs. Run via run_py.sbatch.
K=/mnt/beegfs/elmlund/afan/rnasep_klem; T=/mnt/beegfs/elmlund/afan/routes/toy55; PY=/mnt/nasapps/production/python/3.12/bin/python3; O=$K/factorial/port_sites_cdisc; H=$T/simsites
size(){ d=$1; st=$2; mk=$3; [ -s $d/size/curve.txt ] || $PY $T/site_size.py --em $d/em --start $st --sim ${H}_A,${H}_B --poses $O/node_r/poses_A.txt,$O/node_r/poses_B.txt \
  --mr ${H}_A/meta_rows.txt,${H}_B/meta_rows.txt --nmaster 199996 --mask $mk --out $d/size --cache $T/benchcache_sites --ctfabs 0 --em_lim 15 > $d/size.log 2>&1; }
M=$O/maskref; C=$O/cunion; C2=$O/cunion2
size $M/round1 $M/start_A.txt,$M/start_B.txt $M/round1/mask.mrc & size $M/round2 $M/start_A.txt,$M/start_B.txt $M/round2/mask.mrc &
size $C2/r2mask_r1seed $C2/r2mask_r1seed/start_A.txt,$C2/r2mask_r1seed/start_B.txt $C2/r2mask_r1seed/mask.mrc & wait
size $C2/r1mask_r1seed $C2/r1mask_r1seed/start_A.txt,$C2/r1mask_r1seed/start_B.txt $C2/r1mask_r1seed/mask.mrc & mkdir -p $C/u $C/t; ln -sfn $C/em_union $C/u/em; ln -sfn $C/em_true $C/t/em
size $C/u $C/start_A.txt,$C/start_B.txt $C/union_C.mrc & size $C/t $C/start_A.txt,$C/start_B.txt $T/maps160/simsites/mask_C.mrc & wait
$PY - <<PYEOF
import numpy as np
O = '$O'; R = {'R0 port candidate (10 it)': O + '/fits/r1_s4@r/hi', 'R1 refined round 1': O + '/maskref/round1/size', 'R2 refined round 2': O + '/maskref/round2/size',
     'R2b round-2 mask, round-1 seed': O + '/cunion2/r2mask_r1seed/size', 'R1b round-1 mask, round-1 seed': O + '/cunion2/r1mask_r1seed/size', 'union (truth-chosen)': O + '/cunion/u/size', 'true C mask': O + '/cunion/t/size'}
prev = {'R1 refined round 1': 'R0 port candidate (10 it)', 'R2 refined round 2': 'R1 refined round 1', 'R2b round-2 mask, round-1 seed': 'R1 refined round 1', 'R1b round-1 mask, round-1 seed': 'R1 refined round 1',
        'union (truth-chosen)': 'R0 port candidate (10 it)', 'true C mask': 'R0 port candidate (10 it)'}
P = {}; rk = lambda x: np.argsort(np.argsort(x))
for k, d in R.items():
    try: c = np.loadtxt(d + '/curve.txt'); P[k] = np.loadtxt(d + '/post_master.txt')
    except Exception as e: print(k, 'missing', e); continue
    j = int(np.argmax(c[:, 2])); m = P[k] >= 0
    s = '%-32s hold-out peak %.3f at N %6d, equivalent pure members max %6.0f, class-2 share %.3f' % (k, c[j, 2], c[j, 0], c[:, 6].max(), (P[k][m] >= 0.5).mean())
    if k in prev and prev[k] in P: q = P[prev[k]]; mm = m & (q >= 0); s += ', Spearman with %s %.3f' % (prev[k].split(' ')[0], np.corrcoef(rk(P[k][mm]), rk(q[mm]))[0, 1])
    print(s)
PYEOF
