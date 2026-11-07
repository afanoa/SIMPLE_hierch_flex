#!/bin/bash
# re-cut every saved junk_test run with junk_recut.py (xd2 vs xd3) + the truth make-up of simsites' converged site-C classes; run via run_py.sbatch
K=/mnt/beegfs/elmlund/afan/rnasep_klem; T=/mnt/beegfs/elmlund/afan/routes/toy55; JT=$K/factorial/junk_test; PY=/mnt/nasapps/production/python/3.12/bin/python3; C=$K/stree_port/conv_r1
R(){ $PY $T/junk_recut.py "$@" 2>&1 | grep -v Warning; }
ST="--truth $JT/simsites_truth_jt.txt --jtype jtype --guard C_closed,A_domain"; IT="--truth $K/factorial/armsigm1_truth/arms_truth.txt --guardtail a1_th,a1_ph,a2_th,a2_ph,a3_th,a3_ph,a4_th,a4_ph,a5_th,a5_ph"
R $JT/simsites/K1/junk_scores.txt --tag simsites_K1 $ST; R $JT/simsites/truth8/junk_scores.txt --tag simsites_truth8 $ST
R $JT/simsites/root8/junk_scores.txt --labels $K/rflex_simsites/node_r/refined_f.txt --tag simsites_root8 $ST
R $JT/simsites_convC/convC/junk_scores.txt --labels $JT/simsites_convC/labels_master.txt --tag simsites_convC $ST --out $JT/simsites_convC/convC/junk_flags_xd3.txt
R $JT/igm/K1/junk_scores.txt --tag igm_K1 $IT; R $JT/igm/root8/junk_scores.txt --labels $K/rflex_armsigm1s/node_r/refined_f.txt --tag igm_root8 $IT
R $JT/r1conv/K1/junk_scores.txt --tag r1conv_K1
for p in 1 2 3; do R $JT/r1conv/c$p/junk_scores.txt --labels $C/root_p$p/1_site_rank/site_seeds.txt --tag r1conv_c$p --out $JT/r1conv/c$p/junk_flags_xd3.txt; done
for u in r r2 r4 r6; do R $JT/f692st_$u/K1/junk_scores.txt --tag f692st_${u}_K1 --out $JT/f692st_$u/K1/junk_flags_xd3.txt
  R $JT/f692st_$u/$u/junk_scores.txt --labels $K/rflex_f692st/node_$u/refined_f.txt --tag f692st_${u}_states --out $JT/f692st_$u/$u/junk_flags_xd3.txt; done
$PY - <<PYEOF
import numpy as np
L = np.loadtxt('$JT/simsites_convC/labels_master.txt', dtype=int); T = np.loadtxt('$JT/simsites_truth_jt.txt'); hd = [l.split()[1:] for l in open('$JT/simsites_truth_jt.txt') if l.startswith('#')][0]
c = lambda n: T[:, hd.index(n)]; jt = c('jtype'); cl = (c('C_closed') == 1) & (jt == 0)
for s in (1, 2):
    m = L == s; print('converged site-C class %d (%d): closed genuine %d, open genuine %d, blob junk %d, broken junk %d' % (s, m.sum(), (m & cl).sum(), (m & ~cl & (jt == 0)).sum(), (m & (jt == 1)).sum(), (m & (jt == 2)).sum()))
PYEOF
