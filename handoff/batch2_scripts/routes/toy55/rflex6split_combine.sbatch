#!/bin/bash
#SBATCH -J klem_rflex6comb
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# combine an rflex6split run whose subtrees ran as separate runs (RT_<node>): the main run's leaves.txt (root leaves + the subtrees it
# processed itself) and every RT_<node>/leaves.txt -> RT/combined: leaves.txt, per-row leaf labels (all / non-junk), junk rows, junk flag vs
# truth, scorecards, per-leaf truth report, the bottom-up tree over the leaves. env RT, TAG, CODES (swing:200:1)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem
R=$K_/rflex_${RT:?}; TAG=${TAG:?}; GT=$T/sim$TAG/gt_particles.txt; O=$R/combined; mkdir -p $O; log(){ echo "[$(date '+%F %T') rflex6split-combine $RT] $*"; }
for d in $R $(ls -d ${R}_* 2>/dev/null); do grep -q "rflex6split .*\] DONE" $(ls -t $K_/klem_rflex6split.*.out | xargs grep -l "rflex6split $(basename $d | sed 's/^rflex_//')\]" | head -1) 2>/dev/null || { log "$d not finished"; exit 1; }; done
cat $R/leaves.txt $(ls ${R}_*/leaves.txt 2>/dev/null) > $O/leaves.txt; log "$(grep -c . $O/leaves.txt) leaves from $R and $(ls -d ${R}_* 2>/dev/null | wc -l) subtree runs"
NR=$(wc -l < $R/rows_r.txt)
$PY - <<PYEOF
import numpy as np
L = [l.split() for l in open('$O/leaves.txt') if l.strip()]; out = np.zeros($NR, int); jk = np.zeros($NR, int)
for i, l in enumerate(L, 1):
    lab = np.loadtxt(l[5], dtype=int); m = lab == int(l[6]); assert not (out[m] > 0).any(), 'leaf %s overlaps an earlier leaf' % l[0]
    out[m] = i; jk[m] = int(len(l) > 7 and l[7] == 'junk')
np.savetxt('$O/final_labels.txt', out, fmt='%d'); np.savetxt('$O/final_labels_nojunk.txt', out * (1 - jk), fmt='%d'); np.savetxt('$O/junk_rows.txt', jk, fmt='%d')
open('$O/leaves_list.txt', 'w').write(''.join('%s %s %s %s\n' % (l[0], l[2], l[3], l[4]) for l in L))
print('%d leaves (%d junk) over %d of %d particles: %s' % (len(L), sum(len(l) > 7 for l in L), (out > 0).sum(), $NR, ' '.join('%s:%s%s' % (l[0], l[4], '(junk)' if len(l) > 7 else '') for l in L)))
G = np.loadtxt('$GT', comments='#'); conf = np.zeros($NR, int); conf[G[:, 0].astype(int) - 1] = G[:, 1].astype(int); tj = conf == 0
print('junk flag vs truth: %d of %d true junk particles flagged (recall %.3f); flagged rows %d, of them true junk %.3f' % ((tj & (jk > 0)).sum(), tj.sum(), (tj & (jk > 0)).sum() / max(tj.sum(), 1), jk.sum(), (tj & (jk > 0)).sum() / max(jk.sum(), 1)))
PYEOF
$PY $T/scorecard.py $O/final_labels.txt $GT --codes ${CODES:-swing:200:1} --tag "rflex6split leaves" 2>&1 | grep -v Warning
$PY $T/scorecard.py $O/final_labels_nojunk.txt $GT --codes ${CODES:-swing:200:1} --tag "rflex6split non-junk leaves" 2>&1 | grep -v Warning
$PY $T/leaf_report.py $O $GT 2>&1 | grep -v Warning
$PY $T/state_merge.py $O/final_tree --list $O/leaves_list.txt --link complete --ncut 8,4,2 2>&1 | grep -v Warning | grep "states,\|^== " | head -4
log DONE
