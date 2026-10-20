#!/bin/bash
#SBATCH -J rnasep_tree_realign
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 24:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# One-pass realignment of the real-data tree (run_tree_realign_toy.sbatch on RNase P). The two family pipelines have
# refined every on-grid particle once (10 states per family, local 20->8 A); their poses are merged into ONE meta file
# (family members from their family's refined table, the rest consensus) and the tree is regrown on it from the parent
# down: parent board (junk comps, deviation weights, pose marginalisation) -> families -> child boards at the realigned
# poses -> identification vs the deposited conformers at every level, against the consensus-pose tree (expand1).
# env: E0 (expansion dir whose fam*_pipeline dirs hold refined_ptcl3D.txt), RTAG (_realign1), NFAM (2),
#      parent recipe: PLPS, PSIGMAS, PDEVW (0.05), PPMARG (3), PPNOISE (1), PFOCUS (unset), PBAL (0.5); BAL (0.5) and CPMARG (3; 0 = off) for the children
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; cd $W; log(){ echo "[$(date '+%F %T') rtr] $*"; }
E0=${E0:-$A/rnasep_klem/expand1}; RTAG=${RTAG:-_realign1}; NFAM=${NFAM:-2}; R=$A/rnasep_klem/realign$RTAG; mkdir -p $R
META0=$A/rnasep_coupled/basis_span_all/meta_all.npz; MR=$R/meta_realigned.npz
# 1. one meta file of realigned poses
$PY - <<PYEOF
import numpy as np
M = dict(np.load('$META0')); keys = ('e1', 'e2', 'e3', 'x', 'y'); new = {k: M[k].copy() for k in keys}; src = np.zeros(len(M['e1']), int)
for j in range(1, $NFAM + 1):
    n = 0
    for i, l in enumerate(open('$E0/fam%d_pipeline/refined_ptcl3D.txt' % j)):
        kv = dict(t.split('=', 1) for t in l.split() if '=' in t)
        if float(kv.get('state', 0)) <= 0: continue
        for k in keys: new[k][i] = float(kv[k])
        src[i] = j; n += 1
    print('family %d: %d refined poses' % (j, n))
d = np.degrees(np.arccos(np.clip(np.cos(np.radians(new['e2'] - M['e2'])), -1, 1)))
print('realigned meta: %d particles, %d from the family refinements (|d e2| median %.2f deg), %d keep the consensus pose -> $MR' % (len(src), (src > 0).sum(), np.median(d[src > 0]), (src == 0).sum()))
for k in keys: M[k] = new[k]
np.savez('$MR', **M)
PYEOF
# 2. the parent regrown at the realigned poses
log "parent board at the realigned poses"
env KTAG=$RTAG INIT=npz NFREE=2 META=$MR LPS="${PLPS:-20 20 14 14 10 10 8 8 8}" SIGMAS="${PSIGMAS:-1.5 1.5 1.0 1.0 0.8 0.7 0.7 0.7 0.7}" DEVW=${PDEVW:-0.05} PMARG=${PPMARG:-3} PNOISE=${PPNOISE:-1} ${PFOCUS:+FOCUS=$PFOCUS} BAL=${PBAL:-0.5} bash $W/run_klem_rnasep.sbatch 2>&1 | grep "final:\|states: cross-half\|free1\|Error\|Traceback" | sed 's/^/[parent] /'
echo "== identification of the realigned parent"; $PY $W/deposited_match.py $A/rnasep_klem/klem$RTAG/cells $W/maps160/docked $W/maps160/toy_consensus.mrc | sed 's/^/[parent] /'
# 3. families and child boards at the realigned poses, identification per family
log "families and child boards at the realigned poses"
env PARENT=$A/rnasep_klem/klem$RTAG ETAG=$RTAG NFAM=$NFAM META=$MR BAL=${BAL:-0.5} PMARG=${CPMARG:-3} bash $W/run_expand_rnasep.sbatch 2>&1 | sed 's/^/[expand] /'
$PY - <<PYEOF
import numpy as np
new = [set(np.loadtxt('$A/rnasep_klem/expand$RTAG/family%d_subset.txt' % (j + 1), dtype=int)) for j in range($NFAM)]; old = [set(np.loadtxt('$E0/family%d_subset.txt' % (j + 1), dtype=int)) for j in range($NFAM)]
print('family overlap new x old (Jaccard): ' + ' '.join('n%d/o%d:%.2f' % (a + 1, b + 1, len(new[a] & old[b]) / len(new[a] | old[b])) for a in range($NFAM) for b in range($NFAM)))
PYEOF
log "DONE (second pass per family: env E=$A/rnasep_klem/expand$RTAG FAM=j K=10 META=$MR sbatch run_family_pipeline_rnasep.sbatch)"
