#!/bin/bash
#SBATCH -J klem_twunits
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# tree-weighted FINAL maps of a finished run's delivered units (tw_units.py): unit half-maps from the merge's leaves -> unit tree (state_merge.py
# --list, complete linkage) -> tree2 / emp2 (and raw / own) posterior-mean maps at box 160, lp 10, unsharpened, halves in halves/ -> closed-arm
# readout of every map set (arm_readout.py) and a per-unit comparison with the current maps (unit_maps.py, MERGE/unit_arm.txt).
# env MERGE (merge dir: leaves_int.txt, final_leaves.txt, final_junk.txt), RT (run tag: maps -> view_rflex/<RT>_merged_tree2, _emp2), MSKDIAM (232)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem
MERGE=${MERGE:?}; RT=${RT:?}; O=$K_/tw_units/${RT}_$(basename $MERGE); V=$K_/view_rflex/${RT}_merged; export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16; mkdir -p $O
log(){ echo "[$(date '+%F %T') tw_units $RT] $*"; }
$PY $T/tw_units.py $O --stage halves --leaves $MERGE/leaves_int.txt --units $MERGE/final_leaves.txt || { log "unit half-maps FAILED"; exit 1; }
[ -s $O/tree/merge_tree.txt ] || $PY $T/state_merge.py $O/tree --list $O/units_list.txt --link complete --ncut 2 > $O/tree.log 2>&1 || { log "unit tree FAILED"; tail -3 $O/tree.log; exit 1; }
$PY $T/tw_units.py $O --stage maps --junk $MERGE/final_junk.txt --mapdir $V --mskdiam ${MSKDIAM:-232} 2>&1 | grep -v Warning
[ -s $O/map_dirs.txt ] || { log "maps FAILED"; exit 1; }
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
while read k d; do $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/arm_$k.txt $(ls $d/*.mrc) > /dev/null 2>&1 || log "arm readout $k FAILED"; done < $O/map_dirs.txt
$PY - <<PYEOF
import re, numpy as np
def arm(fn):
    A = {}
    for l in open(fn):
        if l.startswith('#') or not l.strip(): continue
        w = l.split(); A[int(re.findall(r'unit_(\d+)\.mrc', w[0])[0])] = (float(w[1]), float(w[2]), w[3])
    return A
sets = [('current', '$MERGE/unit_arm.txt')] + [(k, '$O/arm_%s.txt' % k) for k in ('raw', 'own', 'tree2', 'emp2')]; A = {k: arm(f) for k, f in sets}
cnt = {int(l.split()[0]): int(l.split()[3]) for l in open('$O/units_list.txt') if l.strip()}; junk = {int(x) for x in open('$MERGE/final_junk.txt').read().split()}
print('closed-arm / open-arm occupancy per unit (call): ' + ' | '.join(k for k, _ in sets))
for u in sorted(cnt, key=lambda u: -A['current'][u][0]):
    print('   %sunit %2d (%6d ptcls): %s' % ('junk ' if u in junk else '     ', u, cnt[u], ' | '.join('%.2f/%.2f %-12s' % A[k][u] for k, _ in sets)))
real = [u for u in cnt if u not in junk]; c = {k: np.array([A[k][u][0] for u in real]) for k, _ in sets}
print('non-junk units: closed-arm corr with current: %s; spread (sd): %s; calls closed/intermediate/open: %s' % (
      ' '.join('%s %.3f' % (k, np.corrcoef(c['current'], c[k])[0, 1]) for k, _ in sets[1:]), ' '.join('%s %.3f' % (k, c[k].std()) for k, _ in sets),
      ' '.join('%s %d/%d/%d' % (k, *[sum(A[k][u][2] == x for u in real) for x in ('closed', 'intermediate', 'open')]) for k, _ in sets)))
PYEOF
log "DONE: $(cat $O/choice.txt)"
