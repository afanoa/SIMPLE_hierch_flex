#!/bin/bash
#SBATCH -J klem_mergeup
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 16:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# family level of merge up (after merge_up.sbatch): every surviving unit of a family is tested against every other across
# branches and region boards, inside the family's periphery (consensus envelope minus the rigid core); a pair merges when its
# difference does not reproduce across the halves beyond a label-shuffle null. No blend pruning and no re-assignment, so
# intermediates and small closed groups stay states; duplicates found by different boards pool. env RTAG, OUT (mergeup), CORE (0.40)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; K=$A/rnasep_klem; R=$K/${RTAG:?}; SC=$R/${SRCC:-consolidated_12_fullw_zoom}
O=$R/${OUT:-mergeup}; L=$O/family; mkdir -p $O/consolidated_family; log(){ echo "[$(date '+%F %T') mergeup-family $RTAG] $*" | tee -a $O/summary.txt; }
SIMP=${SIMP:-$W/simrnasep692}; MRA=${SIMP}_A/meta_rows.txt; MRB=${SIMP}_B/meta_rows.txt; P=$K/dtree_$(basename $SIMP)_$RTAG/poses_famrefined; export SIMP
[ -s $O/region/units/units.txt ] || { log "needs $O/region/units (merge_up.sbatch)"; exit 1; }
for F in 1 2; do [ -s $O/periphery_fam$F.mrc ] || $PY $W/periphery_mask.py ${PARENTS:-$R/parents}/fam$F.mrc $O/periphery_fam$F.mrc --core ${CORE:-0.40} 2>&1 | grep -v Warning | tail -2; done
if [ ! -s $L/units/units.txt ]; then
  FOCUS_FAM=$O/periphery_fam{F}.mrc $PY $W/mu_prep.py family $R $O/region/units $MRA $MRB $L 2>&1 | grep -v Warning | tee -a $O/summary.txt
  J=""; for c in 0 1; do J="$J $(XEM=0 GTX="--arch 0 ${GTXF:-}" sbatch --parsable --export=ALL,LEVEL=$L,POSES=$P,NCHUNK=2,CHUNK=$c,PAR=1 --mem=200G $W/mu_tests.sbatch)"; done
  log "family level: split-test jobs$J (no blend pruning, no re-assignment)"
  for j in $J; do while squeue -h -j $j 2>/dev/null | grep -q .; do sleep 30; done; grep -v 'CHUNK.*DONE' $K/klem_mutest.$j.out; done | tee -a $O/summary.txt
  $PY $W/mu_collect.py $L $O/region/units $MRA $MRB $L/units $SC/arm_states.txt 2>&1 | grep -v Warning | tail -2 | tee -a $O/summary.txt
fi
$PY - <<PYEOF 2>&1 | grep -v Warning | tee -a $O/summary.txt
import numpy as np
U = '$L/units'; st = np.zeros(len(np.loadtxt('$SC/states.txt', dtype=int)), int)
for h, f in (('A', '$MRA'), ('B', '$MRB')): st[np.loadtxt(f).astype(int) - 1] = np.loadtxt('%s/units_%s.txt' % (U, h), dtype=int)
np.savetxt('$O/consolidated_family/states.txt', st, fmt='%d'); rows = [l.split() for l in open(U + '/units.txt') if not l.startswith('#')]
with open('$O/consolidated_family/groups.txt', 'w') as f:
    f.write('# merged-up state (family level): particles, members (consolidated zoom states), closed-arm max over members\n')
    for w in rows: f.write('%s %d %s %s\n' % (w[0], int(w[4]) + int(w[5]), w[6], w[7]))
n = sorted((int(w[4]) + int(w[5]) for w in rows), reverse=True)
print('family level: %d states, particles per state median %d (max %d); closed-arm >= 0.8 among members: %s' % (len(rows), n[len(n) // 2], n[0], ' '.join('%s(%s, %d)' % (w[0], w[7], int(w[4]) + int(w[5])) for w in rows if float(w[7]) >= 0.8) or 'none'))
PYEOF
log "reconstructing the family-level states"; RTAG=$RTAG CDIR=${OUT:-mergeup}/consolidated_family bash $W/mu_rec.sbatch 2>&1 | grep -v "^  *[0-9]* / \|Traceback\|^  File\|^ *\^" | tail -3 | tee -a $O/summary.txt
RTAG=$RTAG CDIR=${OUT:-mergeup}/consolidated_family bash $W/pooled_analysis.sbatch 2>&1 | tee -a $O/summary.txt
bash $W/tidy_halves.sh $O/consolidated_family; log DONE
