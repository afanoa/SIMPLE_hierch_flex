#!/bin/bash
#SBATCH -J klem_mergeup
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 16:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# merge up the tree of an f692e-route run (RTAG): gold-standard split tests on the real halves, zoom level (each region leaf's
# zoom sub-states) then region level (each region board's surviving units); a unit that does not reproduce across the halves
# merges into the surviving group that takes its particles in the cross-half EM. Then every merged-up state is reconstructed
# from all its particles (mu_rec.sbatch) and read out (pooled_analysis.sbatch). env RTAG, SRCC (consolidated_12_fullw_zoom),
# OUT (mergeup), NCHUNK (8), XEM (3); a synthetic set: SIMP, LATP, CACHE, CTFABS 0, GTX "--gt A,B" (truth NMI in the logs)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; K=$A/rnasep_klem; R=$K/${RTAG:?}; SC=$R/${SRCC:-consolidated_12_fullw_zoom}
O=$R/${OUT:-mergeup}; mkdir -p $O/consolidated; log(){ echo "[$(date '+%F %T') mergeup $RTAG] $*" | tee -a $O/summary.txt; }
SIMP=${SIMP:-$W/simrnasep692}; MRA=${SIMP}_A/meta_rows.txt; MRB=${SIMP}_B/meta_rows.txt; P=$K/dtree_$(basename $SIMP)_$RTAG/poses_famrefined; export SIMP
[ -s ${P}_B.txt ] || { mkdir -p $(dirname $P); $PY $W/f692e_half_poses.py $R $SIMP ${LATP:-$W/lat_rnasep692} $P 2>&1 | grep -v Warning | tail -2; }
level(){ local lv=$1 src=$2 L=$O/$1; [ -s $L/units/units.txt ] && return 0
  $PY $W/mu_prep.py $lv $R $src $MRA $MRB $L 2>&1 | grep -v Warning | tee -a $O/summary.txt
  local J=""; for c in $(seq 0 $((${NCHUNK:-8} - 1))); do J="$J $(sbatch --parsable --export=ALL,LEVEL=$L,POSES=$P,XEM=${XEM:-3},NCHUNK=${NCHUNK:-8},CHUNK=$c $W/mu_tests.sbatch)"; done
  log "$lv level: split-test jobs$J"
  for j in $J; do while squeue -h -j $j 2>/dev/null | grep -q .; do sleep 30; done; grep -v 'CHUNK.*DONE' $K/klem_mutest.$j.out; done | tee -a $O/summary.txt
  $PY $W/mu_collect.py $L $src $MRA $MRB $L/units $SC/arm_states.txt 2>&1 | grep -v Warning | tail -2 | tee -a $O/summary.txt; }
log "merge up from $SC"; level zoom $SC; level region $O/zoom/units
$PY - <<PYEOF 2>&1 | grep -v Warning | tee -a $O/summary.txt
import numpy as np
U = '$O/region/units'; st = np.zeros(len(np.loadtxt('$SC/states.txt', dtype=int)), int)
for h, f in (('A', '$MRA'), ('B', '$MRB')): st[np.loadtxt(f).astype(int) - 1] = np.loadtxt('%s/units_%s.txt' % (U, h), dtype=int)
np.savetxt('$O/consolidated/states.txt', st, fmt='%d'); rows = [l.split() for l in open(U + '/units.txt') if not l.startswith('#')]
with open('$O/consolidated/groups.txt', 'w') as f:
    f.write('# merged-up state: particles, members (consolidated zoom states), closed-arm max over members\n')
    for w in rows: f.write('%s %d %s %s\n' % (w[0], int(w[4]) + int(w[5]), w[6], w[7]))
src = sum(1 for l in open('$SC/groups.txt') if not l.startswith('#')); n = [int(w[4]) + int(w[5]) for w in rows]
print('merged up: %d consolidated zoom states -> %d states (%d particles); closed-arm >= 0.8 among members: %s' % (src, len(rows), sum(n), ' '.join('%s(%s)' % (w[0], w[7]) for w in rows if float(w[7]) >= 0.8) or 'none'))
PYEOF
log "reconstructing the merged-up states"; RTAG=$RTAG CDIR=${OUT:-mergeup}/consolidated bash $W/mu_rec.sbatch 2>&1 | grep -v "^  *[0-9]* / \|Traceback\|^  File\|^ *\^" | tail -5 | tee -a $O/summary.txt
RTAG=$RTAG CDIR=${OUT:-mergeup}/consolidated bash $W/pooled_analysis.sbatch 2>&1 | tee -a $O/summary.txt
bash $W/tidy_halves.sh $O/consolidated; log DONE
