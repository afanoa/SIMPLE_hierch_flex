#!/bin/bash
#SBATCH -J toy10_tree_soft
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# The full SOFT tree from a parent board: soft families (weights = family responsibility), children on the weighted
# unions at the realigned poses (feedback recipe, balanced), soft merge; then soft level 3 chained on the children's
# weights. env: PARENT (lat/klem_realign1_bal05), TABLE, BAL (0.5), MSKD, STAG (_soft), THR (0.05)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; cd $W; log(){ echo "[$(date '+%F %T') soft] $*"; }
PARENT=${PARENT:-$W/lat/klem_realign1_bal05}; TABLE=${TABLE:-$W/lat/realign_realign1/ptcl3D_realigned.txt}; STAG=${STAG:-_soft}; E=$W/lat/expand$STAG; mkdir -p $E; CH=""; WL=""
$PY $W/expand_families.py softsplit $PARENT $E --nfam 2 --thr ${THR:-0.05} | sed 's/^/[softsplit] /'
for j in 1 2; do
  log "child board (soft) for family $j ($(grep -c . $E/family${j}_subset.txt) particles)"
  env KTAG=_expand${STAG}_fam$j INIT=pca NSUB=0 SUBSET=$E/family${j}_subset.txt WEIGHTS=$E/family${j}_weights.txt NFREE=0 TABLE=$TABLE PMARG=3 DEVW=0.05 PNOISE=1 BAL=${BAL:-0.5} ${MSKD:+MSKD=$MSKD} bash $W/run_klem2.sbatch 2>&1 | grep "subset\|soft-family\|final:\|states: cross-half\|Error\|Traceback" | sed "s/^/[fam $j] /"
  CH="$CH $W/lat/klem_expand${STAG}_fam$j"; WL="$WL,$E/family${j}_weights.txt"
done
$PY $W/expand_families.py merge $E $CH --gt $W/sim/gt_particles.txt --minrep 0.5 --weights ${WL#,} | sed 's/^/[soft level 2] /'
log "soft level 3 chained on the soft children"
env CB=$W/lat/klem_expand${STAG}_fam TABLE=$TABLE BAL=${BAL:-0.5} ${MSKD:+MSKD=$MSKD} L3=${STAG}_l3 CWPREFIX=$E/family THR=${THR:-0.05} bash $W/run_tree_level3c_toy.sbatch 2>&1 | grep "softsplit\]\|soft-family\|final:\|MERGED\|best-state\|Error\|Traceback\|DONE" | sed 's/^/[l3] /'
log DONE
