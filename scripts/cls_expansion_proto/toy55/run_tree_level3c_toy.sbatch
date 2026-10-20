#!/bin/bash
#SBATCH -J toy10_tree_level3c
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# SOFT level 3 (hierarchical-GTM style): every child board is bisected (spectral), but instead of hard sub-family
# subsets each particle gets a weight per sub-family = its responsibility mass on that sub-family's cells (x the
# child's own weight in a chained soft tree); the grandchildren train on the union above --thr with these weights in
# the M-step, every particle gets a posterior in both, and the merge sends it to the grandchild where weight x max
# responsibility is largest. env: CB (child board prefix), TABLE (poses), BAL, MSKD, L3 (tag), CWPREFIX (child weights
# prefix for a chained soft tree, e.g. lat/expand_soft/family), THR (0.05)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; cd $W; log(){ echo "[$(date '+%F %T') l3c] $*"; }
NFAM=${NFAM:-2}; CB=${CB:-$W/lat/klem_expand_realign1_bal05_fam}; TABLE=${TABLE:-$W/lat/realign_realign1/ptcl3D_realigned.txt}; L3=${L3:-_l3c}; GC=""; WL=""; CH=""
for j in $(seq 1 $NFAM); do
  E=$W/lat/expand${L3}_fam$j; mkdir -p $E
  $PY $W/expand_families.py softsplit $CB$j $E --nfam 2 --thr ${THR:-0.05} ${CWPREFIX:+--parent_weights ${CWPREFIX}${j}_weights.txt} | sed "s/^/[fam $j softsplit] /"
  for s in 1 2; do
    log "grandchild board (soft): family $j sub-family $s ($(grep -c . $E/family${s}_subset.txt) particles)"
    env KTAG=${L3}_fam${j}_$s INIT=pca NSUB=0 SUBSET=$E/family${s}_subset.txt WEIGHTS=$E/family${s}_weights.txt NFREE=0 TABLE=$TABLE FOCUS=${FOCUS:-1} PMARG=${PMARG:-3} DEVW=0.05 PNOISE=1 ${BAL:+BAL=$BAL} ${MSKD:+MSKD=$MSKD} LPS="${LPS:-12 12 10 10 8 8 8 8 8}" SIGMAS="${SIGMAS:-1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7}" bash $W/run_klem2.sbatch 2>&1 | grep "subset\|soft-family\|final:\|states: cross-half\|Error\|Traceback" | sed "s/^/[fam $j.$s] /"
    GC="$GC $W/lat/klem${L3}_fam${j}_$s"; WL="$WL,$E/family${s}_weights.txt"
  done
  CH="$CH $CB$j"
done
$PY $W/expand_families.py merge $W/lat/expand${L3}_fam1 $GC --gt $W/sim/gt_particles.txt --minrep 0.5 --weights ${WL#,} | sed 's/^/[soft level 3] /'
echo "--- the given children (level 2) for reference:"; $PY $W/expand_families.py merge $W/lat/expand${L3}_fam1 $CH --gt $W/sim/gt_particles.txt --minrep 0.5 ${CWPREFIX:+--weights ${CWPREFIX}1_weights.txt,${CWPREFIX}2_weights.txt} | grep MERGED | sed 's/^/[level 2] /'
log DONE
