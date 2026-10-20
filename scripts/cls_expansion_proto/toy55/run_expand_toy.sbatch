#!/bin/bash
#SBATCH -J toy10_expand
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Genetic expansion on toy10: parent run -> families (spectral bisection of the state graph) -> one child klem run per
# family on its particles (finer band, focus recomputed, Fortran kernels) -> children survive by cross-half
# reproducibility -> merged coverage vs the parent. env: PARENT (lat/klem_junk), NFAM (2), LPS, SIGMAS, FOCUS (1), INIT (pca),
# TABLE (pose table for the children, default consensus), PMARG/DEVW/PNOISE (child recipe knobs passed through),
# SPLITMODE (spectral | anneal: the deterministic-annealing criterion decides the number of families, TAU its threshold),
# AUTOGRID=1 sizes every child board to its family (PERCELL particles per cell, 3x3 to 7x7)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; cd $W; PARENT=${PARENT:-$W/lat/klem_junk}; NFAM=${NFAM:-2}; ETAG=${ETAG:-}; E=$W/lat/expand$ETAG; mkdir -p $E
if [ "${SPLITMODE:-spectral}" = anneal ]; then $PY $W/anneal_split.py $PARENT $E --tau ${TAU:-1} --gt $W/sim/gt_particles.txt; NFAM=$(ls $E/family*_subset.txt | wc -l); else $PY $W/expand_families.py split $PARENT $E --nfam $NFAM; fi
CH=""
for j in $(seq 1 $NFAM); do
  echo "[$(date '+%F %T')] child run for family $j"
  NG=""; if [ "${AUTOGRID:-0}" = 1 ]; then NN=$($PY -c "import math; n=$(grep -c . $E/family${j}_subset.txt); print(max(3, min(7, round(math.sqrt(n / ${PERCELL:-1500})))))"); NG="N1=$NN N2=$NN"; echo "  family $j: $(grep -c . $E/family${j}_subset.txt) particles -> ${NN}x${NN} board (~${PERCELL:-1500} particles per cell)"; fi
  env $NG KTAG=_expand${ETAG}_fam$j INIT=${INIT:-pca} NSUB=0 SUBSET=$E/family${j}_subset.txt NFREE=${NFREE:-0} ${FOCUS:+FOCUS=$FOCUS} ${TABLE:+TABLE=$TABLE} ${PMARG:+PMARG=$PMARG} ${DEVW:+DEVW=$DEVW} ${PNOISE:+PNOISE=$PNOISE} LPS="${LPS:-14 14 12 12 10 10 8 8 8}" SIGMAS="${SIGMAS:-1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7}" bash $W/run_klem.sbatch 2>&1 | grep "subset\|focus mask\|final:\|states: cross-half\|Error\|Traceback" | sed "s/^/[fam $j] /"
  CH="$CH $W/lat/klem_expand${ETAG}_fam$j"
done
$PY $W/expand_families.py merge $E $CH --gt $W/sim/gt_particles.txt --minrep ${MINREP:-0.5}
echo "--- parent for reference:"; $PY $W/state_coverage.py $PARENT/klem_latent.txt $W/sim/gt_particles.txt --argmax --tag "parent $(basename $PARENT)" | tail -1
echo "[$(date '+%F %T')] expand DONE"
