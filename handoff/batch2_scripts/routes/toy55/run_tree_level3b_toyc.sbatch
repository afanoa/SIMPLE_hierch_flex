#!/bin/bash
#SBATCH -J toy10_tree_level3b
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Level 3 on a given pair of child boards that share ONE pose table (the realigned tree): every child board is bisected
# by the spectral split, a grandchild board per sub-family at the shared poses (focus 1, pose marginalisation, 12->8 A,
# optional balanced E-step), the four grandchildren merged and scored against the children's own merge.
# env: CB (child board prefix, e.g. lat/klem_expand_realign1_bal05_fam), TABLE (pose table), BAL (balance lambda), L3 (tag)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; cd $W; log(){ echo "[$(date '+%F %T') l3b] $*"; }
NFAM=${NFAM:-2}; CB=${CB:-$W/lat${TAG:-}/klem_expand_realign1_bal05_fam}; TABLE=${TABLE:-$W/lat${TAG:-}/realign_realign1/ptcl3D_realigned.txt}; L3=${L3:-_l3b}; GC=""; CH=""
for j in $(seq 1 $NFAM); do
  E=$W/lat${TAG:-}/expand${L3}_fam$j; mkdir -p $E
  $PY $W/expand_families.py split $CB$j $E --nfam 2 | sed "s/^/[fam $j split] /"
  $PY $W/state_coverage.py $CB$j/klem_latent.txt $W/${SIMDIR:-sim}/gt_particles.txt --argmax --tag "child board family $j" | tail -1
  for s in 1 2; do
    log "grandchild board: family $j sub-family $s"
    env KTAG=${L3}_fam${j}_$s INIT=pca NSUB=0 SUBSET=$E/family${s}_subset.txt NFREE=0 TABLE=$TABLE FOCUS=${FOCUS:-1} PMARG=${PMARG:-3} DEVW=0.05 PNOISE=1 ${BAL:+BAL=$BAL} ${MSKD:+MSKD=$MSKD} LPS="${LPS:-12 12 10 10 8 8 8 8 8}" SIGMAS="${SIGMAS:-1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7}" bash $W/run_klem.sbatch 2>&1 | grep "subset\|final:\|states: cross-half\|Error\|Traceback" | sed "s/^/[fam $j.$s] /"
    GC="$GC $W/lat${TAG:-}/klem${L3}_fam${j}_$s"; $PY $W/state_coverage.py $W/lat${TAG:-}/klem${L3}_fam${j}_$s/klem_latent.txt $W/${SIMDIR:-sim}/gt_particles.txt --argmax --tag "grandchild $j.$s" | tail -1
  done
  CH="$CH $CB$j"
done
$PY $W/expand_families.py merge $W/lat${TAG:-}/expand${L3}_fam1 $GC --gt $W/${SIMDIR:-sim}/gt_particles.txt --minrep 0.5 | sed 's/^/[level 3] /'
echo "--- level 2 (the given children) for reference:"; $PY $W/expand_families.py merge $W/lat${TAG:-}/expand${L3}_fam1 $CH --gt $W/${SIMDIR:-sim}/gt_particles.txt --minrep 0.5 | grep MERGED | sed 's/^/[level 2] /'
log DONE
