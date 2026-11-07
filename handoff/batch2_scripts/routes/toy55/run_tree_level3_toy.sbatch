#!/bin/bash
#SBATCH -J toy10_tree_level3
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# A third level of the tree, where the children still mix near-identical conformers (48/49, 1/10): every feedback child
# board (family j at its refined poses) is split into two sub-families by the same spectral bisection of its state graph,
# and a grandchild board is run per sub-family at the family's refined poses (focus mask on, no free comps, pose
# marginalisation, fine band 12->8 A). The four grandchildren are merged and scored exactly like the children, against the
# children's own merged purity (0.62). env: FB (lat/fb_fam), CB (lat/klem_fb_fam), NFAM (2), LPS, SIGMAS, FOCUS (1), PMARG (3)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; cd $W; log(){ echo "[$(date '+%F %T') l3] $*"; }
NFAM=${NFAM:-2}; FB=${FB:-$W/lat/fb_fam}; CB=${CB:-$W/lat/klem_fb_fam}; L3=${L3:-_l3}; GC=""; CH=""
for j in $(seq 1 $NFAM); do
  E=$W/lat/expand${L3}_fam$j; mkdir -p $E
  $PY $W/expand_families.py split $CB$j $E --nfam 2 | sed "s/^/[fam $j split] /"
  $PY $W/state_coverage.py $CB$j/klem_latent.txt $W/sim/gt_particles.txt --argmax --tag "child board family $j" | tail -1
  for s in 1 2; do
    log "grandchild board: family $j sub-family $s at the family's refined poses"
    env KTAG=${L3}_fam${j}_$s INIT=pca NSUB=0 SUBSET=$E/family${s}_subset.txt NFREE=0 TABLE=$FB$j/refined_state1.txt FOCUS=${FOCUS:-1} PMARG=${PMARG:-3} DEVW=0.05 PNOISE=1 LPS="${LPS:-12 12 10 10 8 8 8 8 8}" SIGMAS="${SIGMAS:-1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7}" bash $W/run_klem.sbatch 2>&1 | grep "subset\|focus mask\|final:\|states: cross-half\|Error\|Traceback" | sed "s/^/[fam $j.$s] /"
    GC="$GC $W/lat/klem${L3}_fam${j}_$s"; $PY $W/state_coverage.py $W/lat/klem${L3}_fam${j}_$s/klem_latent.txt $W/sim/gt_particles.txt --argmax --tag "grandchild $j.$s" | tail -1
  done
  CH="$CH $CB$j"
done
$PY $W/expand_families.py merge $W/lat/expand${L3}_fam1 $GC --gt $W/sim/gt_particles.txt --minrep 0.5 | sed 's/^/[level 3] /'
echo "--- level 2 (the feedback children) for reference:"; $PY $W/expand_families.py merge $W/lat/expand2b $CH --gt $W/sim/gt_particles.txt --minrep 0.5 | grep MERGED | sed 's/^/[level 2] /'
log DONE
