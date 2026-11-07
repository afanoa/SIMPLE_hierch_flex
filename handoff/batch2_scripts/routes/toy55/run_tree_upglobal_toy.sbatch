#!/bin/bash
#SBATCH -J toy10_tree_upglobal
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# "Traversal back up": the four level-3 grandchildren of the balanced realigned tree are rebuilt exactly as before, then
# each scores ALL on-grid particles against its references (one global E-step, --final_all); the merge gives every
# particle the leaf, anywhere in the tree, where its responsibility is largest. Against the hard tree (0.71) this
# measures what re-assignment across the whole tree at the leaves is worth. env: CB, TABLE, BAL (0.5), L3 (_l3g)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; cd $W; log(){ echo "[$(date '+%F %T') upg] $*"; }
CB=${CB:-$W/lat/klem_expand_realign1_bal05_fam}; TABLE=${TABLE:-$W/lat/realign_realign1/ptcl3D_realigned.txt}; L3=${L3:-_l3g}; GC=""
for j in 1 2; do
  E=$W/lat/expand${L3}_fam$j; mkdir -p $E; $PY $W/expand_families.py split $CB$j $E --nfam 2 | sed "s/^/[fam $j split] /"
  for s in 1 2; do
    log "grandchild $j.$s rebuilt, then one E-step over all particles"
    env KTAG=${L3}_fam${j}_$s INIT=pca NSUB=0 SUBSET=$E/family${s}_subset.txt NFREE=0 TABLE=$TABLE FOCUS=1 PMARG=3 DEVW=0.05 PNOISE=1 BAL=${BAL:-0.5} LPS="12 12 10 10 8 8 8 8 8" SIGMAS="1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7" EXTRA="--final_all" bash $W/run_klem2.sbatch 2>&1 | grep "subset\|final E-step\|final:\|states: cross-half\|Error\|Traceback" | sed "s/^/[fam $j.$s] /"
    GC="$GC $W/lat/klem${L3}_fam${j}_$s"
  done
done
$PY $W/expand_families.py merge $W/lat/expand${L3}_fam1 $GC --gt $W/sim/gt_particles.txt --minrep 0.5 --bymaxr | sed 's/^/[global leaves] /'
log DONE
