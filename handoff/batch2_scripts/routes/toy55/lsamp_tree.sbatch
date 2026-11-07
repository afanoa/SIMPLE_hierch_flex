#!/bin/bash
#SBATCH -J klem_lsamp_tree
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=160G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# landscape sampling of every terminal node of a finished tree, one node after the other: lsamp2.py on the node's leaves (nodes = leaf-name
# prefixes before "_s" with >= 2 leaves, or NODES), then (REC 1) plain reconstruct3D of the node's samples at the node's refined poses with
# lsamp_rec.sbatch (lp-10 unsharpened maps in VIEW/<node>, halves in halves/; real data: arm readout, old-closed tracer, radial extent).
# env TREE (dir with leaves_list.txt + final_labels.txt), SC (lsamp_score2 dir), OUT, SIM, MR, CP, POSES, CTFABS, MSKDIAM, ARGS (extra lsamp2.py
# flags, e.g. the GT/codes; commas allowed), NODES (optional list), REC (1), VIEWBASE (rnasep_klem/view_rflex/lsamp_<name>)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; R=${TREE:?}; O=${OUT:?}; mkdir -p $O
VB=${VIEWBASE:-$K_/view_rflex/lsamp_$(basename $O)}; log(){ echo "[$(date '+%F %T') lsamp_tree $(basename $O)] $*"; }
NODES=${NODES:-$(awk '{n = $1; sub(/_[sg][0-9]+$/, "", n); c[n]++} END {for (k in c) if (c[k] >= 2) print k}' $R/leaves_list.txt | sort)}
log "terminal nodes: $NODES"
for u in $NODES; do
  D=$O/node_$u; mkdir -p $D
  if [ ! -s $D/samples.txt ]; then log "node $u: lsamp2.py"
    $PY -u $T/lsamp2.py --scores ${SC:?} --sim ${SIM:?} --mr ${MR:?} --cache ${CP:?} --poses ${POSES:?} --ctfabs ${CTFABS:-0} --mskdiam ${MSKDIAM:-230} --out $D --use ${u}_ \
        --naxes 0 --tag $u ${ARGS:-} 2>&1 | grep --line-buffered -v "Warning\|warn(" > $D/lsamp.log
    grep "sampled axes\|resolvable\|CHOSEN\|steps CV\|check 3\|TRUTH .*span\|TRUTH closed\|no landscape" $D/lsamp.log | sed "s/^/  $u: /" | cut -c1-400; fi
  [ -s $D/samples.txt ] || { log "node $u: no samples"; continue; }
  if [ "${REC:-1}" = 1 ]; then
    PJ=$(dirname $(awk -v u=$u '{n = $1; sub(/_[sg][0-9]+$/, "", n); if (n == u) {print $2; exit}}' $R/leaves_list.txt)); PROJ=$(ls $PJ/*.simple 2>/dev/null | head -1)
    [ -n "$PROJ" ] || PROJ=$(ls -td $(dirname $PJ)/fstageB/*_refine3D/*.simple $(dirname $PJ)/*_refine3D/*.simple 2>/dev/null | head -1)
    PROJ=$PROJ LAB=$D/samples.txt OUT=$D/rec VIEW=$VB/$u NAME=lsamp_$u MSKDIAM=${MSKDIAM:-230} GT=${GT:-} bash $T/lsamp_rec.sbatch 2>&1 | grep -v Warning | sed "s/^/  $u: /" | cut -c1-300
  fi
done
log "DONE"
