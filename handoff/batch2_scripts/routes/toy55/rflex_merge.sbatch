#!/bin/bash
#SBATCH -J klem_rflex_merge
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=110G
#SBATCH -t 08:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# after an rflex run (synthetic, one stack): leaves as integer ids -> their bottom-up tree (state_merge --list) -> tree-weighted leaf maps vs
# truth (tree_wiener4 --list) -> per-row pose table from each leaf's own node refinement -> merge-up along the leaf tree (bu_merge.sbatch:
# splitbench9 cross-half EM, stratified null), scorecard per level. env RT, TAG, MSKDIAM (220), CODES (swing:200:1), PAR (4)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem
R=$K_/rflex_${RT:?}; TAG=${TAG:?}; MD=${MSKDIAM:-220}; P=$R/post; mkdir -p $P; GT=$T/sim$TAG/gt_particles.txt; log(){ echo "[$(date '+%F %T') rflex-merge $RT] $*"; }
awk '{print NR, $3, $4, $5}' $R/leaves.txt > $P/leaves_int.txt; log "$(grep -c . $P/leaves_int.txt) leaves"
[ -s $P/tree/merge_tree.txt ] || $PY $T/state_merge.py $P/tree --list $P/leaves_int.txt --link complete --ncut 8,4,2 2>&1 | grep -v Warning | grep "states,\|^== " | head -4
[ -s $P/tw/summary.txt ] || $PY $T/tree_wiener4.py $P/tw --list $P/leaves_int.txt --labels $R/final_labels.txt --tree $P/tree/merge_tree.txt --gt $GT --truth $T/sim$TAG/maps --mskdiam $MD 2>&1 | grep -v Warning | grep -v "smallest leaf"
# per-row poses: the refined project of the node each leaf comes from (rows outside every leaf: the root's)
if [ ! -s $P/poses_combined.txt ]; then
  for nd in $(awk '{print $6}' $R/leaves.txt | xargs -n1 dirname | sort -u) $R/node_r; do RP=$(ls -td $nd/*_refine3D | head -1)/r.simple; [ -s $P/dump_$(basename $nd).txt ] || simple_exec prg=print_project_field projfile=$RP oritype=ptcl3D mkdir=no 2>/dev/null > $P/dump_$(basename $nd).txt; done
  $PY - <<PYEOF
import numpy as np, os
L = [l.split() for l in open('$R/leaves.txt') if l.strip()]; dumps = {}
root = open('$P/dump_node_r.txt').read().splitlines(); out = list(root); n = len(root)
for l in L:
    nd = os.path.basename(os.path.dirname(l[5])); lab = np.loadtxt(l[5], dtype=int)
    if nd not in dumps: dumps[nd] = open('$P/dump_%s.txt' % nd).read().splitlines()
    for r in np.where(lab == int(l[6]))[0]: out[r] = dumps[nd][r]
open('$P/poses_combined.txt', 'w').write('\n'.join(out) + '\n'); print('pose table: %d rows from %d node refinements' % (n, len(dumps)))
PYEOF
fi
RP0=$(ls -td $R/node_r/*_refine3D | head -1)/r.simple
TAG=$TAG LEAVES=$R/final_labels.txt TREE=$P/tree/merge_tree.txt PROJ=$RP0 POSES=$P/poses_combined.txt OUT=$P/bu GT=$GT CODES=${CODES:-swing:200:1} PAR=${PAR:-4} bash $T/bu_merge.sbatch 2>&1 | grep -v "^\[.*eo:" | tail -40
log DONE
