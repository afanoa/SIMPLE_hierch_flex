#!/bin/bash
#SBATCH -J klem_rflex_final
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# deliverables of a finished rflex tree WITHOUT a cross-branch merge-up (the node split tests did the merging): (1) per-row pose table from
# each leaf's own node refinement -> per-half poses; (2) one unit per leaf (leaves_int / final_leaves / final_junk) -> tree-weighted final maps
# (tw_units.py: tree2 or emp2 by the tree's R^2, box 160, lp 10, unsharpened, halves in halves/) in view_rflex/<RT>_merged_{tree2,emp2};
# (3) score-space embedding against all leaf maps (score_embed.py) in rnasep_klem/score_embed_<RT>. env RT
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$PATH OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 NUMBA_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; R=$K_/rflex_${RT:?}; P=$R/final16; mkdir -p $P && cd $P
H=$T/simrnasep692; MRA=${H}_A/meta_rows.txt; MRB=${H}_B/meta_rows.txt; log(){ echo "[$(date '+%F %T') final $RT] $*"; }
[ -s $R/leaves.txt ] && [ -s $R/final_labels.txt ] || { log "tree not finished (leaves.txt / final_labels.txt missing)"; exit 1; }
awk '{print NR, $3, $4, $5}' $R/leaves.txt > $P/leaves_int.txt; NL=$(grep -c . $P/leaves_int.txt)
{ echo "# unit leaves"; awk '{print NR, NR}' $R/leaves.txt; } > $P/final_leaves.txt; awk '$8 == "junk" {print NR}' $R/leaves.txt > $P/final_junk.txt
log "$NL leaves ($(grep -c . $P/final_junk.txt) junk)"
if [ ! -s $P/poses_B.txt ]; then
  : > $P/leafproj.txt; i=0   # each leaf's project: its refine dir, or (a multi-state leaf in leafmaps/) the node's stage-B refinement
  while read nm mp rest; do i=$((i + 1)); dir=$(dirname $mp)
    if [ "$(basename $dir)" = leafmaps ]; then pj=$(ls -td $(dirname $dir)/fstageB/*_refine3D | head -1); else pj=$dir; fi
    echo "$i $(ls $pj/*.simple | head -1)" >> $P/leafproj.txt; done < $R/leaves.txt
  for f in $(awk '{print $2}' $P/leafproj.txt | sort -u); do k=$(echo $f | md5sum | cut -c1-12); [ -s $P/dump_$k.txt ] || simple_exec prg=print_project_field projfile=$f oritype=ptcl3D mkdir=no > $P/dump_$k.txt 2> /dev/null; done
  $PY - <<PYEOF
import numpy as np, os, hashlib
L = [l.split() for l in open('$R/leaves.txt') if l.strip()]; out = None; PJ = dict(l.split() for l in open('$P/leafproj.txt') if l.strip())
for i, l in enumerate(L, 1):
    k = hashlib.md5((PJ[str(i)] + '\n').encode()).hexdigest()[:12]; d = open('$P/dump_%s.txt' % k).read().splitlines()
    if out is None: out = list(d)
    lab = np.loadtxt(l[5], dtype=int)
    for r in np.where(lab == int(l[6]))[0]: out[r] = d[r]
open('$P/poses_combined.txt', 'w').write('\n'.join(out) + '\n'); print('pose table: %d rows' % len(out))
PYEOF
  for h in A B; do mr=$MRA; [ $h = B ] && mr=$MRB; awk 'NR == FNR {l[FNR] = $0; next} {print l[$1 + 0]}' $P/poses_combined.txt $mr > $P/poses_$h.txt; done; fi
log "per-half poses: $(grep -c . $P/poses_A.txt) / $(grep -c . $P/poses_B.txt) rows"
O=$K_/tw_units/${RT}_final16; mkdir -p $O; V=$K_/view_rflex/${RT}_merged
$PY $T/tw_units.py $O --stage halves --leaves $P/leaves_int.txt --units $P/final_leaves.txt && \
  { [ -s $O/tree/merge_tree.txt ] || $PY $T/state_merge.py $O/tree --list $O/units_list.txt --link complete --ncut 2 > $O/tree.log 2>&1; } && \
  $PY $T/tw_units.py $O --stage maps --junk $P/final_junk.txt --mapdir $V --mskdiam 232 2>&1 | grep -v Warning | tail -4 || log "tree-weighted maps FAILED"
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
[ -s $O/map_dirs.txt ] && while read k d; do $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/arm_$k.txt $(ls $d/*.mrc) 2>&1 | grep -v Warning | tail -1 | sed "s/^/   $k: /"; done < $O/map_dirs.txt
log "score-space embedding"
$PY -u $T/score_embed.py --sim ${H}_A,${H}_B --cache $T/benchcache_rnasep692 --poses $P/poses_A.txt,$P/poses_B.txt --leaves $R/leaves_list.txt \
    --labels $R/final_labels.txt --out $K_/score_embed_$RT --em_lim 15 2>&1 | grep -v "Warning\|warn("
log DONE
