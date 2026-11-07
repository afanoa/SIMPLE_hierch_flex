#!/bin/bash
#SBATCH -J subcls_node2
#SBATCH -p csbdevel
#SBATCH -c 32
#SBATCH --mem=120G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Sub-class route at one node of a recursive flex run (rflex4 layout): the node's members at the poses its flex used (the parent's
# refinement) -> NBINS direction bins over the members (dir_bins.py) -> fresh project (coupled snapshot = cls_expansion pose mean):
# state 1 = members, class = bin -> make_cavgs -> cls_expansion NSUB per bin, vol1 = the node mean -> subcls_embed.py (the node's
# flex latent pooled over the sub-classes: split-half structure readouts, children, truth). Times every stage.
# env RUN (run dir under rnasep_klem, e.g. rflex_s30a), NODE (r22), TAG (truth sim<TAG>; '' = none), STK (default routes/toy55/sim<TAG>/all.mrcs),
#     NBINS (100), NSUB (8), MSKDIAM (220), OUT (default rnasep_klem/subcls_<run>_<node>_b<NBINS>s<NSUB>), CODES (swing:200:230) or TTABLE (per-conformer truth)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
export SIMPLE_QSYS=local OMP_NUM_THREADS=32 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem
B1=$HOME/SIMPLE_focus/snap/focus3; B2=$HOME/SIMPLE_coupled/snap/9eaef62df   # B1 reads the run's projects, B2 = cls_expansion pose mean
RUN=${RUN:?}; NODE=${NODE:?}; TAG=${TAG-swing30}; CODES=${CODES-swing:200:230}; NBINS=${NBINS:-100}; NSUB=${NSUB:-8}; MD=${MSKDIAM:-220}; SMPD=1.464
RD=$K_/$RUN; ND=$RD/node_$NODE; O=${OUT:-$K_/subcls_${RUN#rflex_}_${NODE}_b${NBINS}s${NSUB}}; mkdir -p $O && cd $O
STK=${STK:-$T/sim$TAG/all.mrcs}; GT=$T/sim$TAG/gt_particles.txt
log(){ echo "[$(date '+%F %T') subcls $RUN $NODE] $*"; }
set -- $(awk -v u=$NODE '$1 == u' $RD/queue.txt); PP=${3:-}; ROWS=${4:-}
[ -s "$PP" ] && [ -s "$ROWS" ] || { log "node $NODE not in $RD/queue.txt"; exit 1; }
MEAN=$(ls -td $ND/mean/*_reconstruct3D | head -1)/recvol_state01.mrc; FLEX=$(ls -td $ND/*_flex_pca | head -1)/flex_pca_coordinates.txt
[ -s $MEAN ] && [ -s $FLEX ] && [ -s $STK ] || { log "node mean, flex coordinates or stack missing"; exit 1; }
log "node $NODE: $(awk '$1 > 0' $ROWS | wc -l) members, poses from $PP; NBINS $NBINS, NSUB $NSUB; out $O"
t0=$(date +%s)
if [ ! -s $O/all_binned.txt ]; then
  ( export SIMPLE_PATH=$B1 PATH=$B1/bin:$PATH; simple_exec prg=print_project_field projfile=$PP oritype=ptcl3D mkdir=no > $O/ptcl3D_parent.txt 2> /dev/null )
  $PY - <<PYEOF || { log "node oritab FAILED"; exit 1; }
import numpy as np
lines = [ln.split() for ln in open('$O/ptcl3D_parent.txt') if 'e1=' in ln]; rows = np.loadtxt('$ROWS')
assert len(lines) == len(rows), 'dump rows %d vs node rows %d' % (len(lines), len(rows))
with open('$O/oris_node.txt', 'w') as f:
    for kv, r in zip(lines, rows): f.write(' '.join(t for t in kv if not t.startswith(('state=', 'class='))) + ' state=%d\n' % (1 if r > 0 else 0))
B = 160; c = B // 2; zz, yy, xx = np.mgrid[:B, :B, :B]; rr = np.sqrt((xx - c) ** 2. + (yy - c) ** 2. + (zz - c) ** 2.); R = $MD / 2. / $SMPD
import sys; sys.path.insert(0, '$T'); from model_class import write_mrc
write_mrc('$O/support.mrc', np.clip((R + 3 - rr) / 6., 0, 1).astype(np.float32), $SMPD); print('node oritab %d rows, %d members; support sphere written' % (len(rows), (rows > 0).sum()))
PYEOF
  $PY $A/toyL/dir_bins.py $O/oris_node.txt $NBINS $O/all_binned.txt $O/bins.npz || { log "binning FAILED"; exit 1; }
fi
t1=$(date +%s)
export SIMPLE_PATH=$B2 PATH=$B2/bin:$PATH
if [ ! -s $O/p/p.simple ]; then
  rm -rf $O/p; simple_exec prg=new_project projname=p > $O/newproj.log 2>&1; cd $O/p
  simple_exec prg=import_particles projfile=p.simple stk=$STK oritab=$O/all_binned.txt smpd=$SMPD kv=300 cs=2.7 fraca=0.1 ctf=yes mkdir=no > import.log 2>&1 || { log "import FAILED: $(grep -a -m2 ERROR import.log | cut -c1-160)"; exit 1; }
  # import ignores the oritab's state: the node's members are set by selection (ptcl2D + ptcl3D, no pruning: rows keep their indices)
  simple_exec prg=selection projfile=p.simple oritype=ptcl2D infile=$ROWS mkdir=no > sel.log 2>&1 || { log "selection FAILED: $(grep -a -m2 'ERROR\|WARN' sel.log | cut -c1-160)"; exit 1; }
fi
cd $O/p
if [ ! -s $O/p/cls_expansion_labels.txt ]; then
  log "make_cavgs: $NBINS bin averages"; simple_exec prg=make_cavgs projfile=p.simple ncls=$NBINS mskdiam=$MD refs=cavgs_bins.mrc nparts=1 nthr=32 mkdir=no > MAKE_CAVGS_OUTPUT 2>&1 || { log "make_cavgs FAILED: $(grep -a -m2 'ERROR\|THROW' MAKE_CAVGS_OUTPUT | cut -c1-160)"; exit 1; }
  t2=$(date +%s)
  log "cls_expansion: pose mean (vol1 = node mean), $NSUB sub-classes per bin, lp 20"
  simple_exec prg=cls_expansion projfile=p.simple ncls=$NSUB lp=20 mskdiam=$MD vol1=$MEAN pcg_mskfile=$O/support.mrc nthr=32 mkdir=no > CLS_EXPANSION_OUTPUT 2>&1 || { log "cls_expansion FAILED: $(grep -a -m2 'ERROR\|THROW\|fopen' CLS_EXPANSION_OUTPUT | cut -c1-160)"; exit 1; }
  t3=$(date +%s); log "timing: bins $((t1 - t0)) s, import + make_cavgs $((t2 - t1)) s, cls_expansion $((t3 - t2)) s"
fi
log "sub-classes: $(grep -vc '^#' cls_expansion_labels.txt) labelled particles, $(awk '!/^#/{print $3}' cls_expansion_labels.txt | sort -u | wc -l) sub-classes; first bins: $(grep -m3 'starting: class' CLS_EXPANSION_OUTPUT | awk '{printf "%s ", $NF}')(bins.npz members: $($PY -c "import numpy as np; print(*np.load('$O/bins.npz')['count'][:3])"))"
REFS="rflex_refined_states=$ND/refined_states.txt,rflex_flex_states=$ND/flex_states.txt"; [ -s $ND/refined_states.txt ] || REFS=""; [ -s $RD/coord_$NODE/leaf_states.txt ] && REFS="$REFS,rflex_coord_leaves=$RD/coord_$NODE/leaf_states.txt"
t4=$(date +%s)
TRU=""; [ -n "$TAG" ] && TRU="--gt $GT $([ -n "${CODES:-}" ] && echo --codes $CODES) $([ -n "${TTABLE:-}" ] && echo --ttable $TTABLE) $([ -n "$REFS" ] && echo --refs ${REFS#,})"
$PY $T/subcls_embed2.py --labels $O/p/cls_expansion_labels.txt --flex $FLEX --rows $ROWS $TRU --out $O/embed 2>&1 | grep -v Warning
log "embedding $(( $(date +%s) - t4 )) s; DONE"
