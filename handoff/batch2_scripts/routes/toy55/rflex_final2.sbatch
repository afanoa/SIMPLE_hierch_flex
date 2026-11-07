#!/bin/bash
#SBATCH -J klem_rflex_final2
#SBATCH -p csbdevel
#SBATCH -c 48
#SBATCH --mem=200G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# rflex_final with PLAIN leaf maps: every node's leaves reconstructed together (reconstruct3D, one state per leaf, no pcg focus) at the
# node's stage-B poses, so the delivered maps keep the mobile tip that the focused stage-B maps damp beyond the falloff mask. Then
# (1) tree-weighted final maps (tw_units.py, unit = leaf) in view_rflex/<RT>_final_{tree2,emp2} (raw / own in tw_units/<RT>_plain),
# arm readouts and the radial extent of the most closed maps; (2) score-space embedding against the plain leaf maps (score_embed_<RT>)
# with its figure (viz_<RT>/<RT>_scorespace.png). Node steps run 2 at a time (3 parts x 8 threads each). env RT
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=8
export PY=/mnt/nasapps/production/python/3.12/bin/python3 A=/mnt/beegfs/elmlund/afan RT=${RT:?}; export T=$A/routes/toy55 K_=$A/rnasep_klem
export R=$K_/rflex_$RT P=$K_/rflex_$RT/final_plain; mkdir -p $P && cd $P
H=$T/simrnasep692; export MRA=${H}_A/meta_rows.txt MRB=${H}_B/meta_rows.txt; log(){ echo "[$(date '+%F %T') final2 $RT] $*"; }; export -f log
[ -s $R/leaves.txt ] && [ -s $R/final_labels.txt ] || { log "tree not finished (leaves.txt / final_labels.txt missing)"; exit 1; }
NL=$(grep -c . $R/leaves.txt); awk '$8 == "junk" {print NR}' $R/leaves.txt > $P/final_junk.txt; { echo "# unit leaves"; awk '{print NR, NR}' $R/leaves.txt; } > $P/final_leaves.txt
# per node: its leaves as states 1..n (labels.txt over the project rows) and its stage-B refinement (a single-state leaf's map dir)
$PY - <<PYEOF || { log "node label tables FAILED"; exit 1; }
import numpy as np, re, os, glob
L = [l.split() for l in open('$R/leaves.txt') if l.strip()]; nodes = {}; LB = {}
for i, l in enumerate(L, 1): nodes.setdefault(re.search(r'/node_([^/]+)/', l[1]).group(1), []).append(i)
with open('$P/plain_nodes.txt', 'w') as f:
    for u, ids in nodes.items():
        lab = None; W = '$P/plain_%s' % u; os.makedirs(W, exist_ok=True)
        for j, i in enumerate(ids, 1):
            l = L[i - 1]; LB[l[5]] = LB[l[5]] if l[5] in LB else np.loadtxt(l[5], dtype=int); m = LB[l[5]] == int(l[6])
            lab = np.zeros(len(m), int) if lab is None else lab; assert not lab[m].any(), 'leaf %s overlaps another leaf of node %s' % (l[0], u); lab[m] = j
        sd = [os.path.dirname(L[i - 1][1]) for i in ids if '/leafmaps/' not in L[i - 1][1]]
        rd = sd[0] if sd else max(glob.glob('$R/node_%s/fstageB/*_refine3D' % u), key=os.path.getmtime)
        np.savetxt(W + '/labels.txt', lab, fmt='%d'); open(W + '/stageB_dir.txt', 'w').write(rd + '\n')
        f.write('%s %d %s\n' % (u, len(ids), ','.join(map(str, ids)))); print('node %s: %d leaves, %d particles, poses %s' % (u, len(ids), (lab > 0).sum(), rd))
PYEOF
rec_node(){ local u=$1 n=$2 W=$P/plain_$1 pj
  ls $W/*_reconstruct3D/recvol_state$(printf %02d $n)_odd_unfil.mrc > /dev/null 2>&1 && { log "node $u: plain leaf maps exist"; return 0; }
  pj=$(ls $(cat $W/stageB_dir.txt)/*.simple | head -1); cd $W && cp -f $pj r.simple || return 1
  simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=labels.txt mkdir=no > sel.log 2>&1 && \
    simple_exec prg=update_project projfile=r.simple qsys_name=local nparts=3 nthr=8 mkdir=no > up.log 2>&1 && \
    simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=232 nstates=$n nparts=3 nthr=8 > rec.log 2>&1 \
    || { log "node $u: plain reconstruct3D FAILED: $(grep -a -h -m2 ERROR sel.log up.log rec.log 2>/dev/null | head -2)"; return 1; }
  log "node $u: $n leaf maps reconstructed (plain, poses $pj)"; }
# per-half pose table for the embedding: each row's pose from its leaf's node refinement
poses(){ [ -s $P/poses_B.txt ] && return 0
  : > $P/leafproj.txt; while read u n ids; do pj=$(ls $(cat $P/plain_$u/stageB_dir.txt)/*.simple | head -1); for i in ${ids//,/ }; do echo "$i $pj" >> $P/leafproj.txt; done; done < $P/plain_nodes.txt
  for f in $(awk '{print $2}' $P/leafproj.txt | sort -u); do k=$(echo $f | md5sum | cut -c1-12); [ -s $P/dump_$k.txt ] || simple_exec prg=print_project_field projfile=$f oritype=ptcl3D mkdir=no > $P/dump_$k.txt 2> /dev/null; done
  $PY - <<PYEOF
import numpy as np, hashlib
L = [l.split() for l in open('$R/leaves.txt') if l.strip()]; PJ = dict(l.split() for l in open('$P/leafproj.txt') if l.strip()); out = None; D = {}; LB = {}
for i, l in enumerate(L, 1):
    k = hashlib.md5((PJ[str(i)] + '\n').encode()).hexdigest()[:12]
    if k not in D: D[k] = open('$P/dump_%s.txt' % k).read().splitlines()
    if out is None: out = list(D[k])
    if l[5] not in LB: LB[l[5]] = np.loadtxt(l[5], dtype=int)
    for r in np.where(LB[l[5]] == int(l[6]))[0]: out[r] = D[k][r]
open('$P/poses_combined.txt', 'w').write('\n'.join(out) + '\n'); print('pose table: %d rows' % len(out))
PYEOF
  for h in A B; do mr=$MRA; [ $h = B ] && mr=$MRB; awk 'NR == FNR {l[FNR] = $0; next} {print l[$1 + 0]}' $P/poses_combined.txt $mr > $P/poses_$h.txt; done; }
export -f rec_node poses
log "$NL leaves in $(grep -c . $P/plain_nodes.txt) nodes ($(grep -c . $P/final_junk.txt) junk): plain reconstructions, 2 nodes at a time; pose table in parallel"
poses > $P/poses.log 2>&1 < /dev/null & PJ=$!
xargs -P 2 -L 1 bash -c 'rec_node "$0" "$1"' < $P/plain_nodes.txt; wait $PJ; tail -1 $P/poses.log
log "per-half poses: $(grep -c . $P/poses_A.txt) / $(grep -c . $P/poses_B.txt) rows"
$PY - <<PYEOF || { log "plain half-maps missing"; exit 1; }
import os, glob
L = [l.split() for l in open('$R/leaves.txt') if l.strip()]; loc = {}
for u, n, ids in (l.split() for l in open('$P/plain_nodes.txt') if l.strip()):
    rd = max(glob.glob('$P/plain_%s/*_reconstruct3D' % u), key=os.path.getmtime)
    for j, i in enumerate(map(int, ids.split(',')), 1): loc[i] = tuple('%s/recvol_state%02d_%s_unfil.mrc' % (rd, j, h) for h in ('even', 'odd'))
miss = [L[i - 1][0] for i in range(1, len(L) + 1) if i not in loc or not all(os.path.exists(p) for p in loc[i])]
if miss: print('MISSING plain half-maps: ' + ' '.join(miss)); raise SystemExit(1)
with open('$P/leaves_int.txt', 'w') as f, open('$P/leaves_list.txt', 'w') as g:
    for i, l in enumerate(L, 1): f.write('%d %s %s %s\n' % (i, loc[i][0], loc[i][1], l[4])); g.write('%s %s %s %s\n' % (l[0], loc[i][0], loc[i][1], l[4]))
print('plain half-maps for all %d leaves' % len(L))
PYEOF
export OMP_NUM_THREADS=24 OPENBLAS_NUM_THREADS=24 NUMBA_NUM_THREADS=24
log "score-space embedding against the plain leaf maps (in parallel with the tree-weighted maps)"
$PY -u $T/score_embed.py --sim ${H}_A,${H}_B --cache $T/benchcache_rnasep692 --poses $P/poses_A.txt,$P/poses_B.txt --leaves $P/leaves_list.txt \
    --labels $R/final_labels.txt --out $K_/score_embed_$RT --em_lim 15 > $P/score_embed.log 2>&1 & SJ=$!
O=$K_/tw_units/${RT}_plain; mkdir -p $O; V=$K_/view_rflex/${RT}_final
$PY $T/tw_units.py $O --stage halves --leaves $P/leaves_int.txt --units $P/final_leaves.txt && \
  { [ -s $O/tree/merge_tree.txt ] || $PY $T/state_merge.py $O/tree --list $O/units_list.txt --link complete --ncut 2 > $O/tree.log 2>&1; } && \
  $PY $T/tw_units.py $O --stage maps --junk $P/final_junk.txt --mapdir $V --mskdiam 232 2>&1 | grep -v Warning | tail -4 || log "tree-weighted maps FAILED"
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1); CONS=$T/maps160/toy_consensus.mrc
[ -s $O/map_dirs.txt ] && while read k d; do $PY $T/arm_readout.py $CONS $CL $O/arm_$k.txt $(ls $d/*.mrc) 2>&1 | grep -v Warning | tail -1 | sed "s/^/   $k: /"; done < $O/map_dirs.txt
if [ -s $O/arm_raw.txt ]; then log "most closed plain leaf maps (raw):"
  grep -v '^#' $O/arm_raw.txt | sort -k2 -gr | head -10 | while read p c o call; do i=$(basename $p .mrc); i=$((10#${i##*_})); echo "   $(awk -v i=$i 'NR == i {printf "%-7s %6d ptcls", $1, $5}' $R/leaves.txt): closed-arm $c open-arm $o"; done
  $PY $T/radial_extent.py $CONS $CL --maps $(grep -v '^#' $O/arm_raw.txt | sort -k2 -gr | head -3 | awk '{n = $1; sub(/.*\//, "", n); sub(/\.mrc/, "", n); printf "%s:%s ", n, $1}') 2>&1 | grep -v Warning; fi
wait $SJ; grep -v "Warning\|warn(" $P/score_embed.log | tail -6
mkdir -p $K_/viz_$RT; [ -s $K_/score_embed_$RT/embed.npz ] && [ -s $O/arm_raw.txt ] && $PY $T/plot_score_st.py $K_/score_embed_$RT $O/arm_raw.txt $R/leaves.txt $K_/viz_$RT/${RT}_scorespace.png 2>&1 | grep -v Warning
log DONE
