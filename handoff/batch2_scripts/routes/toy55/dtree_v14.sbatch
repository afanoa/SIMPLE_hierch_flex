#!/bin/bash
#SBATCH -J klem_dtree
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 24:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Data-derived tree on a pair of independent half-sets (env MODE -> synth_<MODE>_A/_B, TAG, KC candidate groups per node (6),
# MINGRP smallest candidate group (400), NMIN particles per half for a child to be split again (2000), MAXDEPTH (4), NSHUF (3),
# PARNODE nodes at a time (2)). Per node, both halves: a board on the node's particles (band finer with depth), <= KC candidate
# groups, equal-N reconstructions from the common consensus poses of the real groups and of NSHUF label shuffles; the largest
# shuffled cross-tree pair statistic is the node's threshold; the real merge (archetype pruning; pairing on the coarse band) gives K*:
# K* >= 2 -> the children recurse, K* = 1 -> leaf. TREE.txt: sizes, statistics, decisions and the truth composition.
# v14: data-derived schedule: node_refine2 (matching low-pass from the half's even/odd FSC every iteration, bounds 25..6 A)
# writes the refined states' FSC 0.5 / 0.143 resolution; every child board (dtree_board6) runs from its parent's FSC 0.5
# to its FSC 0.143 resolution, per half; CANDSRC=refined tests the refinement's re-classified states instead of the board groups.
# v13: per-node pose refinement from depth REFDEPTH (1) with the f692e family recipe (node_refine.sbatch: k-means K on the
# node board -> refine3D_states, local, lpstart 20/14/12 by depth -> 8 A), per half; the node's split test and its children's
# boards use the refined poses (PROJ / POSEA,POSEB); RESUME_ROOT=<tag> reuses that tree's root decision. Use with
# BOARDJOB=dtree_board5.sbatch and BENCHJOB=splitbench_node6.sbatch.
# v12: BOARDJOB (dtree_board4.sbatch + CTFABS/PSCALE for real phase-flipped data), CONFS tokens k@n keep at most n
# particles of conformer k per half, NOGT=1 skips the truth readout of the junk flag.
# v11: junk flag (JUNKDEPTH 0 = root only, JUNKF 0.8): a child whose particles the node's board put on its free components
# (share above JUNKF in both halves) becomes a junk leaf with junk.flag and is not split further.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; K=$A/rnasep_klem
MODE=${MODE:?}; TAG=${TAG:-dtree}; KC=${KC:-6}; MINGRP=${MINGRP:-400}; NMIN=${NMIN:-2000}; MAXDEPTH=${MAXDEPTH:-4}; NSHUF=${NSHUF:-3}; PARNODE=${PARNODE:-2}
RA=$K/synth_${MODE}_A; RB=$K/synth_${MODE}_B; CONS=$W/maps160/toy_consensus.mrc; T=$K/dtree_${MODE}_$TAG; mkdir -p $T
log(){ echo "[$(date '+%F %T') dtree $MODE] $*"; }
waitjobs(){ while squeue -h -j "$(echo $* | tr ' ' ',')" 2>/dev/null | grep -q .; do sleep 15; done; }
rec(){ local D=$1 R=$2 ns=$3; mkdir -p $D/rec && cd $D/rec && cp -f $R/src.simple p.simple && \
  simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$D/states.txt mkdir=no > sel.log 2>&1 && \
  simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=$ns nparts=2 nthr=8 > rec.log 2>&1 && \
  RR=$(ls -td ./*_reconstruct3D | head -1) && mkdir -p $D/maps && for s in $(seq -f "%02g" 1 $ns); do for h in "" _even _odd; do cp -f $RR/recvol_state$s$h.mrc $D/maps/ 2>/dev/null; done; done; }
node(){ local nd=$1 depth=$2; local NA=$RA/$TAG/$nd NB=$RB/$TAG/$nd S=$T/node_$nd.txt; local ja jb thr kstar
  local PA=$RA/src.simple PB=$RB/src.simple POA=$W/lat_${MODE}_A/ptcl3D_r0.txt POB=$W/lat_${MODE}_B/ptcl3D_r0.txt; rm -f $S.ref $S.band
  [ -e $NA/proj.simple ] && PA=$NA/proj.simple; [ -e $NB/proj.simple ] && PB=$NB/proj.simple; [ -e $NA/poses_full.txt ] && POA=$NA/poses_full.txt; [ -e $NB/poses_full.txt ] && POB=$NB/poses_full.txt
  local RESUMED=0; [ -s $NA/children.txt ] && [ -s $NB/children.txt ] && [ -s $S ] && RESUMED=1
  if [ $RESUMED = 1 ]; then kstar=$(grep -a "^K\* = " $NA/bench/bench.log | tail -1 | grep -o "[0-9]*$"); kstar=${kstar:-1}; else
  local nl=$(awk -v a=$(grep -c . $NA/subset.txt) -v b=$(grep -c . $NB/subset.txt) 'BEGIN{m = (a < b ? a : b); n = int(sqrt(m / 1000) + 0.5); if (n < 3) n = 3; if (n > 7) n = 7; print n}')   # lattice side from the node size (~1000 particles per cell row: 3x3 .. 7x7)
  local bra="" brb=""; [ -e $NA/poses_full.txt ] && bra=$(cat $(dirname $(readlink -f $NA/poses_full.txt))/res.txt 2>/dev/null | awk '$1 > 0 && $2 > 0 {print $1 ":" $2}'); [ -e $NB/poses_full.txt ] && brb=$(cat $(dirname $(readlink -f $NB/poses_full.txt))/res.txt 2>/dev/null | awk '$1 > 0 && $2 > 0 {print $1 ":" $2}')
  [ -n "$bra$brb" ] && echo "   board band from the parent refinement (FSC 0.5:0.143, A): A ${bra:-depth default} B ${brb:-depth default}" >> $S.band
  local bj=""; for sd in $(seq 1 ${NSEED:-4}); do bj="$bj $(sbatch --parsable --export=ALL,R=$RA,ND=$NA,DEPTH=$depth,N1=$nl,N2=$nl,SEED=$((sd + ${SEEDOFF:-0})),PROJ=$PA,BANDRES=$bra $W/${BOARDJOB:-dtree_board6.sbatch}) $(sbatch --parsable --export=ALL,R=$RB,ND=$NB,DEPTH=$depth,N1=$nl,N2=$nl,SEED=$((sd + ${SEEDOFF:-0})),PROJ=$PB,BANDRES=$brb $W/${BOARDJOB:-dtree_board6.sbatch})"; done; waitjobs $bj
  for N in $NA $NB; do local best=""; local bv=-1   # the most confident start: highest mean max responsibility (one half's own data: no cross-half information)
    for sd in $(seq $((1 + ${SEEDOFF:-0})) $((${NSEED:-4} + ${SEEDOFF:-0}))); do f=$N/board_s$sd/klem_latent.txt; [ -s $f ] || continue; v=$(awk '!/^#/ {s+=$5; n++} END{printf "%.5f", s/n}' $f); echo "$sd $v" >> $N/board_scores.txt
      if awk -v a=$v -v b=$bv 'BEGIN{exit !(a > b)}'; then bv=$v; best=$sd; fi; done
    [ -n "$best" ] && ln -sfn board_s$best $N/board; done
  if [ ! -s $NA/board/klem_latent.txt ] || [ ! -s $NB/board/klem_latent.txt ]; then echo "== node $nd: board FAILED" > $S; return 1; fi
  if [ $depth -ge ${REFDEPTH:-1} ]; then local lps=${LPSHI:-25}   # per-node refinement; the matching low-pass follows the half's own FSC within ${LPSHI:-25}..${LPSTOP:-6} A
    local rj="$(sbatch --parsable --export=ALL,R=$RA,PROJ=$PA,POSES=$POA,ND=$NA,LPS=$lps,LPSTOP=${LPSTOP:-6} $W/${REFJOB:-node_refine2.sbatch}) $(sbatch --parsable --export=ALL,R=$RB,PROJ=$PB,POSES=$POB,ND=$NB,LPS=$lps,LPSTOP=${LPSTOP:-6} $W/${REFJOB:-node_refine2.sbatch})"; waitjobs $rj
    [ -e $NA/refined/proj.simple ] && { PA=$NA/refined/proj.simple; POA=$NA/refined/poses_full.txt; }; [ -e $NB/refined/proj.simple ] && { PB=$NB/refined/proj.simple; POB=$NB/refined/poses_full.txt; }
    for j in $rj; do grep -a -h "k-means\|\[refined\]\|FAILED\|final FSC" $K/node_refine2.$j.out $K/node_refine.$j.out 2>/dev/null | sed "s/^/   refinement: /"; done > $S.ref; fi
  for H in A B; do if [ $H = A ]; then R=$RA; N=$NA; else R=$RB; N=$NB; fi; NROWS=$($PY -c "import numpy as np; print(len(np.load('$R/meta_toy.npz')['flags']))")
    for kind in cand $( [ "${BACKEND:-bench}" = bench ] || seq -f "shuf%g" 1 $NSHUF); do D=$N/fam0/region$kind/level1; mkdir -p $D
      if [ "${CANDSRC:-board}" = refined ] && [ $kind = cand ] && [ -s $N/refined/states_refined.txt ]; then $PY -c "import numpy as np; s = np.loadtxt('$N/refined/states_refined.txt').astype(int); c = np.bincount(s, minlength=2); keep = [k for k in range(1, len(c)) if c[k] >= $MINGRP]; o = np.zeros_like(s)
for i, k in enumerate(keep): o[s == k] = i + 1
np.savetxt('$D/states.txt', o, fmt='%d'); print('refinement states as candidates: > %d groups (sizes %s; below $MINGRP dropped: %s)' % (len(keep), [int(c[k]) for k in keep], [int(c[k]) for k in range(1, len(c)) if k not in keep]))" > $D/cand.log 2>&1
      else $PY $W/${CANDS:-split_candidates2.py} $N/board $NROWS $KC $D $MINGRP > $D/cand.log 2>&1; fi
      [ $kind != cand ] && $PY -c "import numpy as np; s = np.loadtxt('$D/states.txt').astype(int); on = np.where(s > 0)[0]; s[on] = np.random.default_rng(${kind#shuf}).permutation(s[on]); np.savetxt('$D/states.txt', s, fmt='%d')"
      if [ "${BACKEND:-bench}" = bench ]; then cp -f $D/states.txt $D/states_full.txt; continue; fi
      $PY $W/neq_subsample.py $D/states.txt $D/states_full.txt min > $D/neq.log 2>&1; ns=$($PY -c "import numpy as np; print(int(np.loadtxt('$D/states.txt').max()))")
      ( rec $D $R $ns > /dev/null 2>&1 || echo "rec FAILED $D" >> $S.err ) &
    done; done; wait
  if [ "${BACKEND:-bench}" = bench ]; then local jb2; C=fam0/regioncand/level1/states_full.txt   # K* = the last decision in the bench log (after soft EM)
    jb2=$(sbatch --parsable -J klem_bench --export=ALL,MODE=$MODE,POSEA=$POA,POSEB=$POB,BENCHEXTRA="--seed ${SEEDOFF:-0}" $W/${BENCHJOB:-splitbench_node4.sbatch} $NA/$C $NB/$C $NA/bench); waitjobs $jb2
    kstar=$(grep -a "^K\* = " $NA/bench/bench.log | tail -1 | grep -o "[0-9]*$"); kstar=${kstar:-1}; cp -f $NA/bench/children_A.txt $NA/children.txt; cp -f $NA/bench/children_B.txt $NB/children.txt
    { echo "== node $nd (depth $depth): A $(grep -c . $NA/subset.txt) / B $(grep -c . $NB/subset.txt) particles; candidates A $(grep -o "> [0-9]* groups" $NA/fam0/regioncand/level1/cand.log | head -1) B $(grep -o "> [0-9]* groups" $NB/fam0/regioncand/level1/cand.log | head -1) -> K* $kstar; lattice ${nl}x$nl; boards A $(readlink $NA/board) B $(readlink $NB/board)"
      grep -a "equal N\|shuffle null\|pair statistics\|archetype test\|merges at\|children" $NA/bench/bench.log | sed 's/^/   /'
      if [ "${MODE#55}" != "$MODE" ]; then for H in A B; do if [ $H = A ]; then N=$NA; else N=$NB; fi; echo "   half $H reference clades:"; $PY $W/clade_compose.py $N/children.txt $W/sim${MODE}_$H/gt_particles.txt $W/truthtree_55_200k.npz 0.45,0.40 2>&1 | grep -v Warning | grep -v "reference clades at" | sed 's/^/  /'; done; fi; } > $S
  else
  local MV="FAMS=0 MASK=sphere NPROC=16 MINREC=100 FSCMAX=30 NUIS=${NUIS:-0} STAT=${STAT:-window} MATCHBAND=${MATCHBAND:-low} NULLCAL=0.99"
  for s in $(seq 1 $NSHUF); do env $MV LEAFGLOB=regionshuf$s $PY $W/merge_cross_tree.py $NA $NB $CONS $NA/merge_shuf$s > $NA/merge_shuf$s.log 2>&1; done
  thr=$(cat $NA/merge_shuf*.log | grep -a "pair statistics at the start" | sed 's/.*sorted): //; s/; threshold.*//' | tr ' ' '\n' | grep . | sort -g | tail -1)
  env $MV LEAFGLOB=regioncand ARCH=1 THRMIN=${thr:-0} $PY $W/merge_cross_tree.py $NA $NB $CONS $NA/merge_cand > $NA/merge_cand.log 2>&1
  kstar=$(grep -a "merges -> [0-9]* states" $NA/merge_cand.log | grep -o "> [0-9]*" | grep -o "[0-9]*"); kstar=${kstar:-1}
  for H in A B; do if [ $H = A ]; then N=$NA; else N=$NB; fi; $PY $W/group_labels.py $N/fam0/regioncand/level1/states_full.txt $NA/merge_cand_$H.npz $N/children.txt > /dev/null 2>&1; done
  { echo "== node $nd (depth $depth): A $(grep -c . $NA/subset.txt) / B $(grep -c . $NB/subset.txt) particles; candidates A $(grep -o "> [0-9]* groups" $NA/fam0/regioncand/level1/cand.log | head -1) B $(grep -o "> [0-9]* groups" $NB/fam0/regioncand/level1/cand.log | head -1); shuffled-control threshold ${thr:-?} -> K* $kstar"
    grep -a "pair statistics\|archetype test\|merges ->" $NA/merge_cand.log | sed 's/^/   /'
    for H in A B; do if [ $H = A ]; then N=$NA; else N=$NB; fi; echo "   half $H:"; $PY $W/dtree_compose.py $N/children.txt $W/sim${MODE}_$H/gt_particles.txt 2>&1 | grep -v Warning | sed 's/^/  /'; done; } > $S
  fi
  fi   # RESUMED
  [ -s $S.band ] && cat $S.band >> $S; [ -s $S.ref ] && cat $S.ref >> $S
  local JK=""; if [ "$kstar" -ge 2 ] && [ $depth -le ${JUNKDEPTH:-0} ]; then
    $PY $W/junk_children.py $NA $NB $nd ${JUNKF:-0.8} $( [ "${NOGT:-0}" = 1 ] || echo $W/sim${MODE}_A/gt_particles.txt $W/sim${MODE}_B/gt_particles.txt ) > $S.junk 2>&1; JK=$(tail -1 $S.junk | grep -v "child\|rror"); head -n -1 $S.junk >> $S; fi
  if [ "$kstar" -ge 2 ] && [ $depth -lt $MAXDEPTH ]; then for c in $(seq 1 $kstar); do
    na=$(awk -v c=$c '$1 == c' $NA/children.txt | wc -l); nb=$(awk -v c=$c '$1 == c' $NB/children.txt | wc -l)
    if echo " $JK " | grep -q " $c "; then mkdir -p $RA/$TAG/$nd.$c $RB/$TAG/$nd.$c; touch $RA/$TAG/$nd.$c/junk.flag $RB/$TAG/$nd.$c/junk.flag; echo "$nd.$c" >> $T/junk_leaves.txt
      awk -v c=$c '$1 == c {print NR}' $NA/children.txt > $RA/$TAG/$nd.$c/subset.txt; awk -v c=$c '$1 == c {print NR}' $NB/children.txt > $RB/$TAG/$nd.$c/subset.txt; continue; fi
    if [ $na -ge $NMIN ] && [ $nb -ge $NMIN ]; then mkdir -p $RA/$TAG/$nd.$c $RB/$TAG/$nd.$c; ln -sfn $PA $RA/$TAG/$nd.$c/proj.simple; ln -sfn $PB $RB/$TAG/$nd.$c/proj.simple; ln -sfn $POA $RA/$TAG/$nd.$c/poses_full.txt; ln -sfn $POB $RB/$TAG/$nd.$c/poses_full.txt
      awk -v c=$c '$1 == c {print NR}' $NA/children.txt > $RA/$TAG/$nd.$c/subset.txt; awk -v c=$c '$1 == c {print NR}' $NB/children.txt > $RB/$TAG/$nd.$c/subset.txt
      echo "$nd.$c $((depth + 1))" >> $T/queue_$((depth + 1)).txt
    else echo "   child $nd.$c: A $na / B $nb particles, below $NMIN: leaf" >> $S; fi; done; fi; }
log "start: KC $KC, MINGRP $MINGRP, NMIN $NMIN, MAXDEPTH $MAXDEPTH, NSHUF $NSHUF, PARNODE $PARNODE, junk flag to depth ${JUNKDEPTH:-0} at share ${JUNKF:-0.8}, refinement from depth ${REFDEPTH:-1}, root ${RESUME_ROOT:-new}"
rm -f $T/queue_*.txt $T/node_*.txt $T/junk_leaves.txt; mkdir -p $RA/$TAG/r $RB/$TAG/r; cp -f $RA/subset_all.txt $RA/$TAG/r/subset.txt; cp -f $RB/subset_all.txt $RB/$TAG/r/subset.txt; echo "r 0" > $T/queue_0.txt
if [ -n "${CONFS:-}" ]; then for H in A B; do   # CONFS=a:b:.. (no commas: sbatch --export splits on them): the root holds only these conformers' particles (a known-structure set cut from a larger simulation)
  $PY -c "
import numpy as np
c = np.loadtxt('$W/sim${MODE}_$H/gt_particles.txt', comments='#')[:, 1].astype(int); rng = np.random.default_rng(7); keep = np.zeros(len(c), bool)
for t in '$CONFS'.split(':'):
    k, m = (t.split('@') + [''])[:2]; r = np.where(c == int(k))[0]
    if m: r = np.sort(rng.choice(r, min(int(m), len(r)), replace=False))
    keep[r] = True
np.savetxt('$K/synth_${MODE}_$H/$TAG/r/subset.txt', np.where(keep)[0] + 1, fmt='%d'); print('root $H: %d particles of conformers $CONFS (%s)' % (keep.sum(), ' '.join('%d:%d' % (k, (c[keep] == k).sum()) for k in np.unique(c[keep]))))"; done; fi
if [ -n "${RESUME_ROOT:-}" ]; then for H in A B; do if [ $H = A ]; then R=$RA; else R=$RB; fi; S0=$R/$RESUME_ROOT/r; D0=$R/$TAG/r; mkdir -p $D0
  cp -f $S0/subset.txt $S0/children.txt $D0/; for x in $S0/board_s* $S0/bench $S0/fam0; do [ -e $x ] && ln -sfn $x $D0/$(basename $x); done; ln -sfn $(readlink $S0/board) $D0/board; done
  grep -v "free components" $K/dtree_${MODE}_$RESUME_ROOT/node_r.txt > $T/node_r.txt; log "root decision reused from tree $RESUME_ROOT"; fi
for d in $(seq 0 $MAXDEPTH); do [ -s $T/queue_$d.txt ] || break; log "depth $d: $(wc -l < $T/queue_$d.txt) nodes"
  while read nd dep; do while [ $(jobs -rp | wc -l) -ge $PARNODE ]; do sleep 10; done; node $nd $dep < /dev/null & done < $T/queue_$d.txt; wait
  for f in $(ls $T/node_*.txt 2>/dev/null | sort); do grep -q "depth $d)" $f && cat $f; done; done
log "TREE"; for f in $(ls $T/node_*.txt | sort); do cat $f; done > $T/TREE.txt; cat $T/TREE.txt; log DONE
