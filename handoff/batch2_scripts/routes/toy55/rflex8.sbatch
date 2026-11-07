#!/bin/bash
#SBATCH -J klem_rflex8
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=240G
#SBATCH -t 24:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# recursive flex with data-derived focus masks (the integrated hierarchy, first version; resumable). Per node u (its particles at
# the poses of its parent's refinement): (1) node mean (reconstruct3D); (2) flex_pca on the node: K states from the CTF/pose-aware
# latent (rec_states: half-maps); (3) focus mask from those states (focus_mask.py: cross-half variance after per-shell regression
# on the node mean, significant voxels holding 80 % of it); (4) focused K-state refine3D from the flex states (prob_neigh geom, hard
# assignments; pcg_focus: node mean fixed outside the mask, per-state deviation inside; nspace 5000, box_crop 88); (5) tree of the
# refined maps (state_merge.py, complete linkage): top branches with >= MINP particles recurse while depth < MAXD, the others'
# states are leaves. Final: per-row leaf labels, scorecard + per-leaf truth (GT), the tree over the leaves.
# rflex8: the SPLIT TEST decides every node (splitbench13 on the REAL half-stacks of the node's refined states: junk screen, cross-half EM
# <= 20 A, re-test beyond the EM band with docked group maps, stratified null): children = the supported groups (families >= MINP recurse while
# depth < MAXD; one supported group = the node is one leaf), junk units are junk leaves (decide_real.py). No per-node tree descent.
# rflex7: COARSE-TO-FINE refinement per node: stage A ITA (3) iterations at BOXA 64 / NSPA 2500 / LPA 15 A from the flex states, full-box
# references reconstructed from stage A's assignments and poses, stage B ITB (3) iterations at the regular 88 / 5000 / LP 10 A; flex on all
# the node's particles (NFLEX default 10^9). rflex6: OVER-SPLIT for a later merge-up along the tree: recursion floor MINP 15000 (~ the closed-state size), depth MAXD 4, states
# per node K_node = min(K, max(2, particles / KFLOOR 15000)) so no state starts below ~15k. rflex5: FAST: refinement at BOXC (64) / NSPACE (2500), flex fitted on at most NFLEX (120000) of the node's particles (the rest start
# from random labels; the refinement's first E-step against the flex maps assigns them), sibling nodes run concurrently (PAR 2), every
# heavy step distributed as slurm worker jobs (NP 8 parts x 8 threads); per-node outputs children.txt / node_leaves.txt gathered level by level. rflex4: truncated per node (refine MAXITS 4, flex NPROBE 2 probe iterations); a node whose refinement already wrote its maps is
# taken as refined (resume after truncation). rflex3: children by DESCENT of the node's tree: branches below MINP particles are peeled off as leaves, the descent stops where
# the tree splits into two branches >= MINP (each a family that recurses; a single state >= MINP is a family too).
# rflex2: SRC (any project; default proj<TAG>) and TAG optional (no TAG = real data: no truth readouts; the leaves' closed-arm
# occupancy instead, arm_readout.py). env RT (rnasep_klem/rflex_<RT>), TAG (proj<TAG>, sim<TAG>), K (8), MINP (50000), MAXD (2), MSKDIAM (220), LPF (flex lp, 15),
#     LP (refine lp, 10), MAXITS (8), CODES (swing:200:1), SNAP (focus3), ROOTROWS (per-row file: rows > 0 form the root; default all rows)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=${SNAP:-$HOME/SIMPLE_focus/snap/focus3}; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=${QPART:-csbdevel} OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem
TAG=${TAG:-}; RT=${RT:?}; R=$K_/rflex_$RT; K=${K:-8}; MINP=${MINP:-15000}; MAXD=${MAXD:-4}; MD=${MSKDIAM:-220}; GT=$T/sim$TAG/gt_particles.txt
SRC=${SRC:-$T/proj$TAG/proj$TAG.simple}; H=$T/simrnasep692; MRA=${H}_A/meta_rows.txt; MRB=${H}_B/meta_rows.txt; CP=$T/benchcache_rnasep692; mkdir -p $R && cd $R; log(){ echo "[$(date '+%F %T') rflex $RT] $*"; }
st(){ simple_exec prg=print_project_field projfile=$1 oritype=ptcl3D mkdir=no 2>/dev/null | awk '{s = 0; for (i = 1; i <= NF; i++) if ($i ~ /^state=/) {split($i, a, "="); s = a[2] + 0}; print s}'; }
score(){ [ -n "$TAG" ] || return 0; awk -v tag="$1" 'NR == FNR {if ($0 !~ /^#/) c[$1] = $2; next} {r++; s = $1; if (s <= 0) next; k = c[r]; n[s]++; if (k >= 200 && k <= 230) {x = k - 200; sw[s] += x; sq[s] += x * x; ns[s]++; T1 += x; T2 += x * x; NT++}}
  END {m = 0; for (s in n) if (ns[s] > 0) {mu[s] = sw[s] / ns[s]; ord[++m] = s; W += sq[s] - sw[s] ^ 2 / ns[s]}
    for (i = 1; i <= m; i++) for (j = i + 1; j <= m; j++) if (mu[ord[j]] < mu[ord[i]]) {t = ord[i]; ord[i] = ord[j]; ord[j] = t}
    printf "%s: %d swing particles, eta2 %.3f; per leaf mean +- sd swing (particles, swing share), by mean:", tag, NT, 1 - W / (T2 - T1 ^ 2 / NT); lo = 99; hi = -1
    for (i = 1; i <= m; i++) {s = ord[i]; sd = sqrt(sq[s] / ns[s] - mu[s] ^ 2); printf " %d:%.1f+-%.1f(%d,%.0f%%)", s, mu[s], sd, n[s], 100 * ns[s] / n[s]
      if (n[s] >= 200 && ns[s] / n[s] >= 0.5) {if (mu[s] < lo) lo = mu[s]; if (mu[s] > hi) hi = mu[s]}}
    printf "\n   %s extent (leaves >= 200 particles, >= 50%% swing): %.1f .. %.1f A of 0 .. 30\n", tag, lo, hi}' $GT $2; }
NR=$(simple_exec prg=print_project_field projfile=$SRC oritype=ptcl3D mkdir=no 2>/dev/null | grep -c .); log "project $SRC: $NR rows; K $K, MINP $MINP, MAXD $MAXD"
[ -s $R/queue.txt ] || { if [ -n "${ROOTROWS:-}" ]; then awk '{print ($1 > 0) ? 1 : 0}' $ROOTROWS > $R/rows_r.txt; else awk -v n=$NR 'BEGIN {for (i = 1; i <= n; i++) print 1}' > $R/rows_r.txt; fi; echo "r 0 $SRC $R/rows_r.txt" > $R/queue.txt; : > $R/leaves.txt; }
node(){   # one node: project at the parent's poses -> mean -> flex (subsample) -> focus mask -> focused refinement -> tree descent
  local u=$1 dep=$2 PP=$3 ROWS=$4; local D=$R/node_$u; mkdir -p $D && cd $D
  [ -e $D/node.done ] && return 0
  local NN=$(awk '$1 > 0' $ROWS | wc -l); local K=$(awk -v n=$(awk '$1 > 0' $ROWS | wc -l) -v k=$K -v f=${KFLOOR:-15000} 'BEGIN {m = int(n / f); if (m > k) m = k; if (m < 2) m = 2; print m}')
  if [ ! -s $D/n.simple ]; then cp -f $PP n.tmp.simple && simple_exec prg=selection projfile=n.tmp.simple oritype=ptcl3D infile=$ROWS mkdir=no > sel.log 2>&1 \
    && simple_exec prg=update_project projfile=n.tmp.simple qsys_name=slurm nparts=${NP:-8} nthr=8 mkdir=no > up.log 2>&1 && mv n.tmp.simple n.simple || { log "node $u project FAILED"; return 1; }; fi
  if [ -z "$(ls $D/mean/*_reconstruct3D/recvol_state01.mrc 2>/dev/null)" ]; then mkdir -p $D/mean && cd $D/mean && cp -f $D/n.simple m.simple
    simple_exec prg=reconstruct3D projfile=m.simple pgrp=c1 mskdiam=$MD nstates=1 nparts=${NP:-8} nthr=8 > rec.log 2>&1 || { log "node $u mean FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-150)"; return 1; }; cd $D; fi
  local MEAN=$(ls -td $D/mean/*_reconstruct3D | head -1)/recvol_state01.mrc
  # flex on at most NFLEX of the node's particles
  [ -s $D/flex_rows.txt ] || awk -v nf=${NFLEX:-1000000000} -v nn=$NN 'BEGIN {srand(17); p = (nn > nf) ? nf / nn : 1} {print ($1 > 0 && rand() < p) ? 1 : 0}' $ROWS > $D/flex_rows.txt
  if [ -z "$(ls $D/*_flex_pca/flex_pca_coordinates.txt 2>/dev/null)" ]; then cp -f n.simple f.simple
    simple_exec prg=selection projfile=f.simple oritype=ptcl3D infile=$D/flex_rows.txt mkdir=no > fsel.log 2>&1 || { log "node $u flex selection FAILED"; return 1; }
    log "node $u (depth $dep, $NN particles): flex_pca on $(awk '$1 > 0' $D/flex_rows.txt | wc -l), $K states, lp ${LPF:-15}"
    simple_exec prg=flex_pca projfile=f.simple vol1=$MEAN mskdiam=$MD lp=${LPF:-15} npreimages=$K rec_states=yes min_state_frac=0 umap=no n_probe_iters=${NPROBE:-2} nparts=${NP:-8} nthr=8 > flex.log 2>&1 || { log "node $u flex_pca FAILED: $(grep -a -m2 'ERROR\|THROW' flex.log | cut -c1-200)"; return 1; }; fi
  local FD=$(ls -td $D/*_flex_pca | head -1); local FP=$FD/f.simple; [ -s $FP ] || FP=$D/f.simple
  st $FP > $D/flex_states_raw.txt
  # relabel 1..KS in label order; node particles outside the flex subsample start from a random state (the first E-step assigns them)
  awk -v k=$K 'NR == FNR {if ($1 > 0) u[$1] = 1; next} FNR == 1 {m = 0; for (s in u) ord[++m] = s + 0; for (a = 1; a <= m; a++) for (b = a + 1; b <= m; b++) if (ord[b] < ord[a]) {t = ord[a]; ord[a] = ord[b]; ord[b] = t}; for (a = 1; a <= m; a++) map[ord[a]] = a; srand(23)}
    {print ($1 > 0) ? map[$1] : 0}' $D/flex_states_raw.txt $D/flex_states_raw.txt > $D/flex_states_sub.txt
  local KS=$(awk '$1 > 0' $D/flex_states_sub.txt | sort -n | uniq | wc -l)
  paste $D/flex_states_sub.txt $ROWS | awk -v k=$KS 'BEGIN {srand(23)} {print ($1 > 0) ? $1 : (($2 > 0) ? 1 + int(k * rand()) : 0)}' > $D/flex_states.txt
  mkdir -p $D/refs; : > $D/members.txt; local VOLS="" j=0 s b fm HE HO
  for s in $(awk '$1 > 0' $D/flex_states_raw.txt | sort -n | uniq); do j=$((j + 1)); b=$D/refs/ref_state$(printf %02d $j); fm=$FD/recvol_state$(printf %02d $s)
    cp -f $fm.mrc $b.mrc && cp -f ${fm}_even.mrc ${b}_even.mrc && cp -f ${fm}_odd.mrc ${b}_odd.mrc || { log "node $u: flex state map $fm incomplete"; return 1; }; HE=${fm}_even_unfil.mrc; HO=${fm}_odd_unfil.mrc
    echo "$j $HE $HO $(awk -v s=$j '$1 == s' $D/flex_states_sub.txt | wc -l)" >> $D/members.txt; VOLS="$VOLS vol$j=$b.mrc"; done
  log "node $u: flex states $(awk '{printf "%d ", $4}' $D/members.txt)"
  [ -s $D/mask.mrc ] || $PY $T/focus_mask.py $D/mask.mrc --list $D/members.txt --mskdiam $MD 2>&1 | grep -v Warning | tail -1 | sed "s/^/node $u /"
  if [ -z "$(ls $D/stageB/*_refine3D/recvol_state01_even_unfil.mrc 2>/dev/null)" ]; then
    if [ -z "$(ls $D/stageA/*_refine3D/recvol_state01_even_unfil.mrc 2>/dev/null)" ]; then mkdir -p $D/stageA && cd $D/stageA && cp -f $D/n.simple r.simple
      simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$D/flex_states.txt mkdir=no > rsel.log 2>&1 || { log "node $u state selection FAILED"; return 1; }
      log "node $u: stage A refine3D, $KS states, lp ${LPA:-15}, maxits ${ITA:-3}, box ${BOXA:-64}, nspace ${NSPA:-2500}"
      simple_exec prg=refine3D projfile=r.simple $VOLS nstates=$KS refine=prob_neigh prob_neigh_mode=geom pgrp=c1 mskdiam=$MD lp=${LPA:-15} maxits=${ITA:-3} center=no rec_backend=pcg pcg_mskfile=$D/mask.mrc pcg_focus=yes nspace=${NSPA:-2500} box_crop=${BOXA:-64} smpd_crop=$(awk -v b=${BOXA:-64} 'BEGIN {printf "%.8f", 160 * 1.464 / b}') nparts=${NP:-8} nthr=8 > ref.log 2>&1 || { log "node $u stage A FAILED: $(grep -a -m2 ERROR ref.log | cut -c1-150)"; return 1; }; cd $D; fi
    local RDA=$(ls -td $D/stageA/*_refine3D | head -1); mkdir -p $D/stageB && cd $D/stageB && cp -f $RDA/r.simple b.simple
    if [ -z "$(ls $D/stageB/*_reconstruct3D/recvol_state01.mrc 2>/dev/null)" ]; then
      simple_exec prg=reconstruct3D projfile=b.simple pgrp=c1 mskdiam=$MD nstates=$KS nparts=${NP:-8} nthr=8 > rec.log 2>&1 || { log "node $u stage B references FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-150)"; return 1; }; fi
    local RRB=$(ls -td $D/stageB/*_reconstruct3D | head -1) VB="" q
    for q in $(seq 1 $KS); do VB="$VB vol$q=$RRB/recvol_state$(printf %02d $q).mrc"; done
    log "node $u: stage B refine3D, $KS states, lp ${LP:-10}, maxits ${ITB:-3}, box 88, nspace 5000"
    simple_exec prg=refine3D projfile=b.simple $VB nstates=$KS refine=prob_neigh prob_neigh_mode=geom pgrp=c1 mskdiam=$MD lp=${LP:-10} maxits=${ITB:-3} center=no rec_backend=pcg pcg_mskfile=$D/mask.mrc pcg_focus=yes nspace=5000 box_crop=88 smpd_crop=2.66181803 nparts=${NP:-8} nthr=8 > ref.log 2>&1 || { log "node $u stage B FAILED: $(grep -a -m2 ERROR ref.log | cut -c1-150)"; return 1; }; cd $D
  fi
  local RD=$(ls -td $D/stageB/*_refine3D | head -1); local RP=$RD/b.simple; [ -s $RP ] || RP=$D/stageB/b.simple
  [ -s $D/refined_states.txt ] || st $RP > $D/refined_states.txt
  score "node $u flex states (subsample)" $D/flex_states_sub.txt | tail -2; score "node $u refined states" $D/refined_states.txt | tail -2
  # split test on the real half-stacks: candidates = the refined states in each half's rows, poses from this node's refinement
  local TD=$D/test; mkdir -p $TD
  if [ ! -s $TD/summary.txt ]; then
    [ -s $D/poses_B.txt ] || { simple_exec prg=print_project_field projfile=$RP oritype=ptcl3D mkdir=no 2>/dev/null > $D/poses_master.txt
      awk 'NR == FNR {l[FNR] = $0; next} {print l[$1 + 0]}' $D/poses_master.txt $MRA > $D/poses_A.txt; awk 'NR == FNR {l[FNR] = $0; next} {print l[$1 + 0]}' $D/poses_master.txt $MRB > $D/poses_B.txt; rm -f $D/poses_master.txt; }
    awk 'NR == FNR {s[FNR] = $1; next} {print s[$1 + 0] + 0}' $D/refined_states.txt $MRA > $TD/cand_A.txt; awk 'NR == FNR {s[FNR] = $1; next} {print s[$1 + 0] + 0}' $D/refined_states.txt $MRB > $TD/cand_B.txt
    log "node $u: split test of $KS refined states ($(awk '$1 > 0' $TD/cand_A.txt | wc -l) / $(awk '$1 > 0' $TD/cand_B.txt | wc -l) particles in the half-stacks)"
    $PY -u $T/splitbench13.py --sim ${H}_A,${H}_B --poses $D/poses_A.txt,$D/poses_B.txt --cand $TD/cand_A.txt,$TD/cand_B.txt --out $TD --cache $CP --ctfabs 1 --pair identity --softem 0 --arch 0 --strata 1 \
      --xem 3 --retest_beyond_em 1 --em_lim ${EMLIM:-20} --em_min ${EMMIN:-1000} --dock 2 --junk_proj 0.75 --junk_ssnr 0.35 --summary $TD/summary.txt > $TD/bench.log 2>&1 || { log "node $u split test FAILED: $(grep -a -m1 -A2 'Error\|Traceback' $TD/bench.log | tail -1)"; return 1; }
  fi
  log "node $u: split test $(head -1 $TD/summary.txt | cut -c1-160)"
  $PY $T/decide_real.py $D/refined_states.txt $TD $dep $MRA $MRB --minp $MINP --maxd $MAXD > $D/groups.txt || { log "node $u: decision FAILED"; return 1; }
  : > $D/children.txt; : > $D/node_leaves.txt; local g=0 fg=0 kind np_ mem jm nm
  while read kind np_ mem jm; do
    g=$((g + 1))
    if [ "$kind" = family ]; then fg=$((fg + 1))
      awk -v m="$mem" 'BEGIN {n = split(m, a, "+"); for (i = 1; i <= n; i++) in_[a[i]] = 1} {print ($1 in in_) ? 1 : 0}' $D/refined_states.txt > $R/rows_${u}${fg}.txt
      echo "${u}${fg} $((dep + 1)) $RP $R/rows_${u}${fg}.txt" >> $D/children.txt; log "node $u: supported group $g ($np_ particles, states $mem) -> node ${u}${fg}"
    else
      local J=""; [ "$kind" = junk ] && J=" junk"
      if [ "$mem" = "${mem#*+}" ]; then b=$RD/recvol_state$(printf %02d $mem); echo "${u}_s$mem $b.mrc ${b}_even_unfil.mrc ${b}_odd_unfil.mrc $np_ $D/refined_states.txt $mem$J" >> $D/node_leaves.txt
      else nm=${u}_g$g; mkdir -p $D/leafmaps   # a supported group of several states: one leaf, its maps the particle-weighted average
        awk -v m="$mem" 'BEGIN {n = split(m, a, "+"); for (i = 1; i <= n; i++) in_[a[i]] = 1} {print ($1 in in_) ? 1 : 0}' $D/refined_states.txt > $D/leafmaps/$nm.txt
        $PY - <<PYEOF || { log "node $u: leaf map $nm FAILED"; return 1; }
import numpy as np
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64), np.fromfile(p, np.uint8, 1024 + int(h[23]))
lab = np.loadtxt('$D/refined_states.txt', dtype=int); mem = [int(x) for x in '$mem'.split('+')]; w = np.array([(lab == s).sum() for s in mem], float); w /= w.sum()
for sfx in ('', '_even_unfil', '_odd_unfil'):
    acc = None
    for s, ws in zip(mem, w):
        x, hdr = rd('$RD/recvol_state%02d%s.mrc' % (s, sfx)); acc = ws * x if acc is None else acc + ws * x
    with open('$D/leafmaps/$nm%s.mrc' % sfx, 'wb') as f: f.write(hdr.tobytes()); f.write(acc.astype(np.float32).tobytes())
PYEOF
        echo "$nm $D/leafmaps/$nm.mrc $D/leafmaps/${nm}_even_unfil.mrc $D/leafmaps/${nm}_odd_unfil.mrc $np_ $D/leafmaps/$nm.txt 1$J" >> $D/node_leaves.txt
      fi
      log "node $u: ${kind} ($np_ particles, states $mem) -> leaf"
    fi
  done < $D/groups.txt
  touch $D/node.done
}
[ -s $R/rows_r.txt ] || { if [ -n "${ROOTROWS:-}" ]; then awk '{print ($1 > 0) ? 1 : 0}' $ROOTROWS > $R/rows_r.txt; else awk -v n=$NR 'BEGIN {for (i = 1; i <= n; i++) print 1}' > $R/rows_r.txt; fi; }
echo "r 0 $SRC $R/rows_r.txt" > $R/level0.txt; : > $R/leaves.txt; lv=0
while [ -s $R/level$lv.txt ]; do
  while read u dep PP ROWS; do
    ( node $u $dep $PP $ROWS ) < /dev/null &
    while [ $(jobs -rp | wc -l) -ge ${PAR:-2} ]; do sleep 20; done
  done < $R/level$lv.txt
  wait
  : > $R/level$((lv + 1)).txt
  while read u dep PP ROWS; do
    [ -e $R/node_$u/node.done ] || { log "node $u did not finish"; exit 1; }
    cat $R/node_$u/children.txt >> $R/level$((lv + 1)).txt; cat $R/node_$u/node_leaves.txt >> $R/leaves.txt
  done < $R/level$lv.txt
  lv=$((lv + 1))
done
cd $R
$PY - <<PYEOF
import numpy as np
L = [l.split() for l in open('$R/leaves.txt') if l.strip()]; out = np.zeros($NR, int)
for i, l in enumerate(L, 1): lab = np.loadtxt(l[5], dtype=int); out[lab == int(l[6])] = i
np.savetxt('$R/final_labels.txt', out, fmt='%d'); open('$R/leaves_list.txt', 'w').write(''.join('%s %s %s %s\n' % (l[0], l[2], l[3], l[4]) for l in L))
print('%d leaves over %d particles: %s' % (len(L), (out > 0).sum(), ' '.join('%s:%s' % (l[0], l[4]) for l in L)))
PYEOF
score "rflex leaves" $R/final_labels.txt
[ -n "$TAG" ] && $PY $T/scorecard.py $R/final_labels.txt $GT --codes ${CODES:-swing:200:1} --tag "rflex leaves" 2>&1 | grep -v Warning
if [ -z "$TAG" ]; then awk '{printf "state_%02d %s %s\n", NR, $3, $4}' $R/leaves.txt > $R/leafview_spec.txt; $PY $T/make_view.py $R/leafmaps $R/leafview_spec.txt --box 160 > /dev/null 2>&1   # lp 10, box 160 for the readout
  CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1); $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $R/leaf_arm.txt $R/leafmaps 2>&1 | grep -v Warning | tail -1
  awk 'NR == FNR {n[FNR] = $1; c[FNR] = $5; next} !/^#/ {k = $1; sub(/.*state_/, "", k); sub(/\.mrc/, "", k); printf "   leaf %s (%s, %d ptcls): closed-arm %.2f open-arm %.2f\n", k, n[k + 0], c[k + 0], $2, $3}' $R/leaves.txt $R/leaf_arm.txt | sort -t: -k2 -r | head -12; fi
[ -s $R/final_tree/merge_tree.txt ] || $PY $T/state_merge.py $R/final_tree --list $R/leaves_list.txt --link complete --ncut 8,4,2 2>&1 | grep -v Warning | grep "states,\|^== " | head -4
log DONE
