#!/bin/bash
#SBATCH -J klem_rflex6split
#SBATCH -p csbdevel
#SBATCH -c 48
#SBATCH --mem=128G
#SBATCH -t 24:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# rflex6split = rflex4 with DATA-DERIVED children (resumable). Per node u (its particles at the poses of its parent's refinement): node mean, flex_pca
# K states, focus mask, focused K-state refine3D (as rflex4), then TWO SPLIT TESTS of the K refined states (splitbench16 on the even/odd halves
# of the node refinement, at the root project's consensus poses and at the node's own refined poses; a split needs both), cross-half EM <= EMLIM A, re-test only beyond the EM band with every group map
# docked (3-D translation fitted on the tested band) to its half consensus, view x defocus stratified shuffle null, every pair tested at the equal N of its smaller member; SBX = its extra options). A refined state lacking the node's signal (junk screen before the
# test: per-particle SSNR beyond the EM band < 0.35 of the node's; a low projection with signal is a partial unit, tested) is a junk leaf; every other state
# joins the supported group holding most of its particles (rflex6split_decide.py); one supported group = the node is a leaf; groups >= MINP recurse while depth < MAXD, the others are leaves (a leaf of several
# states: particle-weighted average of their half-maps). Final: per-row leaf labels, junk flags, scorecard (all leaves / non-junk), leaf tree.
# env RT, TAG, SRC, K (8), MINP (50000), MAXD (2), MSKDIAM (220), LPF (15), LP (10), MAXITS (8), NPROBE (flex default), CODES (swing:200:1),
#     SNAP (focus3), ROOTROWS, REUSE (a finished rflex run dir: a node with the same name, identical rows and a reused parent links its
#     mean / flex / refine outputs), EMLIM (20), EMMIN (1000), SBX, SIMD / POSES / CACHE / CTFABS (default sim$TAG,
#     lat_$TAG/ptcl3D_r0.txt, benchcache_full_$TAG, 0)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=${SNAP:-$HOME/SIMPLE_focus/snap/focus3}; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem
TAG=${TAG:-}; RT=${RT:?}; R=$K_/rflex_$RT; K=${K:-8}; MINP=${MINP:-50000}; MAXD=${MAXD:-2}; MD=${MSKDIAM:-220}; GT=$T/sim$TAG/gt_particles.txt
SRC=${SRC:-$T/proj$TAG/proj$TAG.simple}; mkdir -p $R && cd $R; log(){ echo "[$(date '+%F %T') rflex6split $RT] $*"; }
SIMD=${SIMD:-$T/sim$TAG}; CACHE=${CACHE:-$T/benchcache_full_$TAG}; POSES=${POSES:-$T/lat_$TAG/ptcl3D_r0.txt}; GTX=""; [ -n "$TAG" ] && GTX="--gt $GT,$GT"
st(){ simple_exec prg=print_project_field projfile=$1 oritype=ptcl3D mkdir=no 2>/dev/null | awk '{s = 0; for (i = 1; i <= NF; i++) if ($i ~ /^state=/) {split($i, a, "="); s = a[2] + 0}; print s}'; }
eo(){ simple_exec prg=print_project_field projfile=$1 oritype=ptcl3D mkdir=no 2>/dev/null | awk '{e = -1; for (i = 1; i <= NF; i++) if ($i ~ /^eo=/) {split($i, a, "="); e = a[2] + 0}; print e}'; }
score(){ [ -n "$TAG" ] || return 0; awk -v tag="$1" 'NR == FNR {if ($0 !~ /^#/) c[$1] = $2; next} {r++; s = $1; if (s <= 0) next; k = c[r]; n[s]++; if (k >= 200 && k <= 230) {x = k - 200; sw[s] += x; sq[s] += x * x; ns[s]++; T1 += x; T2 += x * x; NT++}}
  END {m = 0; for (s in n) if (ns[s] > 0) {mu[s] = sw[s] / ns[s]; ord[++m] = s; W += sq[s] - sw[s] ^ 2 / ns[s]}
    for (i = 1; i <= m; i++) for (j = i + 1; j <= m; j++) if (mu[ord[j]] < mu[ord[i]]) {t = ord[i]; ord[i] = ord[j]; ord[j] = t}
    printf "%s: %d swing particles, eta2 %.3f; per leaf mean +- sd swing (particles, swing share), by mean:", tag, NT, 1 - W / (T2 - T1 ^ 2 / NT); lo = 99; hi = -1
    for (i = 1; i <= m; i++) {s = ord[i]; sd = sqrt(sq[s] / ns[s] - mu[s] ^ 2); printf " %d:%.1f+-%.1f(%d,%.0f%%)", s, mu[s], sd, n[s], 100 * ns[s] / n[s]
      if (n[s] >= 200 && ns[s] / n[s] >= 0.5) {if (mu[s] < lo) lo = mu[s]; if (mu[s] > hi) hi = mu[s]}}
    printf "\n   %s extent (leaves >= 200 particles, >= 50%% swing): %.1f .. %.1f A of 0 .. 30\n", tag, lo, hi}' $GT $2; }
NR=$(simple_exec prg=print_project_field projfile=$SRC oritype=ptcl3D mkdir=no 2>/dev/null | grep -c .); log "project $SRC: $NR rows; K $K, MINP $MINP, MAXD $MAXD; split test poses $POSES, EM <= ${EMLIM:-20} A, ${SBX---dock 2 --junk_proj 0.75 --junk_ssnr 0.35}${REUSE:+; reuse $REUSE}"
[ -s $POSES ] || { POSES=$R/poses_root.txt; [ -s $POSES ] || simple_exec prg=print_project_field projfile=$SRC oritype=ptcl3D mkdir=no 2>/dev/null > $POSES; }
for f in $SIMD/all.mrcs ${CACHE}_A_c48.npy; do [ -s $f ] || { log "split test input $f missing"; exit 1; }; done
[ -s $R/queue.txt ] || { if [ -n "${ROOTROWS:-}" ]; then awk '{print ($1 > 0) ? 1 : 0}' $ROOTROWS > $R/rows_r.txt; else awk -v n=$NR 'BEGIN {for (i = 1; i <= n; i++) print 1}' > $R/rows_r.txt; fi; echo "r 0 $SRC $R/rows_r.txt" > $R/queue.txt; : > $R/leaves.txt; }
i=0
while true; do
  i=$((i + 1)); L=$(sed -n "${i}p" $R/queue.txt); [ -z "$L" ] && break
  set -- $L; u=$1; dep=$2; PP=$3; ROWS=$4; D=$R/node_$u; mkdir -p $D && cd $D
  [ -s $D/node.done ] && { log "node $u already done"; continue; }
  NN=$(awk '$1 > 0' $ROWS | wc -l)
  # reuse a finished node of REUSE with the same name and rows whose parent was reused (links only: nothing is written into REUSE)
  if [ -n "${REUSE:-}" ] && [ ! -e $D/reused ] && [ -d $REUSE/node_$u ] && cmp -s $ROWS $REUSE/rows_$u.txt; then pu=${u%?}
    if [ "$u" = r ] || [ -e $R/node_$pu/reused ]; then
      for x in mean $(cd $REUSE/node_$u && ls -d *_flex_pca *_refine3D); do ln -sfn $REUSE/node_$u/$x $D/$x; done
      cp -f $REUSE/node_$u/n.simple $REUSE/node_$u/mask.mrc $REUSE/node_$u/refined_states.txt $D/ && touch $D/reused && log "node $u: mean, flex and refinement reused from $REUSE/node_$u"; fi; fi
  # (0) node project at the parent's poses, members state 1
  if [ ! -s $D/n.simple ]; then cp -f $PP n.tmp.simple && simple_exec prg=selection projfile=n.tmp.simple oritype=ptcl3D infile=$ROWS mkdir=no > sel.log 2>&1 \
    && simple_exec prg=update_project projfile=n.tmp.simple qsys_name=local nparts=4 nthr=16 mkdir=no > up.log 2>&1 && mv n.tmp.simple n.simple || { log "node $u project FAILED"; exit 1; }; fi
  # (1) node mean
  if [ -z "$(ls $D/mean/*_reconstruct3D/recvol_state01.mrc 2>/dev/null)" ]; then mkdir -p $D/mean && cd $D/mean && cp -f $D/n.simple m.simple
    simple_exec prg=reconstruct3D projfile=m.simple pgrp=c1 mskdiam=$MD nstates=1 nparts=4 nthr=16 > rec.log 2>&1 || { log "node $u mean FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-150)"; exit 1; }; cd $D; fi
  MEAN=$(ls -td $D/mean/*_reconstruct3D | head -1)/recvol_state01.mrc
  # (2) flex on the node
  if [ -z "$(ls $D/*_flex_pca/flex_pca_coordinates.txt 2>/dev/null)" ]; then cp -f n.simple f.simple
    log "node $u (depth $dep, $NN particles): flex_pca, $K states, lp ${LPF:-15}"
    simple_exec prg=flex_pca projfile=f.simple vol1=$MEAN mskdiam=$MD lp=${LPF:-15} npreimages=$K rec_states=yes min_state_frac=0 umap=no ${NPROBE:+n_probe_iters=$NPROBE} nparts=4 nthr=16 > flex.log 2>&1 || { log "node $u flex_pca FAILED: $(grep -a -m2 'ERROR\|THROW' flex.log | cut -c1-200)"; exit 1; }; fi
  FD=$(ls -td $D/*_flex_pca | head -1); FP=$FD/f.simple; [ -s $FP ] || FP=$D/f.simple
  st $FP > $D/flex_states_raw.txt
  awk 'NR == FNR {if ($1 > 0) u[$1] = 1; next} FNR == 1 {m = 0; for (s in u) ord[++m] = s + 0; for (a = 1; a <= m; a++) for (b = a + 1; b <= m; b++) if (ord[b] < ord[a]) {t = ord[a]; ord[a] = ord[b]; ord[b] = t}; for (a = 1; a <= m; a++) map[ord[a]] = a}
    {print ($1 > 0) ? map[$1] : 0}' $D/flex_states_raw.txt $D/flex_states_raw.txt > $D/flex_states.txt
  mkdir -p $D/refs; : > $D/members.txt; VOLS=""; j=0
  for s in $(awk '$1 > 0' $D/flex_states_raw.txt | sort -n | uniq); do j=$((j + 1)); b=$D/refs/ref_state$(printf %02d $j); fm=$FD/recvol_state$(printf %02d $s)
    if [ -s $fm.mrc ]; then cp -f $fm.mrc $b.mrc && cp -f ${fm}_even.mrc ${b}_even.mrc && cp -f ${fm}_odd.mrc ${b}_odd.mrc || { log "node $u: flex state map $fm incomplete"; exit 1; }; HE=${fm}_even_unfil.mrc; HO=${fm}_odd_unfil.mrc
    else fs=$(printf %03d $s); cp -f $FD/flex_pca_state_$fs.mrc $b.mrc && cp -f $FD/flex_pca_even_state_$fs.mrc ${b}_even.mrc && cp -f $FD/flex_pca_odd_state_$fs.mrc ${b}_odd.mrc || { log "node $u: flex state map $s missing"; exit 1; }; HE=${b}_even.mrc; HO=${b}_odd.mrc; fi
    echo "$j $HE $HO $(awk -v s=$j '$1 == s' $D/flex_states.txt | wc -l)" >> $D/members.txt; VOLS="$VOLS vol$j=$b.mrc"; done
  KS=$j; log "node $u: flex states $(awk '{printf "%d ", $4}' $D/members.txt)"
  # (3) focus mask from the flex states
  [ -s $D/mask.mrc ] || $PY $T/focus_mask.py $D/mask.mrc --list $D/members.txt --mskdiam $MD 2>&1 | grep -v Warning | tail -1
  # (4) focused refinement of the flex states (hard assignments)
  if [ -z "$(ls $D/*_refine3D/recvol_state01_even_unfil.mrc 2>/dev/null)" ]; then cp -f n.simple r.simple
    simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$D/flex_states.txt mkdir=no > rsel.log 2>&1 || { log "node $u state selection FAILED"; exit 1; }
    log "node $u: focused refine3D, $KS states, lp ${LP:-10}, maxits ${MAXITS:-8}"
    simple_exec prg=refine3D projfile=r.simple $VOLS nstates=$KS refine=prob_neigh prob_neigh_mode=geom pgrp=c1 mskdiam=$MD lp=${LP:-10} maxits=${MAXITS:-8} center=no rec_backend=pcg pcg_mskfile=$D/mask.mrc pcg_focus=yes nspace=5000 box_crop=88 smpd_crop=2.66181803 nparts=4 nthr=16 > ref.log 2>&1 || { log "node $u refine3D FAILED: $(grep -a -m2 ERROR ref.log | cut -c1-150)"; exit 1; }; fi
  RD=$(ls -td $D/*_refine3D | head -1); RP=$RD/r.simple; [ -s $RP ] || RP=$D/r.simple
  [ -s $D/refined_states.txt ] || st $RP > $D/refined_states.txt
  [ -s $D/poses_node.txt ] || simple_exec prg=print_project_field projfile=$RP oritype=ptcl3D mkdir=no 2>/dev/null > $D/poses_node.txt
  [ -s $D/eo.txt ] || awk '{e = -1; for (i = 1; i <= NF; i++) if ($i ~ /^eo=/) {split($i, a, "="); e = a[2] + 0}; print e}' $D/poses_node.txt > $D/eo.txt
  score "node $u refined states" $D/refined_states.txt | tail -2
  # (5) split tests: the groups of refined states the halves support, at the root project's consensus poses and at the node's refined poses
  for tv in cons node; do TD=$D/test_$tv; mkdir -p $TD; [ -s $TD/summary.txt ] && continue; TP=$POSES; [ $tv = node ] && TP=$D/poses_node.txt
    awk 'NR == FNR {e[FNR] = $1; next} {print (e[FNR] == 0) ? $1 : 0}' $D/eo.txt $D/refined_states.txt > $TD/cand_A.txt
    awk 'NR == FNR {e[FNR] = $1; next} {print (e[FNR] == 1) ? $1 : 0}' $D/eo.txt $D/refined_states.txt > $TD/cand_B.txt
    OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 $PY -u $T/splitbench16.py --sim $SIMD,$SIMD --poses $TP,$TP --cand $TD/cand_A.txt,$TD/cand_B.txt --out $TD --cache $CACHE --ctfabs ${CTFABS:-0} --pair identity --softem 0 --arch 0 --strata 1 \
      --xem ${XEM:-3} --retest_beyond_em 1 --em_lim ${EMLIM:-20} --em_min ${EMMIN:-1000} ${SBX---dock 2 --junk_proj 0.75 --junk_ssnr 0.35} $GTX --summary $TD/summary.txt > $TD/bench.log 2>&1 &
  done; wait
  log "node $u: split tests of $KS refined states ($(awk '$1 > 0' $D/test_cons/cand_A.txt | wc -l) / $(awk '$1 > 0' $D/test_cons/cand_B.txt | wc -l) particles in the halves)"
  for tv in cons node; do [ -s $D/test_$tv/summary.txt ] || { log "node $u split test $tv FAILED: $(grep -a -m1 -A2 'Error\|Traceback' $D/test_$tv/bench.log | tail -1)"; exit 1; }
    log "node $u: split test ($tv poses) $(head -1 $D/test_$tv/summary.txt)"; grep -a "junk screen unit\|partial units" $D/test_$tv/bench.log | sed "s/^ */   node $u ($tv): /"; done
  $PY $T/rflex6split_decide.py $D/refined_states.txt $D/test_cons,$D/test_node $dep --minp $MINP --maxd $MAXD > $D/groups.txt || { log "node $u: decision FAILED"; exit 1; }
  g=0; fg=0; ng=$(grep -vc '^junk' $D/groups.txt)
  while read kind np_ mem jm; do
    g=$((g + 1))
    if [ "$kind" = family ]; then fg=$((fg + 1))
      awk -v m="$mem" 'BEGIN {n = split(m, a, "+"); for (i = 1; i <= n; i++) in_[a[i]] = 1} {print ($1 in in_) ? 1 : 0}' $D/refined_states.txt > $R/rows_${u}${fg}.txt
      grep -q "^${u}${fg} " $R/queue.txt || echo "${u}${fg} $((dep + 1)) $RP $R/rows_${u}${fg}.txt" >> $R/queue.txt; log "node $u: supported group $g ($np_ particles, states $mem; $jm) -> node ${u}${fg}"
    else
      J=""; [ "$kind" = junk ] && J=" junk"
      if [ "$mem" = "${mem#*+}" ]; then b=$RD/recvol_state$(printf %02d $mem); nm=${u}_s$mem
        grep -q "^$nm " $R/leaves.txt || echo "$nm $b.mrc ${b}_even_unfil.mrc ${b}_odd_unfil.mrc $np_ $D/refined_states.txt $mem$J" >> $R/leaves.txt
      else nm=${u}_g$g; [ $ng = 1 ] && nm=${u}_all; mkdir -p $D/leafmaps
        awk -v m="$mem" 'BEGIN {n = split(m, a, "+"); for (i = 1; i <= n; i++) in_[a[i]] = 1} {print ($1 in in_) ? 1 : 0}' $D/refined_states.txt > $D/leafmaps/$nm.txt
        $PY - <<PYEOF || { log "node $u: leaf map $nm FAILED"; exit 1; }
import numpy as np
lab = np.loadtxt('$D/refined_states.txt', dtype=int); mem = [int(x) for x in '$mem'.split('+')]; w = np.array([(lab == s).sum() for s in mem], float); w /= w.sum()
for sfx in ('', '_even_unfil', '_odd_unfil'):
    fs = ['$RD/recvol_state%02d%s.mrc' % (s, sfx) for s in mem]; hb = open(fs[0], 'rb').read(1024); h = np.frombuffer(hb, np.int32); off = 1024 + int(h[23]); nn = int(h[0]) * int(h[1]) * int(h[2])
    x = sum(wi * np.fromfile(f, np.float32, offset=off, count=nn).astype(np.float64) for wi, f in zip(w, fs)).astype(np.float32)
    hf = np.frombuffer(hb, np.float32).copy(); hf[19], hf[20], hf[21] = x.min(), x.max(), x.mean(); ext = open(fs[0], 'rb').read(off)[1024:]
    open('$D/leafmaps/$nm%s.mrc' % sfx, 'wb').write(hf.tobytes() + ext + x.tobytes())
PYEOF
        grep -q "^$nm " $R/leaves.txt || echo "$nm $D/leafmaps/$nm.mrc $D/leafmaps/${nm}_even_unfil.mrc $D/leafmaps/${nm}_odd_unfil.mrc $np_ $D/leafmaps/$nm.txt 1$J" >> $R/leaves.txt; fi
      log "node $u: group $g ($np_ particles, states $mem; $jm) -> $kind ($nm)"
    fi
  done < $D/groups.txt
  touch $D/node.done
done
cd $R
$PY - <<PYEOF
import numpy as np
L = [l.split() for l in open('$R/leaves.txt') if l.strip()]; out = np.zeros($NR, int); jk = np.zeros($NR, int)
for i, l in enumerate(L, 1): lab = np.loadtxt(l[5], dtype=int); out[lab == int(l[6])] = i; jk[lab == int(l[6])] = int(len(l) > 7 and l[7] == 'junk')
np.savetxt('$R/final_labels.txt', out, fmt='%d'); np.savetxt('$R/final_labels_nojunk.txt', out * (1 - jk), fmt='%d'); np.savetxt('$R/junk_rows.txt', jk, fmt='%d')
open('$R/leaves_list.txt', 'w').write(''.join('%s %s %s %s\n' % (l[0], l[2], l[3], l[4]) for l in L))
print('%d leaves (%d junk) over %d particles: %s' % (len(L), sum(len(l) > 7 for l in L), (out > 0).sum(), ' '.join('%s:%s%s' % (l[0], l[4], '(junk)' if len(l) > 7 else '') for l in L)))
if '$TAG':
    G = np.loadtxt('$GT', comments='#'); conf = np.zeros($NR, int); conf[G[:, 0].astype(int) - 1] = G[:, 1].astype(int); tj = conf == 0
    print('junk flag vs truth: %d of %d true junk particles flagged (recall %.3f); flagged rows %d, of them true junk %.3f' % ((tj & (jk > 0)).sum(), tj.sum(), (tj & (jk > 0)).sum() / max(tj.sum(), 1), jk.sum(), (tj & (jk > 0)).sum() / max(jk.sum(), 1)))
PYEOF
score "rflex6split leaves" $R/final_labels.txt
[ -n "$TAG" ] && $PY $T/scorecard.py $R/final_labels.txt $GT --codes ${CODES:-swing:200:1} --tag "rflex6split leaves" 2>&1 | grep -v Warning
[ -n "$TAG" ] && $PY $T/scorecard.py $R/final_labels_nojunk.txt $GT --codes ${CODES:-swing:200:1} --tag "rflex6split non-junk leaves" 2>&1 | grep -v Warning
[ -s $R/final_tree/merge_tree.txt ] || $PY $T/state_merge.py $R/final_tree --list $R/leaves_list.txt --link complete --ncut 8,4,2 2>&1 | grep -v Warning | grep "states,\|^== " | head -4
[ -n "$TAG" ] && $PY $T/leaf_report.py $R $GT 2>&1 | grep -v Warning
log DONE
