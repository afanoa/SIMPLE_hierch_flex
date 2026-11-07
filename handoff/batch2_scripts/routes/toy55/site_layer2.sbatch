#!/bin/bash
#SBATCH -J klem_site_layer2
#SBATCH -p csbdevel
#SBATCH -c 24
#SBATCH --mem=160G
#SBATCH -t 08:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# FACTORIAL SITE LAYER on a state tree (Phase 1 + 2, resumable): discovery A over the tree's finished nodes (discover_sites.py: components
# of every node's cross-half variation, falloff-prior masks, de-duplicated, ranked by reproducibility) -> per site a 2-class site EM over
# ALL root particles (minus the root's junk leaves) at the root poses (minor_em3.py: cross-half references, per-shell reference scale,
# per-shell noise, learned mixing weights, the classes differ only inside the site mask; class 2 seeded with the discovering node's state
# that expresses the site most) -> class size from the posterior ranking (site_size.py: N* = argmax of the hold-out-band cross-half
# reproducibility of the top-N deviation) -> joint labels, pairwise independence (MI vs permutations), combinations (site_combine.py) ->
# SIMPLE reconstruct3D at the root poses per combination (>= MINC rows) and per site (class 1 / class 2) -> lp 10 unsharpened views in VIEW
# -> readouts: truth scorecard of the combination classes (score_tree_sites.py --labels; synthetic) or closed/open-arm + radial extent (real).
# env RT (tree: rnasep_klem/rflex_<RT>), SRC (root project = poses), H / MRA / MRB / CP / CTFABS (as rflex13s), MSKDIAM (220), TRUTH + TMAPS +
#     TMASKS (synthetic truth), TSIZE (1: discrete truth for the sizing), NSITE (6), EMLIM (15), ITERS (10), PARS (3 site EMs at a time), MINC (1000),
#     OUT, SITES (use this discovery dir), DISC_OPTS (extra discover_sites.py options, e.g. '--frac 0.1 --dice 0.5'), PICK (comma list of candidate names to classify instead of the ranked sites; name@ID seeds class 2 with state ID of the candidate's label file),
#     NODE (classify the particles of this node instead of the root; SRC = a project with that node's poses, e.g. its stage-B b.simple), TOPN (comma list: SIMPLE maps of every site's top-N by posterior), MOTION (1: a root cheap-children site also gets a K-class motion EM seeded by the ordered bins -> motion_coord.txt)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=4
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; RT=${RT:?}; R=$K_/rflex_$RT; SRC=${SRC:?}
H=${H:-$T/simrnasep692}; MRA=${MRA:-${H}_A/meta_rows.txt}; MRB=${MRB:-${H}_B/meta_rows.txt}; CP=${CP:-$T/benchcache_rnasep692}; CTFABS=${CTFABS:-1}; MD=${MSKDIAM:-220}
O=${OUT:-$K_/factorial/${RT}_sites}; V=$K_/view_rflex/factorial_$(basename $O); mkdir -p $O $V; cd $O; log(){ echo "[$(date '+%F %T') site layer $RT] $*"; }
ND=${NODE:-r}; NROWS=$R/rows_$ND.txt   # the node whose particles the site EMs classify (default the root); SRC must then carry that node's poses
NM=$(grep -c . $NROWS); CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1); CONS=$T/maps160/toy_consensus.mrc; MEAN=$(ls -td $R/node_$ND/mean/*_reconstruct3D | head -1)/recvol_state01.mrc
SD=${SITES:-$O/discovery}
# ---- 1 discovery (Phase 1)
if [ ! -s $SD/sites.txt ]; then log "discovery over $(ls -d $R/node_*/node.done 2>/dev/null | wc -l) finished nodes"
  $PY -u $T/discover_sites.py --tree $R --out $SD --mskdiam $MD --allcands 1 ${DISC_OPTS:-} ${TMASKS:+--truth $TMASKS} > $O/discovery.log 2>&1 || { log "discovery FAILED: $(tail -2 $O/discovery.log)"; exit 1; }; fi
grep -a "^SITE\|^TRUTH\|^PRECISION\|^RADIAL" $O/discovery.log 2>/dev/null
if [ -n "${PICK:-}" ]; then   # hand-picked candidates (names node:source:component from candidates.txt, comma-separated) replace the ranked sites
  mkdir -p $O/picked; : > $O/picked/sites.txt; j=0
  for cs in $(echo $PICK | tr ',' ' '); do c=${cs%@*}; L=$(awk -v c=$c '$2 == c' $SD/candidates.txt | head -1); [ -n "$L" ] || { log "PICK $c not among the candidates"; exit 1; }
    j=$((j + 1)); ln -sfn $(echo $L | awk '{print $10}') $O/picked/site_$(printf %02d $j).mrc
    if [ "$cs" != "$c" ]; then echo "$L" | awk -v j=$j -v sd=${cs#*@} '{$1 = j; $8 = sd; NF = 9; print}' >> $O/picked/sites.txt; else echo "$L" | awk -v j=$j '{$1 = j; NF = 9; print}' >> $O/picked/sites.txt; fi; done
  SD=$O/picked; log "picked sites: $(awk '{printf "%s (seed %s) ", $2, $8}' $SD/sites.txt)"; fi
# ---- 2 root poses per half; root rows minus the root's junk leaves
if [ ! -s $O/poses_B.txt ]; then simple_exec prg=print_project_field projfile=$SRC oritype=ptcl3D mkdir=no 2>/dev/null > $O/poses_master.txt
  awk 'NR == FNR {l[FNR] = $0; next} {print l[$1 + 0]}' $O/poses_master.txt $MRA > $O/poses_A.txt; awk 'NR == FNR {l[FNR] = $0; next} {print l[$1 + 0]}' $O/poses_master.txt $MRB > $O/poses_B.txt; rm -f $O/poses_master.txt; fi
$PY - <<PYEOF || { log "root rows FAILED"; exit 1; }
import numpy as np
use = np.loadtxt('$NROWS', dtype=int) > 0; st = np.loadtxt('$R/node_$ND/refined_f.txt', dtype=int); nj = 0
for l in open('$R/node_$ND/groups.txt'):
    w = l.split()
    if w and w[0] == 'junk':
        for s in w[2].split('+'): m = st == int(s); nj += m.sum(); use &= ~m
np.savetxt('$O/rows_use.txt', use.astype(int), fmt='%d'); print('root rows %d, minus %d in root junk leaves -> %d in the site EMs' % (use.sum() + nj, nj, use.sum()))
PYEOF
# ---- 3 per site: seeds -> site EM -> class size
site(){ local j=$1 D=$O/site_$1; mkdir -p $D; local L=$(awk -v j=$j '$1 == j' $SD/sites.txt); local labf=$(echo $L | awk '{print $7}') sid=$(echo $L | awk '{print $8}')
  $PY - <<PYEOF || return 1
import numpy as np
use = np.loadtxt('$O/rows_use.txt', dtype=int) > 0; seed = (np.loadtxt('$labf', dtype=int) == $sid) & use
for h, mr in (('A', '$MRA'), ('B', '$MRB')):
    r = np.loadtxt(mr).astype(int) - 1; np.savetxt('$D/start_%s.txt' % h, np.where(use[r], 1 + seed[r].astype(int), 0), fmt='%d')
tc = ''
if '${TRUTH:-}' and '${TMASKS:-}' and '${TSIZE:-1}' == '1':
    def rd(p):
        h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n)
    M = rd('$SD/site_%02d.mrc' % $j) > 0.5; ov = [((rd(p) > 0.5) & M).sum() / max((rd(p) > 0.5).sum(), 1) for p in '${TMASKS:-}'.split(',')]; col = int(np.argmax(ov))
    T = np.genfromtxt('${TRUTH:-}', comments='#'); v = (T[:, 3 + col] > 0).astype(float); ok = use & (T[:, 2] == 0); st = int(v[seed & ok].mean() > v[ok].mean())
    tc = '--truth ${TRUTH:-} --site %d --state %d' % (col, st); print('site $j ~ true site column %d (mask overlap %s), class 2 ~ state %d (seed %.3f vs all %.3f)' % (col, ' '.join('%.2f' % x for x in ov), st, v[seed & ok].mean(), v[ok].mean()))
open('$D/truthopt.txt', 'w').write(tc + '\n'); print('site $j seed: %d rows of state $sid of $labf' % seed.sum())
PYEOF
  if [ "${MOTION:-0}" = 1 ] && [ "$(basename $labf)" = children_1.txt ] && [ "$(echo $L | awk '{split($2, q, ":"); print q[1]}')" = r ]; then   # K-CLASS MOTION site: the root's ordered cheap bins
    [ -s $D/motion/post_soft_B.txt ] || { $PY -c "
import numpy as np
use = np.loadtxt('$O/rows_use.txt', dtype=int) > 0; ch = np.loadtxt('$labf', dtype=int)
for h, mr in (('A', '$MRA'), ('B', '$MRB')):
    r = np.loadtxt(mr).astype(int) - 1; np.savetxt('$D/mstart_%s.txt' % h, np.where(use[r] & (ch[r] > 0), ch[r], 0), fmt='%d')" && \
      $PY -u $T/minor_em3.py --sim ${H}_A,${H}_B --poses $O/poses_A.txt,$O/poses_B.txt --labels $D/mstart_A.txt,$D/mstart_B.txt --cache $CP --ctfabs $CTFABS --cons_ref $MEAN \
        --mask $SD/site_$(printf %02d $j).mrc --soft 1 --iters ${ITERS:-10} --em_lim ${EMLIM:-15} --post 1 --nomaps 1 --out $D/motion > $D/motion.log 2>&1 || echo "site $j motion EM FAILED: $(tail -2 $D/motion.log)"; }
    [ -s $D/motion/post_soft_B.txt ] && $PY -c "
import numpy as np
c = np.full($NM, np.nan)
for h, mr in (('A', '$MRA'), ('B', '$MRB')):
    r = np.loadtxt(mr).astype(int) - 1; P = np.loadtxt('$D/motion/post_soft_%s.txt' % h); s = np.loadtxt('$D/mstart_%s.txt' % h, dtype=int) > 0; c[r[s]] = (P[s] * np.arange(1, P.shape[1] + 1)).sum(1)
np.savetxt('$D/motion_coord.txt', np.nan_to_num(c, nan=-1.), fmt='%.4f'); print('site $j motion coordinate (posterior-mean bin of %d) written' % P.shape[1])"; fi
  [ -s $D/em/post_soft_B.txt ] || $PY -u $T/minor_em3.py --sim ${H}_A,${H}_B --poses $O/poses_A.txt,$O/poses_B.txt --labels $D/start_A.txt,$D/start_B.txt --cache $CP --ctfabs $CTFABS \
    --cons_ref $MEAN --mask $SD/site_$(printf %02d $j).mrc --soft 1 --iters ${ITERS:-10} --em_lim ${EMLIM:-15} --post 1 --nomaps 1 --out $D/em > $D/em.log 2>&1 || { echo "site $j EM FAILED: $(tail -2 $D/em.log)"; return 1; }
  [ -s $D/top_star.txt ] || $PY -u $T/site_size.py --em $D/em --start $D/start_A.txt,$D/start_B.txt --sim ${H}_A,${H}_B --poses $O/poses_A.txt,$O/poses_B.txt --mr $MRA,$MRB --nmaster $NM \
    --mask $SD/site_$(printf %02d $j).mrc --out $D --cache $CP --ctfabs $CTFABS --em_lim ${EMLIM:-15} $(cat $D/truthopt.txt) > $D/size.log 2>&1 || { echo "site $j sizing FAILED: $(tail -2 $D/size.log)"; return 1; }
  echo "== site $j"; grep -a '^start\|^iteration' $D/em.log | sed 's/deviation FSC.*//' | cut -c1-150 | tail -3; grep -a "^truth\|^N \|^N\*\|^top" $D/size.log; }
NS=$(awk '!/^#/' $SD/sites.txt | head -${NSITE:-6} | wc -l); log "$NS sites -> site EMs over the root ($(awk '$1 > 0' $O/rows_use.txt | wc -l) rows), ${PARS:-3} at a time"
for j in $(awk '!/^#/ {print $1}' $SD/sites.txt | head -${NSITE:-6}); do ( site $j > $O/site_$j.out 2>&1 ) & while [ $(jobs -rp | wc -l) -ge ${PARS:-3} ]; do sleep 20; done; done; wait
for j in $(awk '!/^#/ {print $1}' $SD/sites.txt | head -${NSITE:-6}); do cat $O/site_$j.out; done
# ---- 4 combinations
: > $O/spec_sites.txt; for j in $(awk '!/^#/ {print $1}' $SD/sites.txt | head -${NSITE:-6}); do [ -s $O/site_$j/top_star.txt ] && echo "s$j $O/site_$j/top_star.txt $O/site_$j/post_master.txt" >> $O/spec_sites.txt; done
$PY -u $T/site_combine.py --spec $O/spec_sites.txt --nmaster $NM --out $O/combine --rows $O/rows_use.txt --minc ${MINC:-1000} ${TRUTH:+--truth $TRUTH} 2>&1 | tee $O/combine.log
# ---- 5 SIMPLE maps at the root poses: combinations, then every site's two classes
rec(){ local d=$1 lab=$2; mkdir -p $d; [ -n "$(ls $d/*_reconstruct3D/recvol_state01_even_unfil.mrc 2>/dev/null)" ] && return 0; cd $d && cp -f $SRC r.simple
  simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$lab mkdir=no > sel.log 2>&1 && simple_exec prg=update_project projfile=r.simple qsys_name=local nparts=3 nthr=8 mkdir=no > up.log 2>&1 \
    && simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=$MD nstates=$(awk '$1 > 0' $lab | sort -n | tail -1) nparts=3 nthr=8 > rec.log 2>&1 || log "reconstruct $d FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-150)"; cd $O; }
rec $O/rec_combo $O/combine/combo_labels.txt; : > $O/spec_maps.txt; : > $O/leafmaps.txt; RC=$(ls -td $O/rec_combo/*_reconstruct3D | head -1)
while read st cid n_ bits; do b=$RC/recvol_state$(printf %02d $st); nm=combo_$(echo $bits | tr ' ' '-'); echo "$nm ${b}_even_unfil.mrc ${b}_odd_unfil.mrc" >> $O/spec_maps.txt; echo "$nm ${b}_even_unfil.mrc ${b}_odd_unfil.mrc" >> $O/leafmaps.txt; done < <(grep -v '^#' $O/combine/combos.txt)
for j in $(awk '{print substr($1, 2)}' $O/spec_sites.txt); do rec $O/rec_site_$j $O/combine/site_s${j}_labels.txt; RS=$(ls -td $O/rec_site_$j/*_reconstruct3D | head -1)
  for c in 1 2; do echo "site${j}_class$c $RS/recvol_state0${c}_even_unfil.mrc $RS/recvol_state0${c}_odd_unfil.mrc" >> $O/spec_maps.txt; done; done
if [ -n "${TOPN:-}" ]; then   # per site: SIMPLE maps of the top-N rows by posterior (nested ranking readout, site_r22-style)
  for j in $(awk '{print substr($1, 2)}' $O/spec_sites.txt); do for N in $(echo $TOPN | tr ',' ' '); do
    [ -s $O/site_$j/top$N.txt ] || $PY -c "import numpy as np; p = np.loadtxt('$O/site_$j/post_master.txt'); u = np.loadtxt('$O/rows_use.txt') > 0; p[~u] = -1; o = np.argsort(-p); l = np.zeros(len(p), int); l[o[:$N]] = 1; np.savetxt('$O/site_$j/top$N.txt', l, fmt='%d')"
    rec $O/rec_site${j}_top$N $O/site_$j/top$N.txt; RT_=$(ls -td $O/rec_site${j}_top$N/*_reconstruct3D | head -1); echo "site${j}_top$N $RT_/recvol_state01_even_unfil.mrc $RT_/recvol_state01_odd_unfil.mrc" >> $O/spec_maps.txt; done; done; fi
$PY $T/make_view.py $V $O/spec_maps.txt --box 160 > /dev/null 2>&1; log "maps (lp 10, unsharpened, full maps only): $V"
# ---- 6 readouts
if [ -n "${TRUTH:-}" ]; then
  $PY $T/score_tree_sites.py --tree $R --truth $TRUTH --labels $O/combine/combo_labels.txt --leafmaps $O/leafmaps.txt ${TMAPS:+--maps $TMAPS} ${TMASKS:+--masks $TMASKS} --json $O/score_combos.json 2>&1 | grep -v Warning | tee $O/score_combos.txt
else
  $PY $T/arm_readout.py $CONS $CL $O/arm.txt $(ls $V/*.mrc) 2>&1 | grep -v Warning | tail -1
  $PY $T/radial_extent.py $CONS $CL --maps $(grep -v '^#' $O/arm.txt | sort -k2,2gr | head -8 | awk '{n = $1; sub(/.*\//, "", n); sub(/\.mrc$/, "", n); printf "%s:%s ", n, $1}') 2>&1 | grep -v Warning; fi
log DONE
