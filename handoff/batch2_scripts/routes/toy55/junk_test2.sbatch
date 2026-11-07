#!/bin/bash
#SBATCH -J klem_junk
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# per-node junk test (junk_score.py) per state labelling in MLABS 'name=master_labels ..' (name K1*: one consensus state over the labelled rows;
# file truth:COL:JCOL: truth states, junk rows given random genuine labels) -> flag Jaccard across labellings -> optional REST 'name:label:stageBdir|proj.simple ..'
# rebuilds (KINDS: all / clean / junk / ctrl = equal-size random clean; label 'all' = every state) with arm readouts, FSC, consensus correlation. Poses: POSES 'A,B' or NSIMPLE (a node's n.simple, split by MRA/MRB).
# env H POSES NSIMPLE MRA MRB CP CTFABS EMLIM OUT TRUTH TCOL GUARD GUARDTAIL JTYPE REST
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; O=${OUT:?}; mkdir -p $O; MRA=${MRA:-${H:?}_A/meta_rows.txt}; MRB=${MRB:-${H}_B/meta_rows.txt}
if [ -n "${NSIMPLE:-}" ] && [ ! -s $O/poses_B.txt ]; then simple_exec prg=print_project_field projfile=$NSIMPLE oritype=ptcl3D mkdir=no 2>/dev/null > $O/pm.txt
  for h in A B; do mr=MR$h; awk 'NR == FNR {l[FNR] = $0; next} {print l[$1 + 0]}' $O/pm.txt ${!mr} > $O/poses_$h.txt; done; rm -f $O/pm.txt; fi
POSES=${POSES:-$O/poses_A.txt,$O/poses_B.txt}
NMS=""; for e in ${MLABS:?}; do nm=${e%%=*}; f=${e#*=}; D=$O/$nm; mkdir -p $D; NMS="$NMS $nm"
  [ -s $D/junk_flags_xd.txt ] || { $PY -c "
import numpy as np
f, nm = '$f', '$nm'
if ':' in f:
    tf, cn, jn = f.split(':'); hd = [l.split()[1:] for l in open(tf) if l.startswith('#')][0]; Tt = np.genfromtxt(tf, comments='#'); r = Tt[:, 0].astype(int) - 1
    L = np.zeros(r.max() + 1, int); L[r] = Tt[:, hd.index(cn)]; jk = np.zeros(len(L), bool); jk[r] = Tt[:, hd.index(jn)] > 0; L[jk] = np.random.default_rng(7).choice(L[~jk & (L > 0)], jk.sum())
else: L = np.loadtxt(f, dtype=int)
L = (L > 0).astype(int) if nm.startswith('K1') else L
for h, mr in (('A', '$MRA'), ('B', '$MRB')): r = np.loadtxt(mr).astype(int) - 1; np.savetxt('$D/lab_%s.txt' % h, L[r], fmt='%d')
open('$D/nmaster', 'w').write(str(len(L)))"
  $PY -u $T/junk_score.py --sim ${H}_A,${H}_B --poses ${POSES:?} --labels $D/lab_A.txt,$D/lab_B.txt --mr $MRA,$MRB --nmaster $(cat $D/nmaster) --out $D --cache ${CP:?} --ctfabs ${CTFABS:-1} \
    --em_lim ${EMLIM:-12} --tag $nm ${TRUTH:+--truth $TRUTH --tcol ${TCOL:-junk}} ${GUARD:+--guard $GUARD} ${GUARDTAIL:+--guardtail $GUARDTAIL} ${JTYPE:+--jtype $JTYPE} 2>&1 | grep -v Warning > $D/junk.log; }; cat $D/junk.log; done
$PY - <<PYEOF
import numpy as np
nms = '$NMS'.split()
import os
for cut in ('junk_flags', 'junk_flags_xd', 'junk_flags_xd3'):
    if not all(os.path.exists('$O/%s/%s.txt' % (k, cut)) for k in nms): continue
    F = {k: np.loadtxt('$O/%s/%s.txt' % (k, cut), dtype=int) for k in nms}
    for i, x in enumerate(nms):
        for y in nms[i + 1:]:
            m = (F[x] >= 0) & (F[y] >= 0); a, b = (F[x] == 1) & m, (F[y] == 1) & m
            print('%s: %s vs %s on %d shared node particles: flagged %d / %d, both %d, Jaccard %.3f' % (cut, x, y, m.sum(), a.sum(), b.sum(), (a & b).sum(), (a & b).sum() / max((a | b).sum(), 1)))
PYEOF
[ -n "${REST:-}" ] || exit 0
V=$A/rnasep_klem/view_rflex/junk_$(basename $O); mkdir -p $V; CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1); CONS=$T/maps160/toy_consensus.mrc; : > $O/spec_rest.txt
CUTS=${CUTS:-junk_flags junk_flags_xd}; C0=${CUTS%% *}
for e in $REST; do IFS=: read nm lab SB <<< "$e"; PJ=$SB; [ -d $SB ] && PJ=$SB/b.simple; LF=$(tr ' ' '\n' <<< "$MLABS" | grep "^$nm=" | cut -d= -f2)
  for cut in $CUTS; do for kind in ${KINDS:-all clean junk ctrl}; do [ $kind = all ] && [ $cut != $C0 ] && continue
  R=$O/rest_${nm}_s${lab}_all; [ $kind = all ] || R=$O/rest_${nm}_s${lab}${cut#junk_flags}_$kind; mkdir -p $R
  $PY -c "
import numpy as np
L = np.loadtxt('$LF', dtype=int); F = np.loadtxt('$O/$nm/$cut.txt', dtype=int); s = (L > 0) if '$lab' == 'all' else (L == int('$lab')); k = '$kind'
x = s if k == 'all' else s & (F == 0) if k == 'clean' else s & (F == 1)
if k == 'ctrl': c = np.flatnonzero(s & (F == 0)); x = np.zeros(len(L), bool); x[np.random.default_rng(5).choice(c, min((s & (F == 1)).sum(), len(c)), replace=False)] = True
np.savetxt('$R/rows.txt', x.astype(int), fmt='%d')"
  np_=$(awk '$1 > 0' $R/rows.txt | wc -l); echo "rebuild $nm state $lab $kind (${cut}): $np_ particles"; [ $np_ -ge 1500 ] || continue
  [ $kind = clean ] && cmp -s $R/rows.txt $O/rest_${nm}_s${lab}_all/rows.txt && { echo "  (no junk in this state: clean = all)"; continue; }
  if [ -z "$(ls $R/*_reconstruct3D/recvol_state01_even_unfil.mrc 2>/dev/null)" ] || ! cmp -s $R/rows.txt $R/rows_used.txt; then cd $R; rm -rf $R/*_reconstruct3D; cp -f $PJ r.simple
    simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$R/rows.txt mkdir=no > sel.log 2>&1 && simple_exec prg=update_project projfile=r.simple qsys_name=local nparts=3 nthr=8 mkdir=no > up.log 2>&1 \
      && simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=${MSKDIAM:-232} nstates=1 nparts=3 nthr=8 > rec.log 2>&1 || { echo "reconstruct $R FAILED: $(grep -a -m2 ERROR rec.log)"; cd $O; continue; }; cp $R/rows.txt $R/rows_used.txt; cd $O; fi
  RR=$(ls -td $R/*_reconstruct3D | head -1); echo "  $(basename $R): $(grep -a -m1 -i "0.143" $R/rec.log | tr -s ' ' | cut -c1-120)"
  grep -q "^$(basename $R) " $O/spec_rest.txt || echo "$(basename $R) $RR/recvol_state01_even_unfil.mrc $RR/recvol_state01_odd_unfil.mrc" >> $O/spec_rest.txt; done; done
  [ -d $SB ] && echo "stageB_${nm}_s$lab $SB/recvol_state$(printf %02d $lab)_even_unfil.mrc $SB/recvol_state$(printf %02d $lab)_odd_unfil.mrc" >> $O/spec_rest.txt; done
$PY $T/make_view.py $V $O/spec_rest.txt --box 160 > /dev/null 2>&1; $PY $T/arm_readout.py $CONS $CL $O/arm_rest.txt $(awk -v v=$V '{print v "/" $1 ".mrc"}' $O/spec_rest.txt) 2>&1 | grep -v Warning | tail -1; cat $O/arm_rest.txt
$PY - <<PYEOF
import numpy as np
def rd(p): h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n)
C = rd('$CONS').ravel()
for l in open('$O/spec_rest.txt'): nm = l.split()[0]; print('%s: correlation with the consensus %.3f' % (nm, np.corrcoef(C, rd('$V/%s.mrc' % nm).ravel())[0, 1]))
PYEOF
