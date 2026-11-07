#!/bin/bash
#SBATCH -J klem_splittest
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# data-derived split of the ROOT of a tree pair (env MODE e.g. null_50k, K candidate groups (6), BOARD (parent)): each tree's
# board cells -> K candidate groups, every group reconstructed from the COMMON consensus poses (src.simple), then the
# cross-tree merge (calibrated noise null) decides how many groups are distinct; purity of the surviving groups vs truth
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; MODE=${MODE:?}; K=${K:-6}; BOARD=${BOARD:-parent}; T=split${BOARD}_k$K${SHUF:+_shuf$SHUF}${NEQ:+_neq$NEQ}   # NEQ=min|N: every group reconstructed from the same particle count (no size-dependent filter); SHUF=<seed>: the candidate labels shuffled among the on-grid particles (same group sizes, no real difference)
log(){ echo "[$(date '+%F %T') splittest $MODE] $*"; }
for H in A B; do R=$A/rnasep_klem/synth_${MODE}_$H; D=$R/fam0/region$T/level1; mkdir -p $D/rec
  NROWS=$($PY -c "import numpy as np; print(len(np.load('$R/meta_toy.npz')['flags']))")
  [ -s $D/maps/recvol_state01_even.mrc ] && { echo "maps of $H exist: merge only"; continue; }
  $PY $W/split_candidates.py $R/$BOARD $NROWS $K $D 2>&1 | grep -v Warning
  [ -n "${SHUF:-}" ] && $PY -c "import numpy as np; s = np.loadtxt('$D/states.txt').astype(int); on = np.where(s > 0)[0]; s[on] = np.random.default_rng($SHUF).permutation(s[on]); np.savetxt('$D/states.txt', s, fmt='%d'); print('shuffled', np.bincount(s).tolist())"
  [ -n "${NEQ:-}" ] && $PY $W/neq_subsample.py $D/states.txt $D/states_full.txt $NEQ
  ( cd $D/rec && cp -f $R/src.simple p.simple && simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$D/states.txt mkdir=no > sel.log 2>&1 && \
    simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=$K nparts=4 nthr=8 > rec.log 2>&1 && \
    RR=$(ls -td ./*_reconstruct3D | head -1) && mkdir -p $D/maps && for s in $(seq -f "%02g" 1 $K); do for h in "" _even _odd; do cp -f $RR/recvol_state$s$h.mrc $D/maps/ 2>/dev/null; done; done && echo "rec $H done" || echo "rec $H FAILED: $(grep -a -m2 'ERROR\|THROW' rec.log sel.log | cut -c1-160)" ) &
done; wait; log "reconstructions done"
RA=$A/rnasep_klem/synth_${MODE}_A; RB=$A/rnasep_klem/synth_${MODE}_B
for V in ${VS:-0 1}; do XB=${BANDLO:-40}; echo "== cross-tree merge of the $K root groups (NUIS=$V, band from $XB A, noise-null threshold 99 %)"
  env FAMS=0 LEAFGLOB=region$T MASK=sphere NPROC=16 MINREC=${MINREC:-300} NUIS=$V BANDLO=$XB NULLCAL=0.99 ARCH=${ARCH:-0} THRMIN=${THRMIN:-0} $PY $W/merge_cross_tree.py $RA $RB $W/maps160/toy_consensus.mrc $RA/fam0/region$T/merge_nuis${V}_b$XB${ARCH:+_a$ARCH} 2>&1 | grep -v Warning | grep -a "leaf states\|matched\|pairs at the start\|pair statistics\|null of the statistic\|merges ->\|remaining pairs\|archetype\|Traceback\|Error"
  for H in A B; do R=$A/rnasep_klem/synth_${MODE}_$H
$PY - <<PYEOF
import numpy as np
z = np.load('$RA/fam0/region$T/merge_nuis${V}_b${XB}${ARCH:+_a$ARCH}_$H.npz'); names, labs = list(z['names']), z['labels']
import os; sf = '$R/fam0/region$T/level1/states_full.txt'; st = np.loadtxt(sf if os.path.isfile(sf) else '$R/fam0/region$T/level1/states.txt').astype(int); conf = np.loadtxt('$W/sim${MODE}_$H/gt_particles.txt', comments='#')[:, 1].astype(int)
g = np.zeros(len(st), int)
for nm, lb in zip(names, labs): g[st == int(nm.split('_s')[-1])] = lb
ok = (g > 0) & (conf > 0); out = []
for k in np.unique(g[ok]):
    m = ok & (g == k); c = np.bincount(conf[m]); out.append('%d:(n=%d, top conformer %d %.2f)' % (k, m.sum(), c.argmax(), c.max() / m.sum()))
print('  half $H: %d merged groups: %s' % (len(np.unique(g[ok])), ' '.join(out)))
PYEOF
  done; done
log DONE
