#!/bin/bash
#SBATCH -J toy10_tree_leaf_fb
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Pose feedback from the LEAVES: each level-3 grandchild (the balanced realigned tree, purity 0.71) seeds a 10-state
# refine3D_states (k-means 10 on its posterior means, started from its family's previous refinement project, local
# 20->8 A); the four refined pose tables are merged into one; the four grandchildren are rebuilt at those poses with
# the same recipe (balanced, focus, 12->8 A) and merged. env: L3 (_l3b, the grandchildren to refine), CB (their parent
# children prefix), BAL (0.5), MAXITS (6)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
BK=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$BK PATH=$BK/bin:$BK/scripts:$PATH TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; cd $W; log(){ echo "[$(date '+%F %T') lfb] $*"; }
L3=${L3:-_l3b}; CB=${CB:-$W/lat/klem_expand_realign1_bal05_fam}; RT0=$W/lat/realign_realign1/ptcl3D_realigned.txt; R=$W/lat/leaf_fb$L3; mkdir -p $R; RT=$R/ptcl3D_leaf_realigned.txt; TABS=""
for j in 1 2; do
  PS=$(ls -td $W/lat/state_refine_k10_fam$j/*_refine3D_states | head -1)/p.simple
  for s in 1 2; do
    KD=$W/lat/klem${L3}_fam${j}_$s; log "leaf $j.$s: k-means 10 -> refine3D_states from $PS"
    env KDIR=$KD K=10 NP=10 MAXITS=${MAXITS:-6} RTAG=${L3}_fam${j}_$s PROJSRC=$PS bash $W/run_state_refine_toy.sbatch 2>&1 | grep "k-means\|SEED\|REFINED\|poses\|FAILED" | sed "s/^/[leaf $j.$s] /"
    RD=$(ls -td $W/lat/state_refine_k10${L3}_fam${j}_$s/*_refine3D_states 2>/dev/null | head -1); [ -z "$RD" ] && { log "leaf $j.$s: no refinement output"; continue; }
    simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $R/refined_${j}_$s.txt; TABS="$TABS $R/refined_${j}_$s.txt"
  done
done
$PY - <<PYEOF
import numpy as np
rows = [l for l in open('$RT0') if l.strip() and not l.startswith('#')]; N = len(rows); n = 0
for f in '$TABS'.split():
    t = [l for l in open(f) if l.strip() and not l.startswith('#')]; assert len(t) == N
    for i, l in enumerate(t):
        kv = dict(x.split('=', 1) for x in l.split() if '=' in x)
        if float(kv.get('state', 0)) > 0: rows[i] = l.replace('state=%s' % kv['state'], 'state=1.00000000'); n += 1
open('$RT', 'w').write(''.join(rows)); print('leaf-realigned table: %d rows, %d from the four leaf refinements' % (N, n))
PYEOF
log "poses: $($PY $W/pose_error.py $RT $W/sim/all.txt --tag leaf-realigned)"
log "grandchildren rebuilt at the leaf-refined poses"
env CB=$CB TABLE=$RT BAL=${BAL:-0.5} L3=${L3}_fb bash $W/run_tree_level3b_toy.sbatch 2>&1 | grep "grandchild\|states: cross-half\|MERGED\|best-state\|Error\|Traceback\|DONE" | sed 's/^/[rebuilt] /'
log DONE
