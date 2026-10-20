#!/bin/bash
#SBATCH -J toy10_tree_realign
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# One-pass realignment of the whole tree. The tree has already refined every on-grid particle once (its family's 10-state
# refinement); those poses are merged into ONE table (family members from their family's refined table, the rest
# consensus) and the tree is regrown on it from the parent down: parent board (junk comps) -> families (spectral
# bisection) -> child boards at the realigned poses (the feedback recipe) -> merged purity, against the 0.62 of the
# feedback boards grown on the consensus-pose families. LOOP=1 then takes a second pass: k-means 10 on each regrown child
# -> refine3D_states per family (started from that family's previous refinement project, matched by membership overlap)
# -> feedback boards at the second-pass poses -> merged purity.
# env: FB (prefix of the feedback pose dirs, lat/fb_fam), E0 (old expansion dir with the family subsets), RTAG (_realign1),
#      NFAM (2), LOOP (1), MAXITS (6), BAL (balanced E-step lambda for the parent and the children), PM (children pose marginalisation, 3; 0 = off)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
BK=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$BK PATH=$BK/bin:$BK/scripts:$PATH TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; cd $W; log(){ echo "[$(date '+%F %T') tr] $*"; }
RTAG=${RTAG:-_realign1}; E0=${E0:-$W/lat/expand2b}; FB=${FB:-$W/lat/fb_fam}; NFAM=${NFAM:-2}; R=$W/lat/realign$RTAG; mkdir -p $R; RT=$R/ptcl3D_realigned.txt
# 1. one table of realigned poses
$PY - <<PYEOF
import numpy as np
rows = [l for l in open('$W/lat/ptcl3D_r0.txt') if l.strip() and not l.startswith('#')]; N = len(rows); src = np.zeros(N, int)
for j in range($NFAM):
    fam = np.loadtxt('$E0/family%d_subset.txt' % (j + 1), dtype=int); t = [l for l in open('$FB%d/refined_state1.txt' % (j + 1)) if l.strip() and not l.startswith('#')]
    assert len(t) == N, (len(t), N)
    for p in fam: rows[p - 1] = t[p - 1]; src[p - 1] = j + 1
open('$RT', 'w').write(''.join(rows))
print('realigned table: %d rows; %s; %d rows keep the consensus pose' % (N, ', '.join('%d from family %d' % ((src == j + 1).sum(), j + 1) for j in range($NFAM)), (src == 0).sum()))
PYEOF
log "poses: $($PY $W/pose_error.py $RT $W/sim/all.txt --tag realigned)"
# 2. the parent regrown at the realigned poses
log "parent board at the realigned poses"
env KTAG=$RTAG INIT=pca NSUB=0 NFREE=2 TABLE=$RT ${BAL:+BAL=$BAL} ${MSKD:+MSKD=$MSKD} bash $W/run_klem.sbatch 2>&1 | grep "final:\|states: cross-half\|Error\|Traceback" | sed "s/^/[parent] /"
$PY $W/state_coverage.py $W/lat/klem$RTAG/klem_latent.txt $W/sim/gt_particles.txt --argmax --tag "parent at realigned poses" | tail -1
$PY $W/state_coverage.py $W/lat/klem_junk/klem_latent.txt $W/sim/gt_particles.txt --argmax --tag "parent at consensus poses (klem_junk)" | tail -1
# 3. families and child boards at the realigned poses (feedback recipe: no free comps, pose marginalisation, devw, pnoise)
E=$W/lat/expand$RTAG; mkdir -p $E; $PY $W/expand_families.py split $W/lat/klem$RTAG $E --nfam $NFAM
$PY - <<PYEOF
import numpy as np
new = [set(np.loadtxt('$E/family%d_subset.txt' % (j + 1), dtype=int)) for j in range($NFAM)]; old = [set(np.loadtxt('$E0/family%d_subset.txt' % (j + 1), dtype=int)) for j in range($NFAM)]
print('family overlap new x old (Jaccard): ' + ' '.join('n%d/o%d:%.2f' % (a + 1, b + 1, len(new[a] & old[b]) / len(new[a] | old[b])) for a in range($NFAM) for b in range($NFAM)))
PYEOF
CH=""
for j in $(seq 1 $NFAM); do
  log "child board for family $j at the realigned poses"
  env KTAG=_expand${RTAG}_fam$j INIT=pca NSUB=0 SUBSET=$E/family${j}_subset.txt NFREE=0 TABLE=$RT PMARG=${PM:-3} DEVW=0.05 PNOISE=1 ${BAL:+BAL=$BAL} ${MSKD:+MSKD=$MSKD} bash $W/run_klem.sbatch 2>&1 | grep "subset\|final:\|states: cross-half\|Error\|Traceback" | sed "s/^/[fam $j] /"
  CH="$CH $W/lat/klem_expand${RTAG}_fam$j"
done
$PY $W/expand_families.py merge $E $CH --gt $W/sim/gt_particles.txt --minrep 0.5 | sed 's/^/[regrown tree] /'
echo "--- feedback boards on the consensus-pose families, for reference:"; $PY $W/expand_families.py merge $E0 $W/lat/klem_fb_fam1 $W/lat/klem_fb_fam2 --gt $W/sim/gt_particles.txt --minrep 0.5 | grep MERGED | sed 's/^/[old tree] /'
# 4. second pass: refine the regrown families' states, feed the poses back
if [ "${LOOP:-1}" = 1 ]; then
  CH2=""
  for j in $(seq 1 $NFAM); do
    OLD=$($PY -c "
import numpy as np
new = set(np.loadtxt('$E/family${j}_subset.txt', dtype=int)); print(1 + max(range($NFAM), key=lambda k: len(new & set(np.loadtxt('$E0/family%d_subset.txt' % (k + 1), dtype=int)))))")
    PS=$(ls -td $W/lat/state_refine_k10_fam$OLD/*_refine3D_states | head -1)/p.simple
    log "family $j (old family $OLD): k-means 10 on the regrown child board -> refine3D_states from $PS"
    env KDIR=$W/lat/klem_expand${RTAG}_fam$j K=10 NP=10 MAXITS=${MAXITS:-6} RTAG=${RTAG}_fam$j PROJSRC=$PS bash $W/run_state_refine_toy.sbatch 2>&1 | grep "k-means\|SEED\|REFINED\|poses\|FAILED\|purity\|owning" | sed "s/^/[fam $j refine] /"
    RD=$(ls -td $W/lat/state_refine_k10${RTAG}_fam$j/*_refine3D_states 2>/dev/null | head -1); [ -z "$RD" ] && { log "family $j: no refinement output, skipping the feedback"; continue; }
    T=$W/lat/fb${RTAG}_fam$j; mkdir -p $T
    simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $T/refined_ptcl3D.txt
    awk '{ gsub(/state=[0-9.]+/, "state=1"); print }' $T/refined_ptcl3D.txt > $T/refined_state1.txt
    log "family $j: second-pass poses: $($PY $W/pose_error.py $T/refined_ptcl3D.txt $W/sim/all.txt --tag poses)"
    env KTAG=_fb${RTAG}_fam$j INIT=pca NSUB=0 SUBSET=$E/family${j}_subset.txt NFREE=0 TABLE=$T/refined_state1.txt PMARG=3 DEVW=0.05 PNOISE=1 bash $W/run_klem.sbatch 2>&1 | grep "subset\|final:\|states: cross-half\|Error\|Traceback" | sed "s/^/[fam $j fb2] /"
    CH2="$CH2 $W/lat/klem_fb${RTAG}_fam$j"
  done
  [ -n "$CH2" ] && $PY $W/expand_families.py merge $E $CH2 --gt $W/sim/gt_particles.txt --minrep 0.5 | sed 's/^/[second pass] /'
fi
log DONE
