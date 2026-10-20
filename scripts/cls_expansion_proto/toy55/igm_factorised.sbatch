#!/bin/bash
#SBATCH -J igm_factorised
#SBATCH -p csbdevel
#SBATCH -c 48
#SBATCH --mem=160G
#SBATCH -t 08:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Truth-free factorised route on a hinged-arm set, ONE job (truth only in the readouts): (1) junk = the data-derived tree's
# root children that its root board assigns mostly to the free components (both halves, mapped to full-set rows);
# (2) a whole-particle board on the cleaned set (7x7 + 1 free, halves written) -> motion-pattern regions (derive_regions2);
# (3) one focused board per region (balanced N x N lattice + 1 free, random start, cleaned set, refined full-set poses),
# BPAR at a time; (4) factorised_readout.py. env MODE (armsigm1), TAG (dtree), NLAT (5), KREG (8), MAXREG (8), BPAR (3),
# V1BOARDS (0 = stop after the parent board: the v2-v4 stages need only the cleaned subset and this board)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; K=$A/rnasep_klem; MODE=${MODE:-armsigm1}; TAG=${TAG:-dtree}
O=$W/factorised_$MODE; mkdir -p $O; cd $O; export FB=$K/f692e/klem2_bin SIMPLE_QSYS=local; PS=$K/synth_$MODE/src.simple; NL=${NLAT:-5}
log(){ echo "[$(date '+%F %T') factorised $MODE] $*"; }
# ---- (1) junk from the tree root (no truth): root children whose particles the root board puts mostly on its free components
$PY - <<PYEOF
import numpy as np
N = sum(1 for l in open('$W/sim$MODE/gt_particles.txt') if l.strip() and not l.startswith('#')); junk = np.zeros(N + 1, bool)
for H, off in (('A', -1), ('B', 0)):
    r = '$K/synth_${MODE}_%s/$TAG/r' % H; ch = np.loadtxt(r + '/children.txt', dtype=int); cfg = dict(l.split()[:2] for l in open(r + '/board/board.txt') if l.strip())
    ng = int(cfg['n1']) * int(cfg['n2']); L = np.loadtxt(r + '/board/klem_latent.txt', comments='#'); free = np.zeros(len(ch) + 1, bool); free[L[:, 0].astype(int)] = L[:, 3] > ng
    for c in sorted(set(ch[ch > 0])):
        rows = np.where(ch == c)[0] + 1; f = free[rows].mean(); tag = 'JUNK' if f > 0.5 else ''
        print('half %s root child %d: %d particles, %.2f on the free components %s' % (H, c, len(rows), f, tag))
        if f > 0.5: junk[2 * rows + off] = True
G = np.loadtxt('$W/sim$MODE/gt_particles.txt', comments='#'); tj = np.zeros(N + 1, bool); tj[G[:, 0].astype(int)] = G[:, 1] == 0
keep = np.where(~junk[1:])[0] + 1; np.savetxt('$O/clean_subset.txt', keep, fmt='%d')
print('junk removed: %d particles; truth readout: captures %.3f of the junk, purity %.3f; cleaned set %d particles' % (junk.sum(), (junk & tj).sum() / tj.sum(), (junk & tj).sum() / max(junk.sum(), 1), len(keep)))
PYEOF
[ -s $O/clean_subset.txt ] || { log "junk step FAILED"; exit 1; }
# ---- (2) whole-particle board on the cleaned set -> motion-pattern regions
log "whole-particle board on the cleaned set (7x7 + 1 free, halves)"
[ -s $O/parent/klem_latent.txt ] || SERVER_LIST= NPARTS= NTHR=48 N1=7 N2=7 NFREE=1 INIT=random SEED=1 WRITE_HALVES=1 bash $W/fboard.sh $O/parent $PS $O/clean_subset.txt "20 20 20 14 14 14 10 10 10" "1.5 1.5 1.2 1.0 1.0 0.8 0.7 0.7 0.7" > $O/parent.out 2>&1 || { log "parent board FAILED: $(tail -2 $O/parent.out)"; exit 1; }
[ "${V1BOARDS:-1}" = 0 ] && { log "V1BOARDS=0: cleaned subset and parent board ready"; log DONE; exit 0; }
[ -s $O/regions/region_1.mrc ] || $PY $W/derive_regions2.py $O/parent $O/regions --k ${KREG:-8} --maxreg ${MAXREG:-8} 2>&1 | grep -v Warning | tail -12
NR=$(ls $O/regions/region_*.mrc 2>/dev/null | wc -l); log "$NR data-derived regions"; [ $NR -gt 0 ] || exit 1
# ---- (3) one focused board per region, BPAR at a time (16 threads each)
for f in $(ls $O/regions/region_*.mrc); do r=$(basename $f .mrc); r=${r#region_}; D=$O/arm_$r; mkdir -p $D
  while [ $(jobs -rp | wc -l) -ge ${BPAR:-3} ]; do sleep 10; done
  [ -s $D/board/klem_latent.txt ] || ( SERVER_LIST= NPARTS= NTHR=16 N1=$NL N2=$NL NFREE=1 CTFABS=1 PSCALE=1 FOCUSMASK=$f INIT=random SEED=$r bash $W/fboard.sh $D/board $PS $O/clean_subset.txt "20 14 12 10 10 8 8 8" "1.5 1.2 1.0 1.0 0.8 0.8 0.7 0.7" > $D/board.out 2>&1 || echo "region $r board FAILED" ) &
done; wait; log "region boards done"
# ---- (4) readout (truth)
$PY $W/factorised_readout.py $O $W/sim$MODE/gt_particles.txt $W/truth_$MODE.txt $W/maps160/synth_$MODE/masks 2>&1 | grep -v Warning
log DONE
