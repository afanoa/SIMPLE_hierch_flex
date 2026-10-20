#!/bin/bash
#SBATCH -J igm_factorised3
#SBATCH -p csbdevel
#SBATCH -c 48
#SBATCH --mem=160G
#SBATCH -t 10:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Truth-free factorised route, region step v3 (truth only in readouts), on top of igm_factorised2.sbatch's periphery,
# regions and boards: (1) residual periphery = periphery minus the found regions (dilated); (2) a board focused on the
# residual (7x7 + 1 free, 14 -> 8 A, halves) -> regions inside the residual (derive_regions2 --within) -> one focused board
# each; (3) merge all candidate regions whose board coordinates correlate (region_merge.py, THR); (4) one balanced board
# per merged region; (5) factorised_readout.py (residual stage and merged stage). env MODE (armsigm1), THR (gap), NLAT (5), BPAR (3),
# SKIPMERGE (1 = stop after the residual stage; igm_factorised4.sbatch does the merge with multi-start boards)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; KL=$A/rnasep_klem; MODE=${MODE:-armsigm1}
O1=$W/factorised_$MODE; O2=$W/factorised2_$MODE; O=$W/factorised3_$MODE; mkdir -p $O/resid $O/merged; cd $O; export FB=$KL/f692e/klem2_bin SIMPLE_QSYS=local
PS=$KL/synth_$MODE/src.simple; SUB=$O1/clean_subset.txt; NL=${NLAT:-5}; MK=$W/maps160/synth_$MODE/masks; log(){ echo "[$(date '+%F %T') factorised3 $MODE] $*"; }
board(){ local f=$1 D=$2 s=$3; mkdir -p $D; [ -s $D/board/klem_latent.txt ] || SERVER_LIST= NPARTS= NTHR=16 N1=$NL N2=$NL NFREE=1 CTFABS=1 PSCALE=1 FOCUSMASK=$f INIT=random SEED=$s bash $W/fboard.sh $D/board $PS $SUB "20 14 12 10 10 8 8 8" "1.5 1.2 1.0 1.0 0.8 0.8 0.7 0.7" > $D/board.out 2>&1 || echo "board FAILED $D"; }
# ---- (1) residual periphery
$PY - <<PYEOF
import sys, glob, numpy as np; sys.path.insert(0, '$W'); from model_class import read_mrc, write_mrc
from scipy.ndimage import binary_dilation, gaussian_filter
P = read_mrc('$O2/periphery.mrc'); U = np.zeros(P.shape, bool)
for f in glob.glob('$O2/regions/region_*.mrc'): U |= read_mrc(f) > 0.5
U = binary_dilation(U, iterations=2); R = gaussian_filter(((P > 0.5) & ~U).astype(np.float32), 1.5); write_mrc('$O/residual.mrc', R, 1.464)
v3 = 1.464 ** 3 / 1000.; print('residual periphery: %.0f nm^3 of %.0f nm^3 (found regions cover %.0f nm^3)' % ((R > 0.5).sum() * v3, (P > 0.5).sum() * v3, (U & (P > 0.5)).sum() * v3))
B = P.shape[0]; c = B / 2; zz, yy, xx = np.mgrid[:B, :B, :B]; r = np.sqrt((xx - c) ** 2 + (yy - c) ** 2 + (zz - c) ** 2) * 1.464
sh = [((R > 0.5) & (read_mrc(m) > 0.5) & (r > 55)).sum() / max(((read_mrc(m) > 0.5) & (r > 55) & (P > 0.5)).sum(), 1) for m in sorted(glob.glob('$MK/arm_*.mrc'))]
print('   truth check: share of each arm zone (beyond 55 A, in the periphery) left in the residual %s' % ' '.join('%.2f' % v for v in sh))
PYEOF
# ---- (2) residual board -> residual regions -> their boards
log "board focused on the residual periphery (7x7 + 1 free, 14 -> 8 A)"
[ -s $O/resid/parent/klem_latent.txt ] || SERVER_LIST= NPARTS= NTHR=48 N1=7 N2=7 NFREE=1 CTFABS=1 PSCALE=1 FOCUSMASK=$O/residual.mrc INIT=random SEED=11 WRITE_HALVES=1 bash $W/fboard.sh $O/resid/parent $PS $SUB "14 14 12 12 10 10 8 8 8" "1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7" > $O/resid/parent.out 2>&1 || { log "residual board FAILED: $(tail -2 $O/resid/parent.out)"; exit 1; }
[ -s $O/resid/regions/region_1.mrc ] || $PY $W/derive_regions2.py $O/resid/parent $O/resid/regions --k ${KRES:-4} --maxreg ${MAXRES:-4} --within $O/residual.mrc 2>&1 | grep -v Warning | tail -8
NR=$(ls $O/resid/regions/region_*.mrc 2>/dev/null | wc -l); log "$NR residual regions"
for f in $(ls $O/resid/regions/region_*.mrc 2>/dev/null); do r=$(basename $f .mrc); r=${r#region_}; while [ $(jobs -rp | wc -l) -ge ${BPAR:-3} ]; do sleep 10; done; board $f $O/resid/arm_$r $((20 + r)) & done; wait
echo "== residual stage readout (truth)"; $PY $W/factorised_readout.py $O/resid $W/sim$MODE/gt_particles.txt $W/truth_$MODE.txt $MK 2>&1 | grep -v Warning
[ "${SKIPMERGE:-0}" = 1 ] && { log "SKIPMERGE=1: residual stage done"; log DONE; exit 0; }
# ---- (3) merge all candidates by board-coordinate correlation
CANDS=""; for f in $(ls $O2/regions/region_*.mrc); do r=$(basename $f .mrc); r=${r#region_}; [ -s $O2/arm_$r/board/klem_latent.txt ] && CANDS="$CANDS $f:$O2/arm_$r/board"; done
for f in $(ls $O/resid/regions/region_*.mrc 2>/dev/null); do r=$(basename $f .mrc); r=${r#region_}; [ -s $O/resid/arm_$r/board/klem_latent.txt ] && CANDS="$CANDS $f:$O/resid/arm_$r/board"; done
[ -s $O/merged/regions/region_1.mrc ] || $PY $W/region_merge.py $O/merged/regions ${THR:-gap} $CANDS 2>&1 | grep -v Warning
# ---- (4) one balanced board per merged region
for f in $(ls $O/merged/regions/region_*.mrc); do r=$(basename $f .mrc); r=${r#region_}; while [ $(jobs -rp | wc -l) -ge ${BPAR:-3} ]; do sleep 10; done; board $f $O/merged/arm_$r $((40 + r)) & done; wait
log "merged-region boards done"
echo "== merged stage readout (truth)"; $PY $W/factorised_readout.py $O/merged $W/sim$MODE/gt_particles.txt $W/truth_$MODE.txt $MK 2>&1 | grep -v Warning
log DONE
