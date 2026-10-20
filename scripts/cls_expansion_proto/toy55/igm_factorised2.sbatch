#!/bin/bash
#SBATCH -J igm_factorised2
#SBATCH -p csbdevel
#SBATCH -c 48
#SBATCH --mem=160G
#SBATCH -t 08:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Truth-free factorised route, region step v2 (truth only in readouts): the junk-cleaned set of igm_factorised.sbatch;
# (1) periphery = consensus envelope minus the data-derived rigid core (periphery_mask.py on the full-set refinement map);
# (2) a finer-band whole-particle board (14 -> 8 A, 7x7 + 1 free, halves) started from the 20 -> 10 A board's coordinates;
# (3) motion-pattern regions restricted to the periphery (derive_regions2 --within); (4) one balanced focused board per
# region (as before); (5) factorised_readout.py. env MODE (armsigm1), NLAT (5), KREG (8), MAXREG (8), BPAR (3), CORE (0.40)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; KL=$A/rnasep_klem; MODE=${MODE:-armsigm1}
O1=$W/factorised_$MODE; O=$W/factorised2_$MODE; mkdir -p $O; cd $O; export FB=$KL/f692e/klem2_bin SIMPLE_QSYS=local; PS=$KL/synth_$MODE/src.simple; NL=${NLAT:-5}
log(){ echo "[$(date '+%F %T') factorised2 $MODE] $*"; }
[ -s $O1/clean_subset.txt ] && [ -s $O1/parent/klem_latent.txt ] || { log "needs $O1/clean_subset.txt and $O1/parent (run igm_factorised.sbatch first)"; exit 1; }
# ---- (1) periphery mask from the consensus
$PY $W/periphery_mask.py $W/proj$MODE/1_refine3D/recvol_state01.mrc $O/periphery.mrc --core ${CORE:-0.40} --check $W/maps160/synth_$MODE/masks 2>&1 | grep -v Warning
# ---- (2) finer-band board from the coarse board's coordinates
awk '!/^#/ {printf "%d %.4f %.4f\n", $1, $2, $3}' $O1/parent/klem_latent.txt > $O/init_fine.txt
log "finer-band board (14 -> 8 A) on $(grep -c . $O1/clean_subset.txt) cleaned particles"
[ -s $O/fine/klem_latent.txt ] || SERVER_LIST= NPARTS= NTHR=48 N1=7 N2=7 NFREE=1 INIT=file INITFILE=$O/init_fine.txt WRITE_HALVES=1 bash $W/fboard.sh $O/fine $PS $O1/clean_subset.txt "14 14 12 12 10 10 8 8 8" "1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7" > $O/fine.out 2>&1 || { log "fine board FAILED: $(tail -2 $O/fine.out)"; exit 1; }
# ---- (3) regions inside the periphery
[ -s $O/regions/region_1.mrc ] || $PY $W/derive_regions2.py $O/fine $O/regions --k ${KREG:-8} --maxreg ${MAXREG:-8} --within $O/periphery.mrc 2>&1 | grep -v Warning | tail -12
NR=$(ls $O/regions/region_*.mrc 2>/dev/null | wc -l); log "$NR data-derived regions inside the periphery"; [ $NR -gt 0 ] || exit 1
# ---- (4) one focused board per region
for f in $(ls $O/regions/region_*.mrc); do r=$(basename $f .mrc); r=${r#region_}; D=$O/arm_$r; mkdir -p $D
  while [ $(jobs -rp | wc -l) -ge ${BPAR:-3} ]; do sleep 10; done
  [ -s $D/board/klem_latent.txt ] || ( SERVER_LIST= NPARTS= NTHR=16 N1=$NL N2=$NL NFREE=1 CTFABS=1 PSCALE=1 FOCUSMASK=$f INIT=random SEED=$r bash $W/fboard.sh $D/board $PS $O1/clean_subset.txt "20 14 12 10 10 8 8 8" "1.5 1.2 1.0 1.0 0.8 0.8 0.7 0.7" > $D/board.out 2>&1 || echo "region $r board FAILED" ) &
done; wait; log "region boards done"
# ---- (5) readout (truth)
$PY $W/factorised_readout.py $O $W/sim$MODE/gt_particles.txt $W/truth_$MODE.txt $W/maps160/synth_$MODE/masks 2>&1 | grep -v Warning
log DONE
