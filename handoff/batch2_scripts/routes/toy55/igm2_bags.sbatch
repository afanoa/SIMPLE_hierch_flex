#!/bin/bash
#SBATCH -J igm2_bags
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 16:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# 2D sub-classes as pooled evidence in the per-region 3D boards of an IgM-like set, ONE job (truth only in readouts):
# (1) cls_expansion on ALL 2D classes (run_kernel_cls_all.sh, ~20 sub-classes per class); (2) sub-class + parent labels
# (stage1-build print) and their per-Fab readout; (3) pose metadata (near-true and refined poses), flags = junk-cleaned
# particles; (4) igm2_bagboard.sbatch array over the merged regions: near-true poses with bag beta 2 and 0, refined poses
# with beta 2; (5) factorised_readout.py per run, next to the plain Fortran-board results. env MODE (armsigm2), BETA (2), CONC (12)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; KL=$A/rnasep_klem; export MODE=${MODE:-armsigm2}
D2=$W/route2d_$MODE; E=$D2/expansion_all; O=$W/bags_$MODE; mkdir -p $O; log(){ echo "[$(date '+%F %T') igm2_bags $MODE] $*"; }
# ---- (1) full expansion
if [ ! -s $O/bags.txt ]; then
  [ -e $E/proj.simple ] || { log "cls_expansion on all 2D classes"; CLASSES=$E/classes.txt bash $W/run_kernel_cls_all.sh $E $D2/p/1_abinitio2D/p.simple 220 > $O/expansion_all.log 2>&1; grep -a -E "eligible|RC=|START|END|ERROR" $O/expansion_all.log | head -6; }
  # ---- (2) labels with the build that wrote the project
  ( B=$HOME/SIMPLE_stage1/build; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH; simple_exec prg=print_project_field projfile=$E/proj.simple oritype=ptcl2D mkdir=no > $O/ptcl2D_all.txt 2> $O/ptcl2D_all.err )
  awk '{c=0; p=0; for(i=1;i<=NF;i++){ if($i ~ /^class=/){split($i,a,"="); c=a[2]} if($i ~ /^cluster=/){split($i,b,"="); p=b[2]} }; print int(c) > "'$O'/bags.txt"; print int(p) > "'$O'/parents.txt"}' $O/ptcl2D_all.txt
fi
log "labels: $(grep -c . $O/bags.txt) rows, $(sort -n $O/bags.txt | uniq | wc -l) sub-classes, $(sort -n $O/parents.txt | uniq | wc -l) parents"
echo "== per-Fab readout of the full expansion"; $PY $W/subclass_angle_readout.py $O/bags.txt $W/sim$MODE/gt_particles.txt $W/truth_$MODE.txt 2>&1 | grep -v Warning
# ---- (3) pose metadata, flags = junk-cleaned particles with a sub-class
$PY $W/make_toy_meta.py $W/projtp_${MODE}_3/ptcl3D_tp.txt $O/meta_tp3_raw.npz > /dev/null
$PY - <<PYEOF
import numpy as np
clean = np.loadtxt('$W/factorised_$MODE/clean_subset.txt', dtype=int); bag = np.loadtxt('$O/bags.txt', dtype=int)
for tag, src in (('tp3', '$O/meta_tp3_raw.npz'), ('ref', '$KL/synth_$MODE/meta_toy.npz')):
    M = dict(np.load(src)); n = len(M['flags']); f = np.zeros(n, np.int32); f[clean - 1] = 1; f[:len(bag)][bag[:n] == 0] = 0; M['flags'] = f; np.savez('$O/meta_%s.npz' % tag, **M)
    print('meta %s: %d flagged particles (junk-cleaned, in a sub-class) of %d' % (tag, f.sum(), n))
PYEOF
# ---- (4) board array over the merged regions
CFG=$O/tasks.txt; : > $CFG
for spec in "tp3:$W/factorised3_${MODE}_tp3/merged4/regions:$W/projtp_${MODE}_3/cons_tp.mrc:${BETA:-2} 0" "ref:$W/factorised3_$MODE/merged4/regions:$W/proj$MODE/1_refine3D/recvol_state01.mrc:${BETA:-2}"; do
  IFS=: read tag rdir cons betas <<< "$spec"
  for b in $betas; do R=$O/run_${tag}_beta$b; mkdir -p $R/regions; for f in $(ls $rdir/region_*.mrc); do r=$(basename $f .mrc); r=${r#region_}; ln -sfn $f $R/regions/region_$r.mrc; mkdir -p $R/arm_$r
    echo "$R/arm_$r/board $f $O/meta_$tag.npz $cons $b $O/bags.txt $O/parents.txt" >> $CFG; done; done; done
N=$(grep -c . $CFG); log "$N boards (array, ${CONC:-12} at a time)"
J=$(sbatch --parsable --array=1-$N%${CONC:-12} --export=ALL,CFG=$CFG,MODE=$MODE $W/igm2_bagboard.sbatch); while squeue -h -j $J 2>/dev/null | grep -q .; do sleep 30; done
log "boards done: $(cat $W/igm2_bagboard.${J}_*.out 2>/dev/null | grep -c 'task done') of $N; failed: $(cat $W/igm2_bagboard.${J}_*.out 2>/dev/null | grep -c FAILED)"
# ---- (5) readouts
for R in $O/run_tp3_beta${BETA:-2} $O/run_tp3_beta0 $O/run_ref_beta${BETA:-2}; do echo "== $(basename $R)"; $PY $W/factorised_readout.py $R $W/sim$MODE/gt_particles.txt $W/truth_$MODE.txt $W/maps160/synth_$MODE/masks 2>&1 | grep -v Warning | grep -a -E "CV R|per arm"; done
echo "== Fortran boards (no bags), near-true poses"; grep -a "per arm, best region" $W/igm2_${MODE}_tp3_igm_factorised4.log | tail -1
echo "== Fortran boards (no bags), refined poses"; grep -a "per arm, best region" $W/igm2_${MODE}_igm_factorised4.log | tail -1
log DONE
