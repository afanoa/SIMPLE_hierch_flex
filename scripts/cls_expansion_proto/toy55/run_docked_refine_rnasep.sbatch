#!/bin/bash
#SBATCH -J rnasep_docked_refine
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Per-state refinement of the DOCKED closed set in the paper's order: (1) a rigid-core mask from where the ten deposited
# conformers agree (low variance across conformers inside the density, dilated, soft); (2) refine3D of the docked
# particles with that core mask (local poses, 10 -> LPSTOP A), started from the family poses; (3) plain unmasked
# reconstruction at the core-aligned poses, scored on arm / tip / base; (4) a focal refinement of the P12/P10.1
# periphery with the focal mask from the core-aligned poses, reconstructed and scored the same way.
# env: STATES (one int per project row, 1 = docked), CORE_MSKDIAM (170), LPSTOP (8), MAXITS (6), NP (10), RTAG
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8 TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; P=$A/rnasep_klem/expand1/fam1_pipeline; KF=$P/klem_fb
STATES=${STATES:-$KF/closed_focal_closedfam/rank/states.txt}; FOCAL=${FOCAL:-$KF/closed_focal_closedfam/focal_mask_p12.mrc}; D=$KF/docked_refine${RTAG:-}; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') dr] $*"; }
SC=$KF/closed_focal_closedfam/rank/score.py
$PY - <<PYEOF
import numpy as np, glob, sys
sys.path.insert(0, '$W'); from model_class import read_mrc, write_mrc
from scipy.ndimage import gaussian_filter, binary_dilation, binary_erosion
W = '$W'; SMPD = 1.464; ids = [1, 10, 18, 27, 28, 37, 46, 48, 49, 51]; g = lambda v: gaussian_filter(v, 2.0)
C = np.array([g(read_mrc(glob.glob(W + '/maps160/docked/conf%02d_*.mrc' % i)[0])) for i in ids]); mean = C.mean(0); sd = C.std(0); dens = mean > mean.mean() + 1.5 * mean.std()
rel = sd / np.maximum(mean, 1e-6); thr = np.percentile(rel[dens], 50); core = dens & (rel < thr); core = binary_dilation(binary_erosion(core, iterations=1), iterations=3)
soft = gaussian_filter(core.astype(float), 2.0); soft = np.clip(soft / max(soft.max(), 1e-30), 0, 1); write_mrc('$D/core_mask.mrc', soft.astype(np.float32), SMPD)
n = 160; zz, yy, xx = np.indices((n, n, n)); r = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2) * SMPD
st = np.loadtxt('$STATES', dtype=int); sel = (st == 1).astype(int); np.savetxt('$D/states.txt', sel, fmt='%d')
print('core mask: %d dense voxels, %d in the core (relative sd below the median %.3f), radius of the core 90th pct %.0f A; docked particles %d' % (dens.sum(), core.sum(), thr, np.percentile(r[core], 90), sel.sum()))
PYEOF
PS=$(ls -td $P/*_refine3D_states | head -1)/p.simple; cp -f $PS p.simple; simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$D/states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
log "(2) refine3D on the core: spherical mask of diameter ${CORE_MSKDIAM:-170} A (the arm tip at 95-115 A radius excluded), local poses, 10 -> ${LPSTOP:-8} A"
simple_exec prg=refine3D projfile=p.simple pose_policy=local lpstart=10 lpstop=${LPSTOP:-8} maxits=${MAXITS:-6} mskdiam=${CORE_MSKDIAM:-170} pgrp=c1 nparts=${NP:-10} nthr=8 > refine_core.log 2>&1 || { log "core refine3D FAILED: $(grep -a -m2 ERROR refine_core.log | cut -c1-200)"; tail -3 refine_core.log; exit 1; }
RD=$(ls -td $D/*_refine3D | head -1); cp -f $RD/p.simple p_core.simple
log "(3) plain reconstruction at the core-aligned poses"; mkdir -p rec_core && cd rec_core && cp -f ../p_core.simple p.simple
simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=1 nparts=6 nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED"; exit 1; }
RR=$(ls -td ./*_reconstruct3D | head -1); mkdir -p maps; for h in "" _even _odd; do cp -f $RR/recvol_state01$h.mrc maps/ 2>/dev/null; done; cd $D
$PY $SC $KF/closed_focal_closedfam/rank/maps/recvol_state01.mrc "docked, family poses (before)" $D/rec_core/maps/recvol_state01.mrc "docked, core-aligned poses"
log "(4) skipped: this SIMPLE has no mask-file refine3D (only pcg_mskfile); periphery stays at the core-aligned poses"
log DONE
