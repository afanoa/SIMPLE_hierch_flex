#!/bin/bash
#SBATCH -J rnasep_klem_pose_policy
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Node alignment policy on the real docked node (top N1 of the family-1 focal ranking): klem one-cell boards at the
# family poses, 10 -> 8 A, with klem's strictly local pose step (3x3x3 rotation grid x 2 refinements, Gaussian prior,
# no shifts, no re-assignment) against (a) nothing (baseline, no pose step), (b) the cell reference, everything
# (lever arm), (c) the cell reference masked to the CORE (density x (1 - reproducible variability)), (d) masked to the
# REGION (file masks: the P12/P10.1 focal mask and the density minus it; a one-cell board has no cell variability). Readout per run: closed-arm / tip / base occupancy of the one-cell map (closed_arm_cells.py,
# same shrinkage for all so the comparison is fair) and the pose change vs the family poses.
# env: FD (board dir), N1 (30000), STEP (3), PRIOR (4), EVERY (2), LPS ("10 10 10 8 8 8"), MSKD (232)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; E=$A/rnasep_klem/expand1; P=$E/fam1_pipeline; KD=$P/klem_fb; META=$P/meta_fb.npz
CONS=$A/rnasep_states/25_refine3D_states/rec_final_state01.mrc; LAB=$A/rnasep_softem/states.txt; STK=$(ls /mnt/beegfs/elmlund/rnasep/22_denoise_project/denoise_raw_particles_part??.mrcs | sort | tr "\n" , | sed 's/,$//')
FD=${FD:-$KD/closed_focal}; SUBSET=${SUBSET:-$FD/lever/subset_top30k.txt}; D=$FD/pose_policy; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') pp] $*"; }; export OMP_NUM_THREADS=1
run(){ tag=$1; shift; log "run $tag: $*"; $PY $W/klem.py $STK none none $D/$tag --meta $META --labels $LAB --cons $CONS --states 1,2,3,4 --axis=0,0,0,0 --n1 1 --n2 1 --nfree 0 --nproc 16 --init random --nsub 0 --subset $SUBSET --ctf auto \
    --lps "${LPS:-10 10 10 8 8 8}" --sigmas "0.7 0.7 0.7 0.7 0.7 0.7" --mskdiam ${MSKD:-232} --tag "pp $tag" "$@" 2>&1 | grep "subset\|pose step\|Error\|Traceback" | cut -c1-220
  $PY $W/closed_arm_cells.py $D/$tag --top 1 2>&1 | grep -v Warning | sed "s/^/[$tag] /"; }
$PY - <<PYEOF2
import numpy as np, sys
sys.path.insert(0, '$W'); from model_class import read_mrc, write_mrc
from scipy.ndimage import gaussian_filter, binary_dilation
reg = read_mrc('$FD/focal_mask_p12.mrc'); cons = read_mrc('$CONS'); g = gaussian_filter(cons, 2.0); dens = g > g.mean() + 1.5 * g.std()
core = dens & ~binary_dilation(reg > 0.5, iterations=3); soft = gaussian_filter(core.astype(float), 2.0); soft = np.clip(soft / max(soft.max(), 1e-30), 0, 1)
write_mrc('$D/core_mask.mrc', soft.astype(np.float32), 1.464); print('masks: core %d voxels (density %d minus the dilated P12/P10.1 region), region %d' % ((soft > 0.5).sum(), dens.sum(), (reg > 0.5).sum()))
PYEOF2
run baseline
run lever    --pose_every ${EVERY:-2} --pose_step ${STEP:-3} --pose_prior ${PRIOR:-4} --pose_ref cell
run core     --pose_every ${EVERY:-2} --pose_step ${STEP:-3} --pose_prior ${PRIOR:-4} --pose_ref cell --pose_mask $D/core_mask.mrc
run region   --pose_every ${EVERY:-2} --pose_step ${STEP:-3} --pose_prior ${PRIOR:-4} --pose_ref cell --pose_mask $FD/focal_mask_p12.mrc
$PY - <<PYEOF
import numpy as np, glob
M = np.load('$META'); e = np.stack([M['e1'], M['e2'], M['e3']], 1)
def rotmat(ang, axis):
    r = np.radians(ang); c, s = np.cos(r), np.sin(r); R = np.eye(3)
    if axis == 3: R[0, 0], R[0, 1], R[1, 0], R[1, 1] = c, -s, s, c
    else:         R[0, 0], R[0, 2], R[2, 0], R[2, 2] = c, -s, s, c
    return R
def e2m(v): return rotmat(v[2], 3) @ rotmat(v[1], 2) @ rotmat(v[0], 3)
for tag in ('lever', 'core', 'region'):
    fs = sorted(glob.glob('$D/%s/ptcl3D_pose_it*.txt' % tag), key=lambda f: int(f.split('_it')[-1].split('.')[0]))
    if not fs: print('%s: no pose dump' % tag); continue
    rows = [dict(t.split('=', 1) for t in l.split()) for l in open(fs[-1]) if l.strip()]; pind = np.array([int(r['pind']) for r in rows]) - 1; new = np.array([[float(r['e1']), float(r['e2']), float(r['e3'])] for r in rows])
    ang = np.array([np.degrees(np.arccos(np.clip((np.trace(e2m(new[i]) @ e2m(e[pind[i]]).T) - 1) / 2, -1, 1))) for i in range(len(pind))])
    print('%s: pose change after %s: angle median %.2f deg (p90 %.2f, max %.1f) over %d particles' % (tag, fs[-1].split('/')[-1], np.median(ang), np.percentile(ang, 90), ang.max(), len(ang)))
PYEOF
log DONE
