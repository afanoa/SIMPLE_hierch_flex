#!/bin/bash
#SBATCH -J rnasep_closed_refine
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Can identify + refine give the CLOSED conformer? Three states from the full-run latent (closed-enriched end, middle,
# open end: equal-count bins along the main axis) -> refine3D_states (local 20->8 A) on all 692k -> the refined state
# maps scored by closed-arm occupancy (vs deposited 18 = 1.0, consensus 0.5) and by correlation with the deposited
# conformers. env: KDIR (klem_full_npz), NB (3), MAXITS (6), NP (10)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
BK=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$BK PATH=$BK/bin:$BK/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8 TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; KDIR=${KDIR:-$A/rnasep_klem/klem_full_npz}; NB=${NB:-3}; NP=${NP:-10}
SRC=$A/rnasep_coupled/full692k_2d/rnasep.simple; META=$A/rnasep_coupled/basis_span_all/meta_all.npz; D=$A/rnasep_klem/closed_refine_nb$NB; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') cr] $*"; }
$PY $W/klem_bins.py $KDIR $META $D/states.txt --nb $NB
cp -f $SRC p.simple; simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$D/states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
log "refine3D_states nstates=$NB local 20->8 A maxits ${MAXITS:-6} on the bins of $KDIR"
simple_exec prg=refine3D_states projfile=p.simple nstates=$NB pose_policy=local lpstart=20 lpstop=8 maxits=${MAXITS:-6} mskdiam=220 pgrp=c1 nparts=$NP nthr=8 > refine.log 2>&1 || { log "refine3D_states FAILED: $(grep -a -m2 ERROR refine.log | cut -c1-200)"; exit 1; }
RD=$(ls -td $D/*_refine3D_states | head -1); mkdir -p $D/maps; for w in $(seq -f %02g 1 $NB); do for h in "" _even _odd; do cp -f $RD/recvol_state$w$h.mrc $D/maps/ 2>/dev/null; done; done
simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null | awk '{for(i=1;i<=NF;i++) if($i ~ /^state=/){split($i,a,"="); print int(a[2])}}' | sort -n | uniq -c | awk '{printf "state %d: %d particles\n", $2, $1}'
$PY - <<PYEOF
import numpy as np, glob
from scipy.ndimage import gaussian_filter
SMPD = 1.464; W = '$W'
def read_mrc(p):
    h = np.fromfile(p, dtype=np.int32, count=256); nx, ny, nz = h[0], h[1], h[2]; return np.fromfile(p, dtype=np.float32, offset=1024 + h[23]).reshape(nz, ny, nx).astype(float)
c18 = read_mrc(glob.glob(W + '/maps160/docked/conf18_*.mrc')[0]); cons = read_mrc(W + '/maps160/toy_consensus.mrc'); lp = lambda v: gaussian_filter(v, 2.5); L18, Lc = lp(c18), lp(cons)
dens18 = L18 > L18.mean() + 2.0 * L18.std(); diff = L18 - Lc; arm = (diff > diff.mean() + 3 * diff.std()) & dens18; core = dens18 & (np.abs(diff) < 0.5 * diff.std()); ref = L18[arm].mean() / L18[core].mean()
n = 160; zz, yy, xx = np.indices((n, n, n)); r = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2); msk = 1. / (1. + np.exp((r - 100 / SMPD) / 2.)); F = lambda v: np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v * msk))); band = (r >= n * SMPD / 60.) & (r <= n * SMPD / 20.)
ids = [1, 10, 18, 27, 28, 37, 46, 48, 49, 51]; T = {i: F(read_mrc(glob.glob(W + '/maps160/docked/conf%02d_*.mrc' % i)[0]) - cons)[band] for i in ids}
M = [read_mrc('$D/maps/recvol_state%02d.mrc' % w) for w in range(1, $NB + 1)]; Mm = np.mean(M, 0)
def bc(x, y): return float(np.real(np.sum(np.conj(x) * y)) / np.sqrt(np.sum(np.abs(x) ** 2) * np.sum(np.abs(y) ** 2)))
for w, V in enumerate(M):
    Lv = lp(V); occ = (Lv[arm].mean() / Lv[core].mean()) / ref; Dw = F(V - Mm)[band]; cc = {i: bc(Dw, T[i]) for i in ids}; best = max(cc, key=cc.get)
    print('refined state %d: closed-arm occupancy %.2f (bin before refinement: see latent bins; consensus 0.50, deposited 18 = 1.00) | best conformer %d (%+.2f) | all: %s' % (w + 1, occ, best, cc[best], ' '.join('%d:%+.2f' % (i, cc[i]) for i in ids)))
PYEOF
log DONE
