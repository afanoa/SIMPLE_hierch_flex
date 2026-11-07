#!/bin/bash
#SBATCH -J rnasep_closed_refined_rec
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# The closed-state local refinement's POSES re-measured with a plain reconstruct3D (its own output map is masked and
# filtered, so its 0.60 is not comparable): arm / tip / base occupancy vs the plain reconstruction at the family poses.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8 TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; MR=$A/rnasep_klem/expand1/fam1_pipeline/klem_fb/matched_rec; CR=$MR/closed_refine
D=$CR/rec; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') crr] $*"; }
RD=$(ls -td $CR/*_refine3D | head -1); cp -f $RD/p.simple p.simple; cp -f $CR/states.txt $D/states.txt
simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$D/states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=1 nparts=6 nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-200)"; exit 1; }
RR=$(ls -td $D/*_reconstruct3D | head -1); mkdir -p $D/maps; for h in "" _even _odd; do cp -f $RR/recvol_state01$h.mrc $D/maps/ 2>/dev/null; done
$PY - <<PYEOF
import numpy as np, glob, sys
sys.path.insert(0, '$W'); from model_class import read_mrc, Shells, fftc
from scipy.ndimage import gaussian_filter
W = '$W'; SMPD = 1.464; n = 160; S = Shells(n, SMPD); sh = S.sh; g = lambda v: gaussian_filter(v, 2.5)
cons = read_mrc(W + '/maps160/toy_consensus.mrc'); c18 = read_mrc(glob.glob(W + '/maps160/docked/conf18_*.mrc')[0]); L18, Lc = g(c18), g(cons); dens18 = L18 > L18.mean() + 2.0 * L18.std(); diff = L18 - Lc
arm = (diff > diff.mean() + 3 * diff.std()) & dens18; core = dens18 & (np.abs(diff) < 0.5 * diff.std()); ref = L18[arm].mean() / L18[core].mean()
zz, yy, xx = np.indices((n, n, n)); r = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2) * SMPD; tip = arm & (r > 95); base = arm & (r < 80)
def fsc_res(x, y):
    a_, b_ = fftc(x), fftc(y); num = np.bincount(sh.ravel(), np.real(np.conj(a_) * b_).ravel()); d1 = np.bincount(sh.ravel(), np.abs(a_.ravel()) ** 2); d2 = np.bincount(sh.ravel(), np.abs(b_.ravel()) ** 2)
    f = num / np.maximum(np.sqrt(d1 * d2), 1e-30); s = 1
    while s < 80 and f[s] > 0.143: s += 1
    return n * SMPD / s
for name, path in (('family poses (plain, before)', '$MR/maps/recvol_state03.mrc'), ('closed-refined poses (plain)', '$D/maps/recvol_state01.mrc')):
    v = read_mrc(path); Lv = g(v); c = Lv[core].mean(); occ = (Lv[arm].mean() / c) / ref; ot = (Lv[tip].mean() / c) / (L18[tip].mean() / L18[core].mean()); ob = (Lv[base].mean() / c) / (L18[base].mean() / L18[core].mean())
    try: res = fsc_res(read_mrc(path.replace('.mrc', '_even.mrc')), read_mrc(path.replace('.mrc', '_odd.mrc')))
    except Exception: res = float('nan')
    print('  %-32s: arm %.2f | tip %.2f | base %.2f | FSC %.1f A' % (name, occ, ot, ob, res))
PYEOF
log DONE
