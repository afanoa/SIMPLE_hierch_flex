#!/bin/bash
#SBATCH -J klem_core_frame
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Is the core-aligned frame the docked frame? Reconstruct family 1 at the core-aligned poses (all particles) and dock it
# onto the toy consensus (the frame the 55 deposited maps were docked into): the transform tells whether the strict
# recovery drop (7 -> 3) is a frame offset
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/f692e; O=$R/core_frame; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') coreframe] $*"; }
cp -f $R/core_align/aligned_f1.simple p.simple; simple_exec prg=update_project projfile=p.simple qsys_name=local nparts=6 nthr=16 mkdir=no > up.log 2>&1
log "reconstruct3D of family 1 at the core-aligned poses"; simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=1 nparts=6 nthr=16 > rec.log 2>&1 || { log "reconstruct3D FAILED"; exit 1; }
RR=$(ls -td $O/*_reconstruct3D | head -1); cp -f $RR/recvol_state01.mrc $O/fam1_corealigned.mrc; find "${RR:?}" -name "*_part*.mrc" -delete
log "dock the core-aligned family map onto the toy consensus (docked frame)"; simple_exec prg=dock_volpair vol1=$W/maps160/toy_consensus.mrc vol2=$O/fam1_corealigned.mrc smpd=1.464 lpstart=20 lpstop=8 mskdiam=232 nthr=16 mkdir=no > dock.log 2>&1 || { log "dock_volpair FAILED: $(grep -a -m2 ERROR dock.log | cut -c1-160)"; }
grep -a -i "rotation\|shift\|corr\|euler\|e1\|xform" dock.log | tail -6 | cut -c1-200
$PY - <<PYEOF
import numpy as np, sys
sys.path.insert(0, '$W'); from model_class import read_mrc
a = read_mrc('$W/maps160/toy_consensus.mrc').astype(np.float32); b = read_mrc('$O/fam1_corealigned.mrc').astype(np.float32); c = read_mrc('$R/parents/fam1.mrc').astype(np.float32)
n = a.shape[0]; f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); k2 = kx**2 + ky**2 + kz**2; lp = lambda v: np.real(np.fft.ifftn(np.fft.fftn(v) * np.exp(-3. * k2 / (1.464 / 10.) ** 2)))
def cc(u, v): u, v = lp(u), lp(v); u -= u.mean(); v -= v.mean(); return float((u * v).sum() / np.sqrt((u * u).sum() * (v * v).sum()))
print('10 A correlation, no transform: toy consensus vs family-1 parent (family poses) %.3f; vs family-1 core-aligned %.3f; parent vs core-aligned %.3f' % (cc(a, c), cc(a, b), cc(c, b)))
PYEOF
log DONE
