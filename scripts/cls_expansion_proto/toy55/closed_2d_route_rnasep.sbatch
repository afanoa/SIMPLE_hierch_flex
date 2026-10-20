#!/bin/bash
#SBATCH -J rnasep_closed_2d_route
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# The 2D route to the contact (user: the cls expansion shows the closed state): (1) every expansion sub-class average
# scored for closedness against pose-mean renders of deposited 18 (closed) vs 27 (open) at its members' poses
# (closed_subclasses.py); (2) members of the sub-classes above the 90th percentile of closedness (and > CMIN) form the
# 2D-closed set; (3) that set reconstructed from the DENOISED stacks (klem one-cell board, 6 A, family poses) and from
# the project's stacks (reconstruct3D at the family poses); (4) contact visibility at 8 A vs the deposited maps; montage.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel TOY_SMPD=1.464 OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; P=$A/rnasep_klem/expand1/fam1_pipeline; META=$A/rnasep_coupled/basis_span_all/meta_all.npz
CAVGS=$A/rnasep_coupled/full692k_2d/1_cls_expansion/cls_expansion_cavgs.mrc; [ -s $CAVGS ] || CAVGS=$A/rnasep_coupled/full692k_2d/1_cls_expansion/cls_expansion_cavgs_even.mrc
CONS=$A/rnasep_states/25_refine3D_states/rec_final_state01.mrc; LAB=$A/rnasep_softem/states.txt; STK=$(ls /mnt/beegfs/elmlund/rnasep/22_denoise_project/denoise_raw_particles_part??.mrcs | sort | tr "\n" , | sed 's/,$//')
D=$A/rnasep_klem/closed_2d_route; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') 2d] $*"; }
log "(1) closedness of every sub-class average ($CAVGS)"; $PY $W/closed_subclasses.py $META $CAVGS $W/maps160/docked $W/maps160/toy_consensus.mrc none $D/subcls --closed 18 --open 27 --lo 40 --hi 14 --nmem 300 2>&1 | grep -v "^  [0-9]* / " | tail -25
$PY - <<PYEOF
import numpy as np
R = np.loadtxt('$D/subcls_closedness.txt', comments='#'); sc = R[:, 4]; n = R[:, 2]; sub = R[:, 0].astype(int); q90 = np.percentile(sc, 90); cmin = ${CMIN:-0.2}; pick = sub[(sc > q90) & (sc > cmin)]
M = np.load('$META'); flags = M['flags'].astype(int); lab = M['lab692']; N = len(flags); ref = np.where(flags == 1)[0]; lab = lab if len(lab) == N else np.pad(np.zeros(0, int), 0) if False else lab; st = np.zeros(N, int); members = np.isin(lab, pick); st[(ref[members] if len(lab) == len(ref) else np.where(members)[0])] = 1
np.savetxt('$D/states.txt', st, fmt='%d'); np.savetxt('$D/closed2d_subset.txt', np.where(st == 1)[0] + 1, fmt='%d')
o = np.argsort(-sc); print('closedness: median %.3f, 90th pct %.3f, max %.3f; %d sub-classes picked (score > max(q90, %.2f)) with %d members; bottom decile holds %d members' % (np.median(sc), q90, sc.max(), len(pick), cmin, int(members.sum()), int(np.isin(lab, sub[sc < np.percentile(sc, 10)]).sum()))); print('top 10 sub-classes: ' + ' '.join('%d(n%d,%+.2f)' % (sub[k], n[k], sc[k]) for k in o[:10]))
PYEOF
log "(3a) denoised stacks: klem one-cell board at 6 A, family poses"; export OMP_NUM_THREADS=1
$PY $W/klem.py $STK none none $D/board6 --meta $P/meta_fb.npz --labels $LAB --cons $CONS --states 1,2,3,4 --axis=0,0,0,0 --n1 1 --n2 1 --nfree 0 --nproc 16 --init random --nsub 0 --subset $D/closed2d_subset.txt --ctf auto --lps "6 6" --sigmas "0.7 0.7" --mskdiam 232 --tag "closed2d denoised" 2>&1 | grep "subset\|Error\|Traceback"
export OMP_NUM_THREADS=8; log "(3b) project stacks: reconstruct3D at the family poses"
PS=$(ls -td $P/*_refine3D_states | head -1)/p.simple; cp -f $PS p.simple; simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$D/states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=1 nparts=6 nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-200)"; exit 1; }
RD=$(ls -td $D/*_reconstruct3D | head -1); mkdir -p $D/maps; for h in "" _even _odd; do cp -f $RD/recvol_state01$h.mrc $D/maps/ 2>/dev/null; done
$PY - <<PYEOF
import numpy as np, glob, sys
sys.path.insert(0, '$W'); from model_class import read_mrc, write_mrc, Shells, fftc
from scipy.ndimage import gaussian_filter
W = '$W'; SMPD = 1.464; n = 160; S = Shells(n, SMPD); sh = S.sh
v = read_mrc('$D/board6/cells/cell_00_00.mrc'); nc = v.shape[0]; F = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v))); Fp = np.zeros((n, n, n), complex); c = n // 2; h = nc // 2; Fp[c - h:c - h + nc, c - h:c - h + nc, c - h:c - h + nc] = F
up = np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(Fp)))) * (n / nc) ** 3; write_mrc('$D/closed2d_denoised_6A.mrc', up.astype(np.float32), SMPD)
zz, yy, xx = np.indices((n, n, n)); rr = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2); rA = rr * SMPD
lp = lambda v, res: np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v))) * (rr <= n * SMPD / res)))))
cons = read_mrc(W + '/maps160/toy_consensus.mrc'); c18 = read_mrc(glob.glob(W + '/maps160/docked/conf18_*.mrc')[0]); c27 = read_mrc(glob.glob(W + '/maps160/docked/conf27_*.mrc')[0])
g = lambda v: gaussian_filter(v, 2.5); L18, Lc = g(c18), g(cons); dens18 = L18 > L18.mean() + 2.0 * L18.std(); diff = L18 - Lc; arm = (diff > diff.mean() + 3 * diff.std()) & dens18; core = dens18 & (np.abs(diff) < 0.5 * diff.std()); tip = arm & (rA > 95)
def fsc_res(x, y):
    a_, b_ = fftc(x), fftc(y); num = np.bincount(sh.ravel(), np.real(np.conj(a_) * b_).ravel()); d1 = np.bincount(sh.ravel(), np.abs(a_.ravel()) ** 2); d2 = np.bincount(sh.ravel(), np.abs(b_.ravel()) ** 2)
    f = num / np.maximum(np.sqrt(d1 * d2), 1e-30); s = 1
    while s < 80 and f[s] > 0.143: s += 1
    return n * SMPD / s
raw = read_mrc('$D/maps/recvol_state01.mrc'); st = np.loadtxt('$D/states.txt', dtype=int)
print('2D-closed set: %d particles; raw reconstruction FSC %.1f A' % (int((st == 1).sum()), fsc_res(read_mrc('$D/maps/recvol_state01_even.mrc'), read_mrc('$D/maps/recvol_state01_odd.mrc'))))
print('contact VISIBILITY at 8 A (core median = 1):  %-34s %8s %8s %9s' % ('map', 'tip p99', 'tip p90', 'tip>0.5'))
for name, vv in (('deposited 18 (closed)', c18), ('deposited 27 (open)', c27), ('2D-closed set, raw reconstruct3D', raw), ('2D-closed set, denoised (klem 6A)', up), ('docked 15k (focal), raw', read_mrc('$P/klem_fb/closed_focal_closedfam/rank/maps/recvol_state01.mrc'))):
    L = lp(vv, 8.); L = L / max(np.median(L[core]), 1e-30); print('  %-50s %8.2f %8.2f %9.2f' % (name, np.percentile(L[tip], 99), np.percentile(L[tip], 90), np.mean(L[tip] > 0.5)))
PYEOF
$PY $W/montage_maps.py $W/figs_now/closed2d_lp8.png --lp 8 "$D/closed2d_denoised_6A.mrc:2D-closed set, denoised" "$D/maps/recvol_state01.mrc:2D-closed set, raw" "$(ls $W/maps160/docked/conf18_*.mrc | head -1):deposited 18 (closed)" "$(ls $W/maps160/docked/conf27_*.mrc | head -1):deposited 27 (open)" 2>&1 | grep -v Warning
log DONE
