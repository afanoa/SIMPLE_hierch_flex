#!/bin/bash
#SBATCH -J rnasep_docked_denoised
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Does the contact come out of the DENOISED particles (the stacks the E-step and the 2D expansion read) where the raw
# reconstruct3D map stays faint? The 15k docked particles reconstructed by klem's own cross-half backprojection (one-cell
# board, 6 A band, wide mask) from the denoised stacks at the family poses; the coarse cube Fourier-padded to the 160 box;
# contact visibility at 8 A against the raw reconstruction and the deposited closed/open maps; montage. Also lists the
# stacks the two projects reference.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH TOY_SMPD=1.464 OMP_NUM_THREADS=1
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; P=$A/rnasep_klem/expand1/fam1_pipeline; KF=$P/klem_fb; RK=$KF/closed_focal_closedfam/rank
CONS=$A/rnasep_states/25_refine3D_states/rec_final_state01.mrc; LAB=$A/rnasep_softem/states.txt; STK=$(ls /mnt/beegfs/elmlund/rnasep/22_denoise_project/denoise_raw_particles_part??.mrcs | sort | tr "\n" , | sed 's/,$//')
D=$RK/denoised; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') dd] $*"; }
echo "== stacks referenced by the projects (first rows of the stk segment)"
simple_exec prg=print_project_field projfile=$A/rnasep_coupled/full692k_2d/rnasep.simple oritype=stk mkdir=no 2>/dev/null | head -2 | grep -o 'stk=[^ ]*' | sed 's/^/  2D expansion project: /'
simple_exec prg=print_project_field projfile=$(ls -td $P/*_refine3D_states | head -1)/p.simple oritype=stk mkdir=no 2>/dev/null | head -2 | grep -o 'stk=[^ ]*' | sed 's/^/  refinement project:   /'
$PY -c "
import numpy as np; st = np.loadtxt('$RK/states.txt', dtype=int); p = np.where(st == 1)[0] + 1; np.savetxt('$D/docked_subset.txt', p, fmt='%d'); print('docked particles:', len(p))"
log "one-cell board of the docked particles from the DENOISED stacks, 6 A band"
$PY $W/klem.py $STK none none $D/board6 --meta $P/meta_fb.npz --labels $LAB --cons $CONS --states 1,2,3,4 --axis=0,0,0,0 --n1 1 --n2 1 --nfree 0 --nproc 16 --init random --nsub 0 --subset $D/docked_subset.txt --ctf auto --lps "6 6" --sigmas "0.7 0.7" --mskdiam 232 --tag "docked denoised 6A" 2>&1 | grep "subset\|Error\|Traceback"
$PY - <<PYEOF
import numpy as np, glob, sys
sys.path.insert(0, '$W'); from model_class import read_mrc, write_mrc
from scipy.ndimage import gaussian_filter
W = '$W'; SMPD = 1.464; n = 160
v = read_mrc('$D/board6/cells/cell_00_00.mrc'); nc = v.shape[0]; F = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v))); Fp = np.zeros((n, n, n), complex); c = n // 2; h = nc // 2; Fp[c - h:c - h + nc, c - h:c - h + nc, c - h:c - h + nc] = F
up = np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(Fp)))) * (n / nc) ** 3; write_mrc('$D/docked_denoised_6A.mrc', up.astype(np.float32), SMPD); print('coarse cube %d -> 160 box written' % nc)
zz, yy, xx = np.indices((n, n, n)); rr = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2); rA = rr * SMPD
lp = lambda v, res: np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v))) * (rr <= n * SMPD / res)))))
cons = read_mrc(W + '/maps160/toy_consensus.mrc'); c18 = read_mrc(glob.glob(W + '/maps160/docked/conf18_*.mrc')[0]); c27 = read_mrc(glob.glob(W + '/maps160/docked/conf27_*.mrc')[0])
g = lambda v: gaussian_filter(v, 2.5); L18, Lc = g(c18), g(cons); dens18 = L18 > L18.mean() + 2.0 * L18.std(); diff = L18 - Lc; arm = (diff > diff.mean() + 3 * diff.std()) & dens18; core = dens18 & (np.abs(diff) < 0.5 * diff.std()); tip = arm & (rA > 95)
print('contact VISIBILITY at 8 A (core median = 1):  %-30s %8s %8s %9s' % ('map', 'tip p99', 'tip p90', 'tip>0.5'))
for name, vv in (('deposited 18 (closed)', c18), ('deposited 27 (open)', c27), ('docked 15k, RAW reconstruct3D', read_mrc('$RK/maps/recvol_state01.mrc')), ('docked 15k, DENOISED (klem 6A)', up)):
    L = lp(vv, 8.); L = L / max(np.median(L[core]), 1e-30); print('  %-50s %8.2f %8.2f %9.2f' % (name, np.percentile(L[tip], 99), np.percentile(L[tip], 90), np.mean(L[tip] > 0.5)))
PYEOF
$PY $W/montage_maps.py $W/figs_now/docked_denoised_vs_raw_lp8.png --lp 8 "$D/docked_denoised_6A.mrc:docked 15k, denoised particles" "$RK/maps/recvol_state01.mrc:docked 15k, raw particles" "$(ls $W/maps160/docked/conf18_*.mrc | head -1):deposited 18 (closed)" "$(ls $W/maps160/docked/conf27_*.mrc | head -1):deposited 27 (open)" 2>&1 | grep -v Warning
log DONE
