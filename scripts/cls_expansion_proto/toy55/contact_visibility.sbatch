#!/bin/bash
#SBATCH -J rnasep_contact_vis
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3
$PY - <<PYEOF
import numpy as np, glob, sys
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/routes/toy55'); from model_class import read_mrc
from scipy.ndimage import gaussian_filter
W = '/mnt/beegfs/elmlund/afan/routes/toy55'; SMPD = 1.464; n = 160; zz, yy, xx = np.indices((n, n, n)); rr = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2); rA = rr * SMPD
def lp(v, res): return np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v))) * (rr <= n * SMPD / res)))))
cons = read_mrc(W + '/maps160/toy_consensus.mrc'); c18 = read_mrc(glob.glob(W + '/maps160/docked/conf18_*.mrc')[0]); c27 = read_mrc(glob.glob(W + '/maps160/docked/conf27_*.mrc')[0])
g = lambda v: gaussian_filter(v, 2.5); L18, Lc = g(c18), g(cons); dens18 = L18 > L18.mean() + 2.0 * L18.std(); diff = L18 - Lc; arm = (diff > diff.mean() + 3 * diff.std()) & dens18; core = dens18 & (np.abs(diff) < 0.5 * diff.std()); tip = arm & (rA > 95)
print('contact VISIBILITY at 8 A low-pass: every map scaled so that its core (median of the core voxels) = 1; threshold = half the core level (a typical display contour)')
print('  %-36s %8s %8s %10s %10s' % ('map', 'tip p99', 'tip p90', 'tip>0.5', 'tip>0.35'))
for name, path in (('deposited 18 (closed)', None), ('deposited 27 (open)', 'c27'), ('consensus', 'cons'), ('all 18-matched 40k', '/mnt/beegfs/elmlund/afan/rnasep_klem/expand1/fam1_pipeline/klem_fb/matched_rec/maps/recvol_state03.mrc'), ('docked top 15k', '/mnt/beegfs/elmlund/afan/rnasep_klem/expand1/fam1_pipeline/klem_fb/closed_focal_closedfam/rank/maps/recvol_state01.mrc'), ('undocked bottom 15k', '/mnt/beegfs/elmlund/afan/rnasep_klem/expand1/fam1_pipeline/klem_fb/closed_focal_closedfam/rank/maps/recvol_state03.mrc')):
    v = c18 if path is None else (c27 if path == 'c27' else (cons if path == 'cons' else read_mrc(path))); L = lp(v, 8.); L = L / max(np.median(L[core]), 1e-30)
    print('  %-36s %8.2f %8.2f %10.2f %10.2f' % (name, np.percentile(L[tip], 99), np.percentile(L[tip], 90), np.mean(L[tip] > 0.5), np.mean(L[tip] > 0.35)))
PYEOF
