#!/bin/bash
#SBATCH -J rnasep_sharpen_docked
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; RK=/mnt/beegfs/elmlund/afan/rnasep_klem/expand1/fam1_pipeline/klem_fb/closed_focal_closedfam/rank
$PY - <<PYEOF
import numpy as np, sys
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/routes/toy55'); from model_class import read_mrc, write_mrc
SMPD = 1.464; n = 160; zz, yy, xx = np.indices((n, n, n)); s_inv = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2) / (n * SMPD)
for w, name, res in ((1, 'docked_top15k', 6.5), (3, 'undocked_bottom15k', 7.6)):
    v = read_mrc('/mnt/beegfs/elmlund/afan/rnasep_klem/expand1/fam1_pipeline/klem_fb/closed_focal_closedfam/rank/maps/recvol_state%02d.mrc' % w); F = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v))); F *= np.exp(150. * s_inv ** 2 / 4.)
    fc = 1. / res; F *= 0.5 * (1 + np.cos(np.pi * np.clip((s_inv - 0.85 * fc) / (0.3 * fc), 0, 1))); vs = np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(F))))
    write_mrc('/mnt/beegfs/elmlund/afan/rnasep_klem/expand1/fam1_pipeline/klem_fb/closed_focal_closedfam/rank/conf18_%s_sharpB150.mrc' % name, vs.astype(np.float32), SMPD); print('wrote', name)
PYEOF
