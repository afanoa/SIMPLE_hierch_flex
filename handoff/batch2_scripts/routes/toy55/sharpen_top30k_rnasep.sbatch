#!/bin/bash
#SBATCH -J rnasep_sharpen_top30k2
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# B-sharpened (-150, low-passed at the half-map FSC, unmasked) maps of the family-1-wide docked ranking for download.
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; FD=${FD:-$A/rnasep_klem/expand1/fam1_pipeline/klem_fb/closed_focal}; mkdir -p $FD/lever/for_download
$PY - <<PYEOF
import numpy as np, sys
sys.path.insert(0, '$W'); from model_class import read_mrc, write_mrc, Shells, fftc
SMPD = 1.464; n = 160; S = Shells(n, SMPD); sh = S.sh; zz, yy, xx = np.indices((n, n, n)); rr = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2); s_inv = rr / (n * SMPD)
def fsc_res(x, y):
    a_, b_ = fftc(x), fftc(y); num = np.bincount(sh.ravel(), np.real(np.conj(a_) * b_).ravel()); d1 = np.bincount(sh.ravel(), np.abs(a_.ravel()) ** 2); d2 = np.bincount(sh.ravel(), np.abs(b_.ravel()) ** 2)
    f = num / np.maximum(np.sqrt(d1 * d2), 1e-30); s = 1
    while s < 80 and f[s] > 0.143: s += 1
    return n * SMPD / s
for st, name in (('01', 'conf18_docked_top30k_sharpB150.mrc'), ('02', 'conf18_next20k_sharpB150.mrc'), ('03', 'undocked_bottom30k_sharpB150.mrc')):
    src = '$FD/lever/maps/recvol_state%s.mrc' % st; v = read_mrc(src); F = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v))); F *= np.exp(150. * s_inv ** 2 / 4.)
    res = fsc_res(read_mrc(src.replace('.mrc', '_even.mrc')), read_mrc(src.replace('.mrc', '_odd.mrc'))); fc = 1. / res; F *= 0.5 * (1 + np.cos(np.pi * np.clip((s_inv - 0.85 * fc) / (0.3 * fc), 0, 1)))
    write_mrc('$FD/lever/for_download/' + name, np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(F)))).astype(np.float32), SMPD); print('%s: low-pass %.1f A' % (name, res))
PYEOF
