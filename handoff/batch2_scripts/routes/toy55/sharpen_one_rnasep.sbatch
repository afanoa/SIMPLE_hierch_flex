#!/bin/bash
#SBATCH -J rnasep_sharpen_one
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# B-sharpen one unmasked map, low-passed at its half-map FSC. env: SRC (recvol_stateNN.mrc with _even/_odd beside it), DST, BFAC (150), LPFIX (0 = the FSC resolution, else a fixed low-pass in A)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; mkdir -p $(dirname $DST)
$PY - <<PYEOF
import numpy as np, sys
sys.path.insert(0, '$W'); from model_class import read_mrc, write_mrc, Shells, fftc
SMPD = 1.464; n = 160; S = Shells(n, SMPD); sh = S.sh; zz, yy, xx = np.indices((n, n, n)); rr = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2); s_inv = rr / (n * SMPD)
def fsc_res(x, y):
    a_, b_ = fftc(x), fftc(y); num = np.bincount(sh.ravel(), np.real(np.conj(a_) * b_).ravel()); d1 = np.bincount(sh.ravel(), np.abs(a_.ravel()) ** 2); d2 = np.bincount(sh.ravel(), np.abs(b_.ravel()) ** 2)
    f = num / np.maximum(np.sqrt(d1 * d2), 1e-30); s = 1
    while s < 80 and f[s] > 0.143: s += 1
    return n * SMPD / s
src = '$SRC'; v = read_mrc(src); F = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v))); F *= np.exp(float('${BFAC:-150}') * s_inv ** 2 / 4.)
res = fsc_res(read_mrc(src.replace('.mrc', '_even.mrc')), read_mrc(src.replace('.mrc', '_odd.mrc'))); res = max(res, float('${LPFIX:-0}')); fc = 1. / res; F *= 0.5 * (1 + np.cos(np.pi * np.clip((s_inv - 0.85 * fc) / (0.3 * fc), 0, 1)))
write_mrc('$DST', np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(F)))).astype(np.float32), SMPD); print('%s: low-pass %.1f A' % ('$DST', res))
PYEOF
