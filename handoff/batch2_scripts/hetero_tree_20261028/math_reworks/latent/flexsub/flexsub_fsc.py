#!/usr/bin/env python3
"""flexsub_fsc.py: the sharpness of a stage B's state references: per state the even/odd FSC of the unfiltered half maps inside a
soft sphere of the mask diameter, its 0.5 / 0.143 crossings and its mean over 40-20, 20-15 and 15-10 A, beside the state's size.
usage: flexsub_fsc.py MSKDIAM NAME=STAGEB_DIR [NAME=STAGEB_DIR ...]"""
import glob, os, re, sys
import numpy as np
sys.path.insert(0, os.environ.get('PPINFO', '/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo'))
from mrcio import read_mrc


def soft_sphere(n, rad, hw=6.):
    x = np.arange(n) - n // 2
    r = np.sqrt(x[:, None, None] ** 2 + x[None, :, None] ** 2 + x[None, None, :] ** 2)
    w = np.clip(0.5 * (1. + np.cos(np.pi * (r - (rad - hw)) / (2. * hw))), 0., 1.)
    w[r <= rad - hw] = 1.
    w[r >= rad + hw] = 0.
    return w


def fsc(a, b):
    n = a.shape[0]
    A, B_ = np.fft.fftn(a), np.fft.fftn(b)
    f = np.fft.fftfreq(n) * n
    s = np.rint(np.sqrt(f[:, None, None] ** 2 + f[None, :, None] ** 2 + f[None, None, :] ** 2)).astype(int).ravel()
    num = np.bincount(s, (A * B_.conj()).real.ravel(), n)
    den = np.sqrt(np.bincount(s, (np.abs(A) ** 2).ravel(), n) * np.bincount(s, (np.abs(B_) ** 2).ravel(), n))
    return (num / np.maximum(den, 1e-30))[1:n // 2]


def crossing(c, res, t):
    for i in range(len(c)):
        if c[i] < t:
            return res[i - 1] if i > 0 else res[0]
    return res[-1]


msk = float(sys.argv[1])
for arg in sys.argv[2:]:
    name, d = arg.split('=', 1)
    for e in sorted(glob.glob(os.path.join(d, 'recvol_state??_even_unfil.mrc'))):
        st = int(re.search(r'state(\d\d)', e).group(1))
        o = e.replace('_even_', '_odd_')
        if not os.path.exists(o):
            continue
        ve, smpd = read_mrc(e)
        vo, _ = read_mrc(o)
        n = ve.shape[0]
        m = soft_sphere(n, msk / 2. / smpd)
        c = fsc(ve * m, vo * m)
        res = n * smpd / np.arange(1, n // 2)
        band = lambda lo, hi: float(c[(res <= lo) & (res > hi)].mean())
        print('%-10s state %d: FSC 0.5 %.2f A, 0.143 %.2f A; mean FSC 40-20 A %.4f, 20-15 A %.4f, 15-10 A %.4f'
              % (name, st, crossing(c, res, 0.5), crossing(c, res, 0.143), band(40, 20), band(20, 15), band(15, 10)))
