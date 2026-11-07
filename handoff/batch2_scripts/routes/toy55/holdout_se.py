#!/usr/bin/env python3
"""usage: holdout_se.py name:mask.mrc:curve.txt[:prev_name] ... [--em_lim 15 --margin 2 --crop 48 --box 160 --smpd 1.464 --mskdiam 232]
SE of site_size.py's hold-out peak: Fisher SE = (1 - r^2) / sqrt(n_eff - 3), n_eff = (hold-out band coefficients / 2, Hermitian) x
V_eff / c^3, V_eff = (sum M^2)^2 / sum M^4 = the participation volume of the site mask M on the crop grid (both maps are multiplied by M);
plus each round's peak gain over the round it refines in units of sqrt(SE_new^2 + SE_old^2)."""
import sys, numpy as np
kw = {'em_lim': 15., 'margin': 2, 'crop': 48, 'box': 160, 'smpd': 1.464, 'mskdiam': 232.}; args = []
it = iter(sys.argv[1:])
for x in it:
    if x.startswith('--'): kw[x[2:]] = type(kw[x[2:]])(next(it))
    else: args.append(x)
c, n = int(kw['crop']), int(kw['box']); kb = min(c // 2 - 1, int(round(n * kw['smpd'] / kw['em_lim']))); kt = c // 2 - 1; lo = kb + int(kw['margin'])
f3 = np.fft.fftfreq(c) * c / n; kx, ky, kz = np.meshgrid(f3, f3, f3, indexing='ij'); sh = np.minimum(np.round(np.sqrt(kx ** 2 + ky ** 2 + kz ** 2) * n).astype(int), c // 2)
nband = int(((sh >= lo) & (sh <= kt)).sum()); sl = np.r_[0:c // 2, n - c // 2:n]
g3 = np.arange(c) - c // 2; r3 = np.sqrt(g3[:, None, None] ** 2 + g3[None, :, None] ** 2 + g3[None, None, :] ** 2) * n / c; R0 = kw['mskdiam'] / 2 / kw['smpd']; msk3 = np.clip((R0 + 4 - r3) / 8., 0, 1)
def rd(p): h = np.fromfile(p, np.int32, 256); m = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=m ** 3).reshape(m, m, m).astype(np.float64)
def crop_real(v): F = np.fft.fftn(np.fft.ifftshift(v)); return np.fft.fftshift(np.real(np.fft.ifftn(F[np.ix_(sl, sl, sl)]))) * (c / v.shape[0]) ** 3
print('hold-out band shells %d..%d of crop %d (%.1f-%.1f A): %d coefficients' % (lo, kt, c, n * kw['smpd'] / lo, n * kw['smpd'] / kt, nband))
R = {}
for a_ in args:
    p = a_.split(':'); nm, mk, cv = p[0], p[1], p[2]; prev = p[3] if len(p) > 3 else ''
    M = np.clip(crop_real(rd(mk)), 0, 1) * msk3; veff = (M ** 2).sum() ** 2 / max((M ** 4).sum(), 1e-30); neff = 0.5 * nband * veff / c ** 3
    C = np.loadtxt(cv); j = int(np.argmax(C[:, 2])); r = C[j, 2]; se = (1 - r ** 2) / np.sqrt(max(neff - 3, 1.)); R[nm] = (r, se)
    s = '%-10s mask participation volume %.0f voxels (%.2f %% of the crop box), n_eff %.0f; hold-out peak %.3f at N %d, SE %.3f' % (nm, veff, 100 * veff / c ** 3, neff, r, C[j, 0], se)
    if prev in R: r0, s0 = R[prev]; s += '; gain over %s %+.3f = %+.1f combined SE' % (prev, r - r0, (r - r0) / np.sqrt(se ** 2 + s0 ** 2))
    print(s)
