#!/usr/bin/env python3
"""lsamp_swingpos.py: truth check of map readouts on the swing30 sets: every map in each DIR (state_NN.mrc, box 160) placed on the true swing
direction (conf230 - conf200, lp 15 A, inside the sphere, least squares after removing the mean conformer) in A of swing, next to the true
mean swing of the map's particles (DIR/labels.txt, master rows). usage: lsamp_swingpos.py MAPSDIR GT DIR [DIR ...]"""
import sys, glob, os, numpy as np
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
md, gt = sys.argv[1], sys.argv[2]; n = 160; smpd = 1.464; f = np.fft.fftfreq(n, d=smpd); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n, d=smpd), indexing='ij')
G = np.exp(-0.5 * (kx ** 2 + ky ** 2 + kz ** 2) * (15. / 1.2) ** 2); lp = lambda v: np.fft.irfftn(np.fft.rfftn(v) * G, s=v.shape)
g = np.arange(n) - n // 2; R = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2); sph = R < 0.48 * n
T = {k: lp(rd('%s/conf%d_swing_m160.mrc' % (md, 200 + k))) for k in range(31)}; mean_t = sum(T.values()) / 31.; D = T[30] - T[0]
def std(v): v = v - v[sph].mean(); return v / np.sqrt((v[sph] ** 2).mean())
Tz = {k: std(v) for k, v in T.items()}; Dz = Tz[30] - Tz[0]; base = Tz[0]
gtab = np.loadtxt(gt, comments='#'); conf = np.zeros(int(gtab[:, 0].max()) + 1, int); conf[gtab[:, 0].astype(int)] = gtab[:, 1]
for d in sys.argv[3:]:
    lab = np.loadtxt(d + '/labels.txt', dtype=int); cf = conf[1:len(lab) + 1]; out = []
    for fn in sorted(glob.glob(d + '/state_*.mrc')):
        k = int(os.path.basename(fn)[6:8]); v = std(lp(rd(fn))); pos = 30. * np.sum(((v - base) * Dz)[sph]) / np.sum((Dz * Dz)[sph])
        m = lab == k; sw = (cf[m] >= 200) & (cf[m] <= 230); out.append('%d:%.1f/%.1f(%d%%)' % (k, pos, (cf[m][sw] - 200).mean() if sw.any() else np.nan, 100 * sw.mean()))
    print('%-12s map position / true mean swing (swing share), A: %s' % (os.path.basename(d), ' '.join(out)))
