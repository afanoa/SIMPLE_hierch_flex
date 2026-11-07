#!/usr/bin/env python3
"""Real-space variation mask of a set of maps (e.g. children ordered along a motion): each map's even / odd halves low-passed to LP A and
scaled (least squares inside the support) to the weighted mean of the set; v(x) = sum_j w_j (E_j - E)(O_j - O) / sum w (cross-half:
noise cancels), Gaussian-smoothed (SIG px); mask = voxels above FRAC of the maximum inside the sphere, dilated DIL px, soft edge SOFT px.
Optional readout: share of the closed-arm / open-arm regions (arm_readout.py's definition from CONS and CLOSED) inside the mask.
usage: arm_mask.py OUT.mrc --list LIST (name even odd particles) [--lp 15 --frac 0.2 --sig 1.5 --dil 3 --soft 3 --mskdiam 232
       --cons consensus.mrc --closed closed.mrc]"""
import argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--list', required=True); ap.add_argument('--lp', type=float, default=15.)
ap.add_argument('--frac', type=float, default=0.2); ap.add_argument('--sig', type=float, default=1.5); ap.add_argument('--dil', type=int, default=3)
ap.add_argument('--soft', type=float, default=3.); ap.add_argument('--mskdiam', type=float, default=232.); ap.add_argument('--cons', default=''); ap.add_argument('--closed', default='')
a = ap.parse_args()
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
def wr(p, x, smpd):
    h = np.zeros(256, np.int32); h[0:3] = x.shape; h[3] = 2; h[7:10] = x.shape; hf = h.view(np.float32); hf[10:13] = np.array(x.shape) * smpd; hf[13:16] = 90.
    h[16:19] = [1, 2, 3]; hf[19], hf[20], hf[21] = x.min(), x.max(), x.mean(); h[52] = int.from_bytes(b'MAP ', 'little'); h[53] = 0x4444
    with open(p, 'wb') as f: f.write(h.tobytes()); f.write(x.astype(np.float32).tobytes())
L = [l.split() for l in open(a.list) if l.strip() and not l.startswith('#')]; w = np.array([float(l[3]) for l in L]); w /= w.sum()
n = int(np.fromfile(L[0][1], np.int32, 1)[0]); smpd = float(np.fromfile(L[0][1], np.float32, 13)[10]) / n
f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n), indexing='ij'); k = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2) / smpd
lpf = np.exp(-(k * a.lp) ** 2 * np.log(2)); lpr = lambda x: np.fft.irfftn(np.fft.rfftn(x) * lpf, s=x.shape)
g = np.arange(n) - n // 2; R = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2); sph = R <= a.mskdiam / 2 / smpd
E = np.array([lpr(rd(l[1])) for l in L]); O = np.array([lpr(rd(l[2])) for l in L]); Em, Om = np.tensordot(w, E, 1), np.tensordot(w, O, 1)
sup = ((Em + Om) / 2 > 0.1 * ((Em + Om) / 2)[sph].max()) & sph
for j in range(len(L)):   # least-squares scale of every map to the mean inside the support
    E[j] *= (E[j][sup] @ Em[sup]) / max(E[j][sup] @ E[j][sup], 1e-30); O[j] *= (O[j][sup] @ Om[sup]) / max(O[j][sup] @ O[j][sup], 1e-30)
Em, Om = np.tensordot(w, E, 1), np.tensordot(w, O, 1); v = np.tensordot(w, (E - Em) * (O - Om), 1)
gw = np.minimum(np.arange(n), n - np.arange(n)); Rw = np.sqrt(gw[:, None, None] ** 2 + gw[None, :, None] ** 2 + gw[None, None, :] ** 2)
conv = lambda x, kk: np.fft.irfftn(np.fft.rfftn(x) * np.fft.rfftn(kk / kk.sum()), s=x.shape)
vs = conv(v, np.exp(-Rw ** 2 / (2 * a.sig ** 2))) * sph; core = vs > a.frac * vs.max()
m = (conv(core.astype(float), (Rw <= a.dil).astype(float)) > 0.5 / (Rw <= a.dil).sum()).astype(float); m = np.clip(conv(m, np.exp(-Rw ** 2 / (2 * a.soft ** 2))), 0, 1) * sph
wr(a.out, m, smpd); wr(a.out.replace('.mrc', '_variance.mrc'), vs, smpd)
line = 'real-space variation mask %s: %.1f %% of the molecule support (%.1f %% of the sphere)' % (a.out, 100 * ((m > 0.5) & sup).sum() / sup.sum(), 100 * (m > 0.5).sum() / sph.sum())
if a.cons and a.closed:   # arm regions as arm_readout.py (10 A, standardised, consensus vs deposited closed), resampled to this box
    def crop_to(x):
        nb = x.shape[0]
        if nb == n: return x
        F = np.fft.fftshift(np.fft.fftn(x)); lo = nb // 2 - n // 2; return np.real(np.fft.ifftn(np.fft.ifftshift(F[lo:lo + n, lo:lo + n, lo:lo + n]))) * (n / nb) ** 3
    G10 = np.exp(-0.5 * (k * (10. / 1.2)) ** 2); l10 = lambda x: np.fft.irfftn(np.fft.rfftn(x) * G10, s=x.shape)
    cs = l10(crop_to(rd(a.cons))); cl = l10(crop_to(rd(a.closed))); s2 = cs > 0.2 * cs.max(); out_ = cs < 0.02 * cs.max()
    std = lambda x: (x - x[out_].mean()) / max(np.sqrt(np.mean((x - x[out_].mean())[s2] ** 2)), 1e-12); d = std(cl) - std(cs)
    thr = 0.5 * np.percentile(np.abs(d[s2 | (cl > 0.2 * cl.max())]), 99); mC, mO = d > thr, d < -thr; mb = m > 0.5
    line += '; closed-arm region covered %.2f, open-arm region covered %.2f, share of the mask on the arm regions %.2f' % ((mb & mC).sum() / max(mC.sum(), 1), (mb & mO).sum() / max(mO.sum(), 1), (mb & (mC | mO)).sum() / max(mb.sum(), 1))
print(line)
