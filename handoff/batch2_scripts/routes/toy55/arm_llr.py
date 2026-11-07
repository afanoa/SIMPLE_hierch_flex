#!/usr/bin/env python3
"""Per-image likelihood ceiling for the hinged-arm set at KNOWN poses: every image of a random subset is scored against
the Fourier slices (CTF, shift, band 60..lp A, residual-whitened) of ALL synthetic conformer maps at its true pose,
whole image and inside the projected region mask. Reports (1) the per-image conformer accuracy (chance 1/nconf),
(2) per arm the per-image log-likelihood ratio truth vs conformers whose arm differs by > FAR deg: mean, std, the
detectability mean/std and the images that must be pooled for 3 sigma. The Euler / shift / CTF sign conventions are
calibrated on the data (mean truth-slice correlation over the subset).
usage: arm_llr.py SIMDIR TRUTH_ARMS MASK SYNTHMAPS OUT [--nsub 3000 --crop 80 --lp 10 --hp 60 --far 30 --smpd 1.464]"""
import sys, os, glob, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc
ap = argparse.ArgumentParser(); ap.add_argument('sim'); ap.add_argument('truth'); ap.add_argument('mask'); ap.add_argument('maps'); ap.add_argument('out')
ap.add_argument('--nsub', type=int, default=3000); ap.add_argument('--crop', type=int, default=80); ap.add_argument('--lp', type=float, default=10.); ap.add_argument('--hp', type=float, default=60.)
ap.add_argument('--far', type=float, default=30.); ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--mskdiam', type=float, default=230.); ap.add_argument('--kv', type=float, default=300.)
ap.add_argument('--cs', type=float, default=2.7); ap.add_argument('--fraca', type=float, default=0.1); ap.add_argument('--seed', type=int, default=0); ap.add_argument('--lsscale', type=int, default=1); ap.add_argument('--oris', default='none'); a = ap.parse_args()   # --oris: pose table to score at (e.g. refined) instead of the simulation truth
os.makedirs(a.out, exist_ok=True); rng = np.random.default_rng(a.seed)
# ---- truth, poses, CTF
G = np.loadtxt(a.sim + '/gt_particles.txt', comments='#'); conf = G[:, 1].astype(int); N = len(conf)
if a.truth != 'none':
    ang = {int(w[0]): np.array([float(x) for x in w[1:]]) for w in (l.split() for l in open(a.truth)) if w and not w[0].startswith('#')}
    confs = sorted(ang); narm = len(ang[confs[0]]) // 2; TH = np.array([ang[c] for c in confs])   # (nconf, 2*narm): theta, phi per arm
else:
    confs = sorted(set(conf[conf > 0].tolist())); narm = 0; TH = None
keys = ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast')
P = np.zeros((N, len(keys)), np.float64)
for i, l in enumerate(open(a.sim + '/all.txt' if a.oris == 'none' else a.oris)):
    kv = dict(t.split('=', 1) for t in l.split() if '=' in t); P[i] = [float(kv[k]) for k in keys]
sel = np.sort(rng.choice(np.where(conf > 0)[0], min(a.nsub, int((conf > 0).sum())), replace=False)); n_s = len(sel)
h = np.frombuffer(open(a.sim + '/all.mrcs', 'rb').read(1024), np.int32); nx, ny, nz, nsym = [int(v) for v in (h[0], h[1], h[2], h[23])]
imgs = np.array(np.memmap(a.sim + '/all.mrcs', np.float32, 'r', offset=1024 + nsym, shape=(nz, ny, nx))[sel], np.float32); n = nx; c = a.crop
print('%d of %d particles (%d conformers, %d arms), box %d -> crop %d, band %.0f..%.0f A' % (n_s, N, len(confs), narm, n, c, a.hp, a.lp), flush=True)
# ---- 2D / 3D Fourier geometry (as avgs_states3d.py: origin-centred transforms, crop to c, shells in units of the full box)
x = np.arange(n) - n // 2; rr2 = np.hypot(x[:, None], x[None, :]); R = a.mskdiam / 2 / a.smpd; msk2 = np.clip((R + 4 - rr2) / 8., 0, 1).astype(np.float32)
sl = np.r_[0:c // 2, n - c // 2:n]
def fft2c(img): return np.fft.fft2(np.fft.ifftshift(img * msk2))[np.ix_(sl, sl)]
fc = np.fft.fftfreq(c) * c / n; fy2, fx2 = np.meshgrid(fc, fc, indexing='ij'); sh2 = np.round(np.hypot(fx2, fy2) * n).astype(int)
klo, khi = int(round(n * a.smpd / a.hp)), int(round(n * a.smpd / a.lp)); band = (sh2 >= klo) & (sh2 <= khi); nsh = khi + 2
F = np.array([fft2c(im) for im in imgs]); del imgs
def rot(e1d, e2d, e3d):
    a1, a2, a3 = np.radians([e1d, e2d, e3d]); c1, s1, c2, s2, c3, s3 = np.cos(a1), np.sin(a1), np.cos(a2), np.sin(a2), np.cos(a3), np.sin(a3)
    return np.array([[c1, -s1, 0], [s1, c1, 0], [0, 0, 1]]) @ np.array([[c2, 0, s2], [0, 1, 0], [-s2, 0, c2]]) @ np.array([[c3, -s3, 0], [s3, c3, 0], [0, 0, 1]])
kk0 = np.stack([fx2, fy2, np.zeros_like(fx2)], -1).reshape(-1, 3) * n
def slice_index(Rm, convT, yflip):
    kk = kk0 * np.array([1., yflip, 1.]); p = kk @ Rm.T if convT else kk @ Rm; q = np.round(p).astype(int)
    ok = (np.abs(q) <= c // 2 - 1).all(1); q = q % c; return q[:, 2] * c * c + q[:, 1] * c + q[:, 0], ok
x3 = np.arange(n) - n // 2; rr3 = np.sqrt(x3[:, None, None] ** 2 + x3[None, :, None] ** 2 + x3[None, None, :] ** 2); msk3 = np.clip((R + 4 - rr3) / 8., 0, 1).astype(np.float32)
sl3 = np.r_[0:c // 2, n - c // 2:n]
def vfft(v): return np.fft.fftn(np.fft.ifftshift(v * msk3))[np.ix_(sl3, sl3, sl3)].ravel().astype(np.complex64)
mapf = lambda k: ('%s/conf%02d_synth.mrc' % (a.maps, k)) if os.path.isfile('%s/conf%02d_synth.mrc' % (a.maps, k)) else sorted(glob.glob('%s/conf%02d_*.mrc' % (a.maps, k)) + glob.glob('%s/conf%02d.mrc' % (a.maps, k)))[0]
VF = np.array([vfft(read_mrc(mapf(k)).astype(np.float32)) for k in confs]); print('conformer maps transformed', flush=True)
MF = vfft(read_mrc(a.mask).astype(np.float32) / max(read_mrc(a.mask).max(), 1e-12))
# ---- CTF (SIMPLE: kv, cs, fraca; defocus in microns, astigmatism angle in degrees) and shifts, with free signs
lam = 12.2643247 / np.sqrt(a.kv * 1e3 * (1. + a.kv * 1e3 * 0.978466e-6)); kfreq = np.hypot(fx2, fy2) / a.smpd   # 1/A
def ctf2d(dfx, dfy, angast, sgn):
    ang = np.arctan2(fy2, fx2); df = 0.5 * ((dfx + dfy) + (dfx - dfy) * np.cos(2. * (ang - np.radians(angast)))) * 1e4
    chi = np.pi * lam * df * kfreq ** 2 - 0.5 * np.pi * a.cs * 1e7 * lam ** 3 * kfreq ** 4
    return sgn * (np.sqrt(1. - a.fraca ** 2) * np.sin(chi) + a.fraca * np.cos(chi))
def ramp(sx, sy, sgn): return np.exp(sgn * 2j * np.pi * (fx2 * sx + fy2 * sy) * (n / n)).astype(np.complex64)   # fx2 in cycles per full box
ktrue = np.array([confs.index(k) for k in conf[sel]])
def model(i, j, conv):   # slice of conformer j for particle i under convention conv
    convT, s1, yf, s3, ssh, sctf = conv; e1, e2, e3, sx, sy, dfx, dfy, aa = P[sel[i]]
    idx, ok = slice_index(rot(s1 * e1, e2, s3 * e3), convT, yf); S = np.where(ok, VF[j][idx], 0).reshape(c, c)
    return S * ctf2d(dfx, dfy, aa, sctf) * ramp(sx, sy, ssh)
def corr_band(A, B, bm): a_, b_ = A[bm], B[bm]; return float(np.real(np.vdot(a_, b_)) / (np.linalg.norm(a_) * np.linalg.norm(b_) + 1e-30))
# ---- convention calibration on the truth slices (low band, mean over particles)
band_lo = (sh2 >= 2) & (sh2 <= khi); ncal = min(400, n_s); best = None
for convT in (False, True):
    for s1 in (1., -1.):
        for yf in (1., -1.):
            for s3 in (1., -1.):
                for ssh in (1., -1.):
                    for sctf in (1., -1.):
                        conv = (convT, s1, yf, s3, ssh, sctf); cc = np.mean([corr_band(model(i, ktrue[i], conv), F[i], band_lo) for i in range(ncal)])
                        if best is None or cc > best[0]: best = (cc, conv)
conv = best[1]; print('convention: transpose=%s e1sign=%+d yflip=%+d e3sign=%+d shiftsign=%+d ctfsign=%+d (mean truth-slice corr %.4f over %d particles)' % (conv[0], int(conv[1]), int(conv[2]), int(conv[3]), int(conv[4]), int(conv[5]), best[0], ncal), flush=True)
for conv_alt in [(conv[0], conv[1], conv[2], conv[3], conv[4], -conv[5]), (conv[0], conv[1], conv[2], conv[3], -conv[4], conv[5]), (conv[0], conv[1], conv[2], -conv[3], conv[4], conv[5])]:
    print('  flipping one sign (ctf/shift/e3): corr %.4f' % np.mean([corr_band(model(i, ktrue[i], conv_alt), F[i], band_lo) for i in range(ncal)]))
# ---- residual whitening from the truth slices, then scores of every conformer, whole image and inside the projected mask
num = np.zeros(nsh); den = np.zeros(nsh)
for i in range(n_s):
    r = np.abs(F[i] - model(i, ktrue[i], conv)) ** 2; num += np.bincount(sh2.ravel(), r.ravel(), minlength=nsh)[:nsh]; den += np.bincount(sh2.ravel(), None, minlength=nsh)[:nsh]
sig2 = num / np.maximum(den, 1); W = np.where(band, 1. / np.sqrt(np.maximum(sig2[np.minimum(sh2, nsh - 1)], 1e-30)), 0.).astype(np.float32)
sig_sig = num / np.maximum(den, 1)
print('residual power per shell (truth): %s' % ' '.join('%.2g' % v for v in sig2[klo:khi + 1:max(1, (khi - klo) // 8)]), flush=True)
SSE = np.zeros((n_s, len(confs)), np.float32); SSEm = np.zeros((n_s, len(confs)), np.float32); SC = np.ones((n_s, len(confs)), np.float32); mfrac = np.zeros(n_s)
for i in range(n_s):
    e1, e2, e3, sx, sy, dfx, dfy, aa = P[sel[i]]; idx, ok = slice_index(rot(conv[1] * e1, e2, conv[3] * e3), conv[0], conv[2])
    cph = ctf2d(dfx, dfy, aa, conv[5]) * ramp(sx, sy, conv[4]); Fw = F[i] * W
    pm = np.real(np.fft.ifft2(np.where(ok, MF[idx], 0).reshape(c, c) * ramp(sx, sy, conv[4]))); pm = np.clip(pm / max(pm.max(), 1e-12), 0, 1); mfrac[i] = pm.mean()
    for j in range(len(confs)):
        Pw = np.where(ok, VF[j][idx], 0).reshape(c, c) * cph * W   # per-image LS amplitude of the map (a free scale, as in any refinement), else the map as is
        if a.lsscale: SC[i, j] = max(np.real(np.vdot(Pw, Fw)) / max(np.real(np.vdot(Pw, Pw)), 1e-30), 0.)
        D = Fw - SC[i, j] * Pw; SSE[i, j] = 0.5 * np.sum(np.abs(D) ** 2) / (c * c)
        d = np.fft.ifft2(D); SSEm[i, j] = 0.5 * np.sum(np.abs(d * pm) ** 2)
    if i % 500 == 0: print('  scored %d / %d' % (i, n_s), flush=True)
np.savez(a.out + '/arm_llr.npz', sel=sel, ktrue=ktrue, SSE=SSE, SSEm=SSEm, SC=SC, TH=TH, confs=confs, mfrac=mfrac)
print('LS scale at the truth: median %.3f, 10-90 pct %.3f-%.3f (lsscale %d)' % (np.median(SC[np.arange(n_s), ktrue]), np.percentile(SC[np.arange(n_s), ktrue], 10), np.percentile(SC[np.arange(n_s), ktrue], 90), a.lsscale))
def report(S, tag):
    top = np.bincount(np.argmin(S, 1), minlength=len(confs)); print('%s: most-chosen conformer takes %.2f of the images; spread of the per-conformer mean score over conformers %.3g vs over images %.3g' % (tag, top.max() / n_s, S.mean(0).std(), S.std(0).mean()))
    acc = np.mean(np.argmin(S, 1) == ktrue); print('%s: per-image conformer accuracy %.3f (chance %.3f), rank of the truth median %d of %d' % (tag, acc, 1. / len(confs), np.median([(S[i] < S[i, ktrue[i]]).sum() + 1 for i in range(n_s)]), len(confs)))
    LLR = S - S[np.arange(n_s), ktrue][:, None]   # positive = the truth fits better (nats, whitened)
    dp = []
    for t in range(len(confs)):
        rows = np.where(ktrue == t)[0]
        if len(rows) < 8: continue
        for j in range(len(confs)):
            if j == t: continue
            v = LLR[rows, j]; dp.append(v.mean() / max(v.std(), 1e-12))
    dp = np.array(dp); print('  conformer pairs: per-image d\' median %.3f, 10th pct %.3f -> images for 3 sigma: %d (median pair), %d (hard pairs)' % (np.median(dp), np.percentile(dp, 10), int(np.ceil((3 / max(np.median(dp), 1e-9)) ** 2)), int(np.ceil((3 / max(np.percentile(dp, 10), 1e-9)) ** 2))))
    for m in (1, 10, 100, 1000):
        hit = []
        for t in range(len(confs)):
            rows = np.where(ktrue == t)[0]
            if len(rows) < 8: continue
            for _ in range(20):
                g = rng.choice(rows, m, replace=True); hit.append(int(np.argmin(S[g].sum(0)) == t))
        print('  pooled %4d images of one conformer: conformer accuracy %.2f (chance %.3f)' % (m, np.mean(hit), 1. / len(confs)))
    for k in range(narm):
        dth = np.abs(TH[:, 2 * k][None, :] - TH[ktrue, 2 * k][:, None]); far = dth > a.far; near = (dth < 5.) & (np.arange(len(confs))[None, :] != ktrue[:, None])
        li = np.array([LLR[i, far[i]].mean() if far[i].any() else np.nan for i in range(n_s)]); mu, sd = np.nanmean(li), np.nanstd(li)
        ln = np.array([LLR[i, near[i]].mean() if near[i].any() else np.nan for i in range(n_s)])
        print('  arm %d: LLR truth vs arm-%d-far (> %.0f deg) per image %.3g +- %.3g (d\' %.3f, %d images for 3 sigma); vs arm-%d-near (< 5 deg, other arms free) %.3g' % (k + 1, k + 1, a.far, mu, sd, mu / max(sd, 1e-12), int(np.ceil((3 * sd / max(abs(mu), 1e-12)) ** 2)), k + 1, np.nanmean(ln)))
print('projected mask covers %.3f of the image area (mean)' % mfrac.mean())
report(SSE, 'WHOLE image'); report(SSEm, 'MASKED to the region')
