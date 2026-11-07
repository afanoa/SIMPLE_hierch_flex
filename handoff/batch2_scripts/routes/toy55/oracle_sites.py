#!/usr/bin/env python3
"""Oracle per-particle ceiling of every SITE of simsites: each image (subset of the non-junk rows) is scored at its pose (true or refined)
against the Fourier slices of the TRUE combination maps (CTF, shift, band 2 shells .. LP A, per-image least-squares map scale, per-shell
noise from the residuals at the truth, independent coefficients only = one Hermitian half-plane), and per site s the log-likelihood ratio
g = log p(x | site s in state 1) - log p(x | state 0), the other sites at their true states. Reports per site and band: d' between the
true states, AUC, the Bayes accuracy with the true prior, and the precision of the top-N rows by posterior for the site's MINORITY state
(the ranking readout a site classifier can at best reach). Conventions as splitbench16 / minor_em3 (--conv 1,1,1,1,1,1, signed CTF).
usage: oracle_sites.py --sim SIMDIR --poses POSES --truth sites_truth.txt --maps MAPSDIR --out OUT.npz [--nsub 40000 --lp 10,20 --crop 48]"""
import argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--sim', required=True); ap.add_argument('--poses', required=True); ap.add_argument('--truth', required=True); ap.add_argument('--maps', required=True)
ap.add_argument('--out', required=True); ap.add_argument('--nsub', type=int, default=40000); ap.add_argument('--lp', default='10,20'); ap.add_argument('--crop', type=int, default=48)
ap.add_argument('--mskdiam', type=float, default=230.); ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--kv', type=float, default=300.); ap.add_argument('--cs', type=float, default=2.7)
ap.add_argument('--fraca', type=float, default=0.1); ap.add_argument('--conv', default='1,1,1,1,1,1'); ap.add_argument('--ctfabs', type=int, default=0); ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--topn', default='1000,2500,5000,10000'); a = ap.parse_args(); rng = np.random.default_rng(a.seed)
cv = [float(v) for v in a.conv.split(',')]; convT, s1, yf, s3, ssh, sctf = bool(cv[0]), cv[1], cv[2], cv[3], cv[4], cv[5]
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
hdr = [l.split()[1:] for l in open(a.truth) if l.startswith('#')]; names = hdr[0][3:]; T = np.genfromtxt(a.truth, comments='#'); combo = T[:, 1].astype(int); X = T[:, 3:].astype(int)
cmb = sorted(set(combo[combo > 0].tolist())); cmap = {c: X[np.where(combo == c)[0][0]] for c in cmb}; inv = {tuple(v): c for c, v in cmap.items()}
rows = np.sort(rng.choice(np.where(combo > 0)[0], min(a.nsub, int((combo > 0).sum())), replace=False)); NS = len(rows)
P = np.array([[float(kv(l).get(k, 0)) for k in ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast')] for l in open(a.poses) if l.strip()])[rows]
h0 = np.frombuffer(open(a.sim + '/all.mrcs', 'rb').read(1024), np.int32); n, nz, nsym = int(h0[0]), int(h0[2]), int(h0[23]); c = a.crop
I = np.memmap(a.sim + '/all.mrcs', np.float32, 'r', offset=1024 + nsym, shape=(nz, n, n))
x = np.arange(n) - n // 2; R0 = a.mskdiam / 2 / a.smpd; msk2 = np.clip((R0 + 4 - np.hypot(x[:, None], x[None, :])) / 8., 0, 1).astype(np.float32); sl = np.r_[0:c // 2, n - c // 2:n]
fc = np.fft.fftfreq(c) * c / n; fy2, fx2 = np.meshgrid(fc, fc, indexing='ij'); kk0 = np.stack([fx2, fy2, np.zeros_like(fx2)], -1).reshape(-1, 3) * n
r2 = np.round(np.hypot(kk0[:, 0], kk0[:, 1])).astype(int); half = (kk0[:, 0] > 1e-9) | ((np.abs(kk0[:, 0]) < 1e-9) & (kk0[:, 1] > 1e-9))
lps = [float(v) for v in a.lp.split(',')]; kbmax = min(c // 2 - 1, int(round(n * a.smpd / min(lps)))); cols = np.where((r2 >= 2) & (r2 <= kbmax) & half)[0]; shq = r2[cols]
kf = (np.hypot(fx2, fy2) / a.smpd).ravel()[cols]; ang = np.arctan2(fy2, fx2).ravel()[cols]; lam = 12.2643247 / np.sqrt(a.kv * 1e3 * (1. + a.kv * 1e3 * 0.978466e-6))
def rotm(e1, e2, e3):
    a1, a2, a3 = np.radians(e1), np.radians(e2), np.radians(e3); c1, si1, c2, si2, c3, si3 = np.cos(a1), np.sin(a1), np.cos(a2), np.sin(a2), np.cos(a3), np.sin(a3)
    Z1 = np.zeros((len(e1), 3, 3)); Z1[:, 0, 0] = c1; Z1[:, 0, 1] = -si1; Z1[:, 1, 0] = si1; Z1[:, 1, 1] = c1; Z1[:, 2, 2] = 1
    Y2 = np.zeros((len(e1), 3, 3)); Y2[:, 0, 0] = c2; Y2[:, 0, 2] = si2; Y2[:, 2, 0] = -si2; Y2[:, 2, 2] = c2; Y2[:, 1, 1] = 1
    Z3 = np.zeros((len(e1), 3, 3)); Z3[:, 0, 0] = c3; Z3[:, 0, 1] = -si3; Z3[:, 1, 0] = si3; Z3[:, 1, 1] = c3; Z3[:, 2, 2] = 1
    return Z1 @ Y2 @ Z3
g3 = np.arange(c) - c // 2; r3 = np.sqrt(g3[:, None, None] ** 2 + g3[None, :, None] ** 2 + g3[None, None, :] ** 2) * n / c; msk3 = np.clip((R0 + 4 - r3) / 8., 0, 1)
def vfour(v):   # box-n map -> c^3 Fourier (cropped, masked like the references of minor_em3)
    F = np.fft.fftn(np.fft.ifftshift(v)); vc = np.fft.fftshift(np.real(np.fft.ifftn(F[np.ix_(sl, sl, sl)]))) * (c / n) ** 3
    return np.r_[np.fft.fftn(np.fft.ifftshift(vc * msk3)).ravel(), 0.].astype(np.complex64)
VF = {k: vfour(rd('%s/conf%02d_sites.mrc' % (a.maps, k))) for k in cmb}
Y = np.zeros((NS, len(cols)), np.complex64); CT = np.zeros((NS, len(cols)), np.float32); IX = np.zeros((NS, len(cols)), np.int32)
for b in range(0, NS, 5000):
    rr = rows[b:b + 5000]; pp = P[b:b + 5000]; blk = np.array(I[rr], np.float32) * msk2; F = np.fft.fft2(np.fft.ifftshift(blk, axes=(1, 2)))[:, sl][:, :, sl].reshape(len(rr), -1)[:, cols]
    df = 0.5 * ((pp[:, 5] + pp[:, 6])[:, None] + (pp[:, 5] - pp[:, 6])[:, None] * np.cos(2. * (ang[None] - np.radians(pp[:, 7])[:, None]))) * 1e4
    chi = np.pi * lam * df * kf[None] ** 2 - 0.5 * np.pi * a.cs * 1e7 * lam ** 3 * kf[None] ** 4; C = (sctf * (np.sqrt(1. - a.fraca ** 2) * np.sin(chi) + a.fraca * np.cos(chi))).astype(np.float32)
    if a.ctfabs: C = np.abs(C)
    ramp = np.exp(ssh * 2j * np.pi * (fx2.ravel()[cols][None] * pp[:, 3, None] + fy2.ravel()[cols][None] * pp[:, 4, None])).astype(np.complex64)
    Y[b:b + len(rr)] = F * np.conj(ramp); CT[b:b + len(rr)] = C   # image shifted back to the slice frame (minor_em3: Wv = C conj(ramp) F)
    Rm = rotm(s1 * pp[:, 0], pp[:, 1], s3 * pp[:, 2]); kk = kk0[cols] * np.array([1., yf, 1.])
    p3 = np.einsum('qj,nij->nqi', kk, Rm) if convT else np.einsum('qj,nji->nqi', kk, Rm); q = np.round(p3).astype(np.int32)
    ok = (np.abs(q) <= c // 2 - 1).all(2); q = q % c; ixb = q[..., 2] * c * c + q[..., 1] * c + q[..., 0]; ixb[~ok] = c ** 3; IX[b:b + len(rr)] = ixb
print('%d particles scored, crop %d, %d independent coefficients per image (shells 2..%d)' % (NS, c, len(cols), kbmax), flush=True)
ktrue = combo[rows]; res = {}
for lp in lps:
    kb = int(round(n * a.smpd / lp)); cm = shq <= kb; ns = kb + 1
    # per-shell noise from the residuals at the truth (per-image LS scale)
    def sse(k_of_row, wsh):
        out = np.zeros(NS);
        for k in cmb:
            m = np.where(k_of_row == k)[0]
            for b in range(0, len(m), 5000):
                r_ = m[b:b + 5000]; Pk = VF[k][IX[r_]][:, cm] * CT[r_][:, cm]; Yb = Y[r_][:, cm]; w = wsh[shq[cm]][None]
                beta = np.maximum((np.real(Pk * np.conj(Yb)) * w).sum(1) / np.maximum((np.abs(Pk) ** 2 * w).sum(1), 1e-30), 0.)
                out[r_] = (np.abs(Yb - beta[:, None] * Pk) ** 2 * w).sum(1)
        return out
    num = np.zeros(ns); cnt = np.zeros(ns)
    for k in cmb:
        m = np.where(ktrue == k)[0]
        for b in range(0, len(m), 5000):
            r_ = m[b:b + 5000]; Pk = VF[k][IX[r_]][:, cm] * CT[r_][:, cm]; Yb = Y[r_][:, cm]
            beta = np.maximum(np.real(Pk * np.conj(Yb)).sum(1) / np.maximum((np.abs(Pk) ** 2).sum(1), 1e-30), 0.)
            num += np.bincount(np.tile(shq[cm], len(r_)), (np.abs(Yb - beta[:, None] * Pk) ** 2).ravel(), ns); cnt += np.bincount(np.tile(shq[cm], len(r_)), minlength=ns)
    wsh = 1. / np.maximum(num / np.maximum(cnt, 1), 1e-30)
    for j, nm in enumerate(names):
        st = np.array([cmap[k] for k in ktrue]); s0 = st.copy(); s0[:, j] = 0; s1_ = st.copy(); s1_[:, j] = 1
        k0 = np.array([inv[tuple(v)] for v in s0]); k1 = np.array([inv[tuple(v)] for v in s1_]); g = sse(k0, wsh) - sse(k1, wsh)   # log p(x|1) - log p(x|0) (nats)
        y = st[:, j] == 1; p1 = y.mean(); mu1, mu0 = g[y].mean(), g[~y].mean(); dp = (mu1 - mu0) / np.sqrt(0.5 * (g[y].var() + g[~y].var()))
        o = np.argsort(np.r_[g[y], g[~y]]); rk = np.empty(len(o)); rk[o] = np.arange(len(o)); auc = (rk[:y.sum()].sum() - y.sum() * (y.sum() - 1) / 2) / (y.sum() * (~y).sum())
        post = 1. / (1. + np.exp(-(g + np.log(p1 / (1 - p1))))); acc = np.mean((post > 0.5) == y)
        minor = 1 if p1 < 0.5 else 0; sc = g if minor == 1 else -g; o = np.argsort(-sc); full = len(combo[combo > 0])
        prec = ['top %d: %.3f' % (N_, (st[o[:max(1, int(round(N_ * NS / full)))], j] == minor).mean()) for N_ in [int(v) for v in a.topn.split(',')]]
        res['%s_lp%g' % (nm, lp)] = dict(dprime=float(dp), auc=float(auc), acc=float(acc), prior=float(p1), g=g, y=y)
        print('lp %4.1f site %-9s state-1 share %.3f: d\' %.2f, AUC %.3f, Bayes accuracy %.3f (prior-only %.3f); precision of the minority state (%d) among the top-N of %d rows by posterior: %s' %
              (lp, nm, p1, dp, auc, acc, max(p1, 1 - p1), minor, full, ', '.join(prec)), flush=True)
np.savez(a.out, rows=rows, ktrue=ktrue, **{k: v['g'] for k, v in res.items()})
