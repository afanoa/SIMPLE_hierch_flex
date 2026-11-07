#!/usr/bin/env python3
"""lsamp_score.py (landscape sampling, step 1): every tree particle scored against every leaf map at its own leaf-node pose, against the
half-map it did not contribute to (eo 0 -> the leaf's odd map, eo 1 -> its even map): Gaussian log-likelihood ratio vs the zero model over
shells 2..kb (kb = EM_LIM), per-shell reference scale fitted on the own-leaf references, per-shell noise from the own-leaf residuals
(score_embed.py's model). Also the NOISE covariance of the score vector, G(k, l) = sum_f w_f CTF_f^2 Re(P_k conj P_l) (the covariance a
particle with no heterogeneity would show: its noise projected on the reference slices), averaged over particles per eo half.
Layouts: --sim S_A,S_B --mr MRA,MRB (half stacks; meta_rows = master row of every half row; caches CP_A_c48.npy / CP_B_c48.npy) or
--sim S --mr full (one stack in master-row order, cache CP_A_c48.npy holds every row).
Writes OUT/scores.npy (master rows x K, float32, nan outside the tree) and OUT/meta.npz (eo, own, names, counts, gram0, gram1, em_lim).
usage: lsamp_score.py --sim .. --mr .. --cache CP --poses POSES (one key=value line per master row) --leaves LIST (name even odd count)
       --labels master_labels (leaf 1..K, 0 = none) --out DIR [--em_lim 20 --ctfabs 0 --exclude 3,7 (leaf ids not scored against)]"""
import os, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--sim', required=True); ap.add_argument('--mr', required=True); ap.add_argument('--cache', required=True)
ap.add_argument('--poses', required=True); ap.add_argument('--leaves', required=True); ap.add_argument('--labels', required=True); ap.add_argument('--out', required=True)
ap.add_argument('--em_lim', type=float, default=20.); ap.add_argument('--crop', type=int, default=48); ap.add_argument('--smpd', type=float, default=1.464)
ap.add_argument('--kv', type=float, default=300.); ap.add_argument('--cs', type=float, default=2.7); ap.add_argument('--fraca', type=float, default=0.1)
ap.add_argument('--ctfabs', type=int, default=0); ap.add_argument('--mskdiam', type=float, default=230.); ap.add_argument('--ngram', type=int, default=1500, help='particles per block used for the noise covariance')
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
def log(*s): print(*s, flush=True)
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
lam = 12.2643247 / np.sqrt(a.kv * 1e3 * (1. + a.kv * 1e3 * 0.978466e-6))
def rotm(e1, e2, e3):
    a1, a2, a3 = np.radians(e1), np.radians(e2), np.radians(e3); c1, si1, c2, si2, c3, si3 = np.cos(a1), np.sin(a1), np.cos(a2), np.sin(a2), np.cos(a3), np.sin(a3)
    Z1 = np.zeros((len(e1), 3, 3)); Z1[:, 0, 0] = c1; Z1[:, 0, 1] = -si1; Z1[:, 1, 0] = si1; Z1[:, 1, 1] = c1; Z1[:, 2, 2] = 1
    Y2 = np.zeros((len(e1), 3, 3)); Y2[:, 0, 0] = c2; Y2[:, 0, 2] = si2; Y2[:, 2, 0] = -si2; Y2[:, 2, 2] = c2; Y2[:, 1, 1] = 1
    Z3 = np.zeros((len(e1), 3, 3)); Z3[:, 0, 0] = c3; Z3[:, 0, 1] = -si3; Z3[:, 1, 0] = si3; Z3[:, 1, 1] = c3; Z3[:, 2, 2] = 1
    return Z1 @ Y2 @ Z3
SIMS = a.sim.split(','); full = a.mr == 'full'
n = int(np.frombuffer(open(SIMS[0] + '/all.mrcs', 'rb').read(1024), np.int32)[0]); c = a.crop
fc = np.fft.fftfreq(c) * c / n; fy2, fx2 = np.meshgrid(fc, fc, indexing='ij'); kk0 = np.stack([fx2, fy2, np.zeros_like(fx2)], -1).reshape(-1, 3) * n
f3 = np.fft.fftfreq(c) * c / n; kx, ky, kz = np.meshgrid(f3, f3, f3, indexing='ij'); sh = np.minimum(np.round(np.sqrt(kx ** 2 + ky ** 2 + kz ** 2) * n).astype(int), c // 2).ravel()
ns = c // 2 + 1; kb = min(c // 2 - 1, int(round(n * a.smpd / a.em_lim))); r2 = np.round(np.hypot(kk0[:, 0], kk0[:, 1])).astype(int); cols = np.where((r2 >= 2) & (r2 <= kb))[0]
kf = (np.hypot(fx2, fy2) / a.smpd).ravel()[cols]; ang = np.arctan2(fy2, fx2).ravel()[cols]
LAB = np.loadtxt(a.labels, dtype=int); NM = len(LAB); L = [l.split() for l in open(a.leaves) if l.strip() and not l.startswith('#')]; K = len(L)
def crop_four(v):   # box-n real-space map -> the c^3 FFT-layout Fourier map of this frame (as score_embed.py)
    nb = v.shape[0]; sb = np.r_[0:c // 2, nb - c // 2:nb]
    F = np.fft.fftn(np.fft.ifftshift(v)); vc = np.fft.fftshift(np.real(np.fft.ifftn(F[np.ix_(sb, sb, sb)]))); return np.r_[np.fft.fftn(np.fft.ifftshift(vc)).ravel(), 0.].astype(np.complex64)
REF = [(crop_four(rd(l[1])), crop_four(rd(l[2]))) for l in L]
log('%d leaf maps cropped to %d^3; score band shells 2..%d (%.1f A), %d pixels per image; %s layout' % (K, c, kb, n * a.smpd / kb, len(cols), 'full-stack' if full else 'half-stack'))
P_all = []; EO_all = np.full(NM, -1, int)
for i, l in enumerate(open(a.poses)):   # one line per master row
    d = kv(l); P_all.append([float(d.get(k, 0)) for k in ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast')]); EO_all[i] = int(float(d.get('eo', -1)))
P_all = np.array(P_all); assert len(P_all) == NM, 'poses %d lines vs %d label rows' % (len(P_all), NM)
SC = np.full((NM, K), np.nan, np.float32); GR = [np.zeros((K, K)), np.zeros((K, K))]; GN = [0, 0]
halves = [('A', np.arange(NM))] if full else [(t, np.loadtxt(a.mr.split(',')[h]).astype(int) - 1) for h, t in enumerate('AB')]
def build_cache(simdir, f):   # every particle of the stack: soft disc of mskdiam, 2-D FFT, the c x c low-frequency block (splitbench's cache)
    h = np.frombuffer(open(simdir + '/all.mrcs', 'rb').read(1024), np.int32); nx, nz, nsym = int(h[0]), int(h[2]), int(h[23])
    I = np.memmap(simdir + '/all.mrcs', np.float32, 'r', offset=1024 + nsym, shape=(nz, nx, nx)); x = np.arange(nx) - nx // 2; rr = np.hypot(x[:, None], x[None, :])
    m2 = np.clip((a.mskdiam / 2 / a.smpd + 4 - rr) / 8., 0, 1).astype(np.float32); sl = np.r_[0:c // 2, nx - c // 2:nx]; out = np.lib.format.open_memmap(f, 'w+', np.complex64, (nz, c, c))
    for b in range(0, nz, 4000): blk = np.array(I[b:b + 4000], np.float32) * m2; out[b:b + len(blk)] = np.fft.fft2(np.fft.ifftshift(blk, axes=(1, 2)))[:, sl][:, :, sl]
    out.flush(); del out; log('built cache %s (%d particles)' % (f, nz))
for h_, (tag, mr) in enumerate(halves):   # mr[r] = master row of cache row r
    cf = '%s_%s_c%d.npy' % (a.cache, tag, c)
    if not os.path.isfile(cf): build_cache(SIMS[h_], cf)
    rows = np.where(LAB[mr] > 0)[0]; Fc = np.load(cf, mmap_mode='r'); N = len(rows); own = LAB[mr[rows]] - 1
    P = P_all[mr[rows]]; EO = EO_all[mr[rows]]; bad = (EO < 0).sum()
    if bad: log('  %d tree particles without an eo flag: scored against the even maps' % bad); EO = np.maximum(EO, 0)
    S = np.zeros((N, K), np.float32)
    for b in range(0, N, 20000):
        rr = rows[b:b + 20000]; pp = P[b:b + 20000]; F = np.asarray(Fc[rr]).reshape(len(rr), -1)[:, cols]
        df = 0.5 * ((pp[:, 5] + pp[:, 6])[:, None] + (pp[:, 5] - pp[:, 6])[:, None] * np.cos(2. * (ang[None] - np.radians(pp[:, 7])[:, None]))) * 1e4
        chi = np.pi * lam * df * kf[None] ** 2 - 0.5 * np.pi * a.cs * 1e7 * lam ** 3 * kf[None] ** 4; C = (np.sqrt(1. - a.fraca ** 2) * np.sin(chi) + a.fraca * np.cos(chi)).astype(np.float32)
        if a.ctfabs: C = np.abs(C)
        ramp = np.exp(2j * np.pi * (fx2.ravel()[cols][None] * pp[:, 3, None] + fy2.ravel()[cols][None] * pp[:, 4, None])).astype(np.complex64)
        Wv = C * np.conj(ramp) * F; Wd = C ** 2; F2 = (np.abs(F) ** 2).astype(np.float32)
        Rm = rotm(pp[:, 0], pp[:, 1], pp[:, 2]); q = np.round(np.einsum('qj,nij->nqi', kk0[cols], Rm)).astype(np.int32)
        ok = (np.abs(q) <= c // 2 - 1).all(2); q = q % c; ix = q[..., 2] * c * c + q[..., 1] * c + q[..., 0]; ix[~ok] = c ** 3; shq = np.minimum(np.r_[sh, 0][ix], ns - 1)
        eo = EO[b:b + 20000]; ob = own[b:b + 20000]; Pown = np.zeros_like(Wv)
        for g in range(K):
            for e_ in (0, 1):
                mm = (ob == g) & (eo == e_)
                if mm.any(): Pown[mm] = REF[g][1 - e_][ix[mm]]
        bn = np.bincount(shq.ravel(), np.real(Pown * np.conj(Wv)).ravel(), ns); bd = np.bincount(shq.ravel(), (Wd * np.abs(Pown) ** 2).ravel(), ns)
        beta = np.clip(bn / np.maximum(bd, 1e-30), 0., 4.).astype(np.float32); bpx = beta[shq]
        res = F2 - 2. * np.real(bpx * Pown * np.conj(Wv)) + Wd * (bpx * np.abs(Pown)) ** 2
        wsh = (1. / np.maximum(np.bincount(shq.ravel(), res.ravel(), ns) / np.maximum(np.bincount(shq.ravel(), minlength=ns), 1), 1e-12)).astype(np.float32); wi = wsh[shq]
        ng = min(a.ngram, len(rr)); PG = np.zeros((ng, len(cols), K), np.complex64)
        for g in range(K):
            Pg = np.where((eo == 0)[:, None], REF[g][1][ix], REF[g][0][ix]) * bpx
            S[b:b + len(rr), g] = 0.5 * ((2. * np.real(Pg * np.conj(Wv)) - Wd * np.abs(Pg) ** 2) * wi).sum(1); PG[:, :, g] = Pg[:ng]
        wq = (wi * Wd)[:ng]   # noise covariance of the scores of these particles (per eo half)
        for e_ in (0, 1):
            m = eo[:ng] == e_
            if m.any():
                X = PG[m] * np.sqrt(wq[m])[:, :, None]; GR[e_] += np.einsum('npk,npl->kl', X.real, X.real) + np.einsum('npk,npl->kl', X.imag, X.imag); GN[e_] += int(m.sum())
    SC[mr[rows]] = S
    rank = (S > S[np.arange(N), own][:, None]).sum(1)
    log('half %s: %d particles scored against %d leaves; own leaf ranks first for %.1f %% (top 3 %.1f %%; chance %.1f %%)' % (tag, N, K, 100 * (rank == 0).mean(), 100 * (rank < 3).mean(), 100. / K))
G0, G1 = GR[0] / max(GN[0], 1), GR[1] / max(GN[1], 1)
np.save(a.out + '/scores.npy', SC); own_all = np.where(LAB > 0, LAB - 1, -1)
np.savez(a.out + '/meta.npz', eo=EO_all, own=own_all, names=np.array([l[0] for l in L]), counts=np.array([int(l[3]) for l in L]), gram0=G0, gram1=G1, em_lim=a.em_lim, kb=kb)
# empirical vs noise covariance of the centred scores (signal = excess over noise)
tr = np.isfinite(SC).all(1); Cn = np.eye(K) - 1. / K
for e_, G in ((0, G0), (1, G1)):
    m = tr & (EO_all == e_); X = SC[m] - SC[m].mean(1, keepdims=True); X -= X.mean(0); Cs = X.T @ X / max(m.sum() - 1, 1); Gc = Cn @ G @ Cn
    log('eo %d: %d particles; trace of the centred score covariance %.4g, of the noise covariance %.4g (ratio %.2f)' % (e_, m.sum(), np.trace(Cs), np.trace(Gc), np.trace(Cs) / max(np.trace(Gc), 1e-30)))
log('wrote %s/scores.npy (%d x %d) and meta.npz' % (a.out, NM, K))
