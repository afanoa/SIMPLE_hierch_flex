#!/usr/bin/env python3
"""land_pscore.py: per-particle log-likelihoods against K path maps with the pose integrated locally per map. Gaussian likelihood of each
image (crop CROP, shells 2..kb = EM_LIM, per-shell scale beta and per-shell noise from the mean path map, per-particle noise level, |CTF|
for phase-flipped stacks) at its tree pose theta0 and, per map, its quadratic expansion in a 5-D pose perturbation delta (3 small rotations
about the particle axes, 2 shifts): l(delta) ~ l0 + g.delta - 1/2 delta' H delta, with the rotation derivatives from the map's Fourier
gradient (central-slice) and the shift derivatives from the phase ramp. With a Gaussian pose prior (ROT_SD deg, SH_SD px):
  lprof = l0 + 1/2 g' (H + S^-1)^-1 g             (pose optimised under the prior)
  lmarg = lprof - 1/2 log det(I + S H)             (pose marginalised: the Laplace evidence)
Also the pixel-split (kx + ky parity) contrasts of three map pairs and their model variances (truth-free likelihood temperature).
usage: land_pscore.py --sim S | SA,SB --mr full | MRA,MRB --cache CP --poses P --rowflags F --maps LIST --out OUT.npz [--ctfabs 0 --em_lim 10]"""
import os, argparse, numpy as np
ap = argparse.ArgumentParser()
for k in ('sim', 'mr', 'cache', 'poses', 'rowflags', 'maps', 'out'): ap.add_argument('--' + k, required=True)
ap.add_argument('--ctfabs', type=int, default=0); ap.add_argument('--em_lim', type=float, default=10.); ap.add_argument('--crop', type=int, default=48); ap.add_argument('--mskdiam', type=float, default=230.)
ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--kv', type=float, default=300.); ap.add_argument('--cs', type=float, default=2.7); ap.add_argument('--fraca', type=float, default=0.1)
ap.add_argument('--rot_sd', type=float, default=2.5); ap.add_argument('--sh_sd', type=float, default=1.0); ap.add_argument('--block', type=int, default=8000); a = ap.parse_args()
def log(*s): print(*s, flush=True)
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
def rd(p):
    h = np.fromfile(p, np.int32, 256); n_ = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n_ ** 3).reshape(n_, n_, n_).astype(np.float64)
lamb = 12.2643247 / np.sqrt(a.kv * 1e3 * (1. + a.kv * 1e3 * 0.978466e-6))
def rotm(e1, e2, e3):
    a1, a2, a3 = np.radians(e1), np.radians(e2), np.radians(e3); c1, si1, c2, si2, c3, si3 = np.cos(a1), np.sin(a1), np.cos(a2), np.sin(a2), np.cos(a3), np.sin(a3)
    Z1 = np.zeros((len(e1), 3, 3)); Z1[:, 0, 0] = c1; Z1[:, 0, 1] = -si1; Z1[:, 1, 0] = si1; Z1[:, 1, 1] = c1; Z1[:, 2, 2] = 1
    Y2 = np.zeros((len(e1), 3, 3)); Y2[:, 0, 0] = c2; Y2[:, 0, 2] = si2; Y2[:, 2, 0] = -si2; Y2[:, 2, 2] = c2; Y2[:, 1, 1] = 1
    Z3 = np.zeros((len(e1), 3, 3)); Z3[:, 0, 0] = c3; Z3[:, 0, 1] = -si3; Z3[:, 1, 0] = si3; Z3[:, 1, 1] = c3; Z3[:, 2, 2] = 1
    return Z1 @ Y2 @ Z3
SIMS = a.sim.split(','); full = a.mr == 'full'; n = int(np.frombuffer(open(SIMS[0] + '/all.mrcs', 'rb').read(1024), np.int32)[0]); c = a.crop
fc = np.fft.fftfreq(c) * c / n; fy2, fx2 = np.meshgrid(fc, fc, indexing='ij'); kk0 = np.stack([fx2, fy2, np.zeros_like(fx2)], -1).reshape(-1, 3) * n
kb = min(c // 2 - 1, int(round(n * a.smpd / a.em_lim))); r2 = np.round(np.hypot(kk0[:, 0], kk0[:, 1])).astype(int); cols = np.where((r2 >= 2) & (r2 <= kb))[0]; psh = r2[cols]; ns = c // 2 + 1
kf = (np.hypot(fx2, fy2) / a.smpd).ravel()[cols]; ang = np.arctan2(fy2, fx2).ravel()[cols]; fxp, fyp = fx2.ravel()[cols], fy2.ravel()[cols]
par = (np.mod(np.round(kk0[cols, 0]).astype(int) + np.round(kk0[cols, 1]).astype(int), 2) == 1); k3 = kk0[cols]
VA = [np.c_[np.zeros(len(cols)), np.zeros(len(cols)), k3[:, 1]], np.c_[np.zeros(len(cols)), np.zeros(len(cols)), -k3[:, 0]], np.c_[-k3[:, 1], k3[:, 0], np.zeros(len(cols))]]   # e_a x k
R0 = a.mskdiam / 2 / a.smpd; g3 = np.arange(c) - c // 2; r3 = np.sqrt(g3[:, None, None] ** 2 + g3[None, :, None] ** 2 + g3[None, None, :] ** 2) * n / c
msk3 = np.clip((R0 + 4 - r3) / 8., 0, 1)
def crop_real(v):
    nb = v.shape[0]; sb = np.r_[0:c // 2, nb - c // 2:nb]; Fv = np.fft.fftn(np.fft.ifftshift(v)); return np.fft.fftshift(np.real(np.fft.ifftn(Fv[np.ix_(sb, sb, sb)]))) * msk3
def four(vc): return np.r_[np.fft.fftn(np.fft.ifftshift(vc)).ravel(), 0.].astype(np.complex64)
MF = [l.split()[0] for l in open(a.maps) if l.strip()]; K = len(MF); V, GR = [], []
X3 = [g3[None, None, :] * np.ones((c, c, c)), g3[None, :, None] * np.ones((c, c, c)), g3[:, None, None] * np.ones((c, c, c))]   # q0 <-> array axis 2, q1 <-> 1, q2 <-> 0
for f in MF:
    vc = crop_real(rd(f)); V.append(four(vc)); GR.append([four(vc * (-2j * np.pi * X3[j] / c)) for j in range(3)])
REF = np.mean(V, 0); log('%d path maps cropped to %d^3 with Fourier gradients; band shells 2..%d (%.1f A), %d pixels' % (K, c, kb, n * a.smpd / kb, len(cols)))
NM = sum(1 for _ in open(a.rowflags)); flag = np.loadtxt(a.rowflags).astype(int) != 0
P_all = np.zeros((NM, 8)); EO = np.full(NM, -1)
for r, l in enumerate(open(a.poses)):
    if r < NM and flag[r]: d = kv(l); P_all[r] = [float(d.get(k, 0)) for k in ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast')]; EO[r] = int(float(d.get('eo', -1)))
Sinv = np.diag(1. / np.r_[[np.radians(a.rot_sd) ** 2] * 3, [a.sh_sd ** 2] * 2]); logdetS = -np.linalg.slogdet(Sinv)[1]
PAIRS = [(0, K - 1), (0, K // 2), (K // 2, K - 1)]
halves = [('A', np.where(flag)[0], np.where(flag)[0])] if full else []
if not full:
    for h, tg in enumerate('AB'):
        mr = np.loadtxt(a.mr.split(',')[h]).astype(int) - 1; sel = np.where(flag[mr])[0]; halves.append((tg, sel, mr[sel]))
OUT = {k_: [] for k_ in ('rows', 'half', 'eo', 'l0', 'lprof', 'lmarg', 'l0A', 'l0B', 'ppn', 'dA', 'dB', 'vA', 'vB')}
for h_, (tg, crow, mrow) in enumerate(halves):
    Fc = np.load('%s_%s_c%d.npy' % (a.cache, tg, c), mmap_mode='r'); N = len(crow); log('half %s: %d particles' % (tg, N))
    for b0 in range(0, N, a.block):
        cr = crow[b0:b0 + a.block]; mr_ = mrow[b0:b0 + a.block]; pp = P_all[mr_]; nb_ = len(cr); F = np.asarray(Fc[cr]).reshape(nb_, -1)[:, cols]
        df = 0.5 * ((pp[:, 5] + pp[:, 6])[:, None] + (pp[:, 5] - pp[:, 6])[:, None] * np.cos(2. * (ang[None] - np.radians(pp[:, 7])[:, None]))) * 1e4
        chi = np.pi * lamb * df * kf[None] ** 2 - 0.5 * np.pi * a.cs * 1e7 * lamb ** 3 * kf[None] ** 4; C = (np.sqrt(1. - a.fraca ** 2) * np.sin(chi) + a.fraca * np.cos(chi)).astype(np.float32)
        if a.ctfabs: C = np.abs(C)
        ramp = np.exp(2j * np.pi * (fxp[None] * pp[:, 3, None] + fyp[None] * pp[:, 4, None])).astype(np.complex64)
        Wv = C * np.conj(ramp) * F; Wd = C ** 2; F2 = (np.abs(F) ** 2).astype(np.float32)
        Rm = rotm(pp[:, 0], pp[:, 1], pp[:, 2]); q = np.einsum('qj,nij->nqi', k3, Rm); qi = np.round(q).astype(np.int32); okq = (np.abs(qi) <= c // 2 - 1).all(2); qi = qi % c
        ix = qi[..., 2] * c * c + qi[..., 1] * c + qi[..., 0]; ix[~okq] = c ** 3
        Pr = REF[ix]; bn = np.bincount(np.broadcast_to(psh, Pr.shape).ravel(), np.real(Pr * np.conj(Wv)).ravel(), ns); bd = np.bincount(np.broadcast_to(psh, Pr.shape).ravel(), (Wd * np.abs(Pr) ** 2).ravel(), ns)
        beta = np.clip(bn / np.maximum(bd, 1e-30), 0., None).astype(np.float32); bp = beta[psh][None]
        res = F2 - 2. * np.real(bp * Pr * np.conj(Wv)) + Wd * (bp * np.abs(Pr)) ** 2
        wsh = 1. / np.maximum(np.bincount(np.broadcast_to(psh, res.shape).ravel(), res.ravel(), ns) / np.maximum(np.bincount(psh, minlength=ns) * nb_, 1), 1e-30)
        wi0 = wsh[psh][None].astype(np.float32); ppn = (res * wi0).mean(1); wi = (wi0 / np.maximum(ppn, 1e-6)[:, None]).astype(np.float32)
        DQ = [np.einsum('pj,nij->npi', va, Rm).astype(np.float32) for va in VA]   # rotation of q per small particle-frame rotation
        l0 = np.zeros((nb_, K), np.float32); lp_ = np.zeros_like(l0); lm = np.zeros_like(l0); lA = np.zeros_like(l0); lB = np.zeros_like(l0); PK = {}
        for k in range(K):
            P0 = bp * V[k][ix]; Gq = [GR[k][j][ix] for j in range(3)]
            J = [bp * (Gq[0] * DQ[ax][..., 0] + Gq[1] * DQ[ax][..., 1] + Gq[2] * DQ[ax][..., 2]) for ax in range(3)] + [2j * np.pi * fxp[None] * P0, 2j * np.pi * fyp[None] * P0]
            tt = 0.5 * wi * (2. * np.real(np.conj(P0) * Wv) - Wd * np.abs(P0) ** 2); l0[:, k] = tt.sum(1); lA[:, k] = tt[:, ~par].sum(1); lB[:, k] = tt[:, par].sum(1)
            r = Wv - Wd * P0; g = np.stack([(wi * np.real(np.conj(Jx) * r)).sum(1) for Jx in J], 1)
            H = np.zeros((nb_, 5, 5))
            for u in range(5):
                for v in range(u, 5): H[:, u, v] = H[:, v, u] = (wi * Wd * np.real(np.conj(J[u]) * J[v])).sum(1)
            M_ = H + Sinv[None]; sol = np.linalg.solve(M_, g[..., None])[..., 0]; lp_[:, k] = l0[:, k] + 0.5 * (g * sol).sum(1)
            lm[:, k] = lp_[:, k] - 0.5 * (np.linalg.slogdet(M_)[1] + logdetS)
            if any(k in pr for pr in PAIRS): PK[k] = P0
        dA = np.stack([lA[:, j] - lA[:, k] for j, k in PAIRS], 1); dB = np.stack([lB[:, j] - lB[:, k] for j, k in PAIRS], 1)
        vA = np.stack([(wi * Wd * np.abs(PK[j] - PK[k]) ** 2)[:, ~par].sum(1) for j, k in PAIRS], 1); vB = np.stack([(wi * Wd * np.abs(PK[j] - PK[k]) ** 2)[:, par].sum(1) for j, k in PAIRS], 1)
        for k_, v_ in (('rows', mr_), ('half', np.full(nb_, h_)), ('eo', EO[mr_]), ('l0', l0), ('lprof', lp_), ('lmarg', lm), ('l0A', lA), ('l0B', lB), ('ppn', ppn), ('dA', dA), ('dB', dB), ('vA', vA), ('vB', vB)):
            OUT[k_].append(v_)
        log('  block %d-%d done' % (b0, b0 + nb_))
OUT = {k_: np.concatenate(v_) for k_, v_ in OUT.items()}
ts = np.median((OUT['dA'] - OUT['dB']) ** 2 / np.maximum(OUT['vA'] + OUT['vB'], 1e-30), 0); ts_m = np.mean((OUT['dA'] - OUT['dB']) ** 2, 0) / np.mean(OUT['vA'] + OUT['vB'], 0)
log('pixel-split temperature (variance ratio of A-B contrasts to the model) per pair: mean-ratio %s | median %s' % (' '.join('%.3f' % x for x in ts_m), ' '.join('%.3f' % x for x in ts)))
np.savez(a.out, maps=np.array(MF), tau_split=ts_m, **OUT); log('wrote %s (%d particles x %d maps); DONE' % (a.out, len(OUT['rows']), K))
