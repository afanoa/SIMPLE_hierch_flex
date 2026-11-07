#!/usr/bin/env python3
"""Score-space embedding of a finished tree: every particle scored against every leaf map at its own leaf's refined pose, each particle
against the half-map it did not contribute to (its eo flag: eo 0 -> the odd map, eo 1 -> the even map), Gaussian log-likelihood ratio
vs the zero model over shells 2..EM_LIM with a per-shell reference scale fitted on the own-leaf references and per-shell noise from the
own-leaf residuals (as minor_em.py). The score vector minus its particle mean -> leading principal components -> UMAP of a subsample.
Writes OUT/scores_{A,B}.npy (float16, half rows x leaves, 0 outside), OUT/embed.npz (xy, leaf, row, pcs, evr, own-leaf rank).
usage: score_embed.py --sim A,B --cache C --poses PA,PB --leaves LIST (name even odd count) --labels master_labels --out DIR [--em_lim 15 --nsub 80000]"""
import os, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--sim', required=True); ap.add_argument('--cache', required=True); ap.add_argument('--poses', required=True)
ap.add_argument('--leaves', required=True); ap.add_argument('--labels', required=True); ap.add_argument('--out', required=True); ap.add_argument('--em_lim', type=float, default=15.)
ap.add_argument('--nsub', type=int, default=80000); ap.add_argument('--npc', type=int, default=10); ap.add_argument('--crop', type=int, default=48); ap.add_argument('--mskdiam', type=float, default=230.)
ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--kv', type=float, default=300.); ap.add_argument('--cs', type=float, default=2.7); ap.add_argument('--fraca', type=float, default=0.1)
ap.add_argument('--ctfabs', type=int, default=1); a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
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
SIMS = a.sim.split(','); n = int(np.frombuffer(open(SIMS[0] + '/all.mrcs', 'rb').read(1024), np.int32)[0]); c = a.crop
x = np.arange(n) - n // 2; R0 = a.mskdiam / 2 / a.smpd; sl = np.r_[0:c // 2, n - c // 2:n]
fc = np.fft.fftfreq(c) * c / n; fy2, fx2 = np.meshgrid(fc, fc, indexing='ij'); kk0 = np.stack([fx2, fy2, np.zeros_like(fx2)], -1).reshape(-1, 3) * n
f3 = np.fft.fftfreq(c) * c / n; kx, ky, kz = np.meshgrid(f3, f3, f3, indexing='ij'); sh = np.minimum(np.round(np.sqrt(kx ** 2 + ky ** 2 + kz ** 2) * n).astype(int), c // 2).ravel()
ns = c // 2 + 1; kb = min(c // 2 - 1, int(round(n * a.smpd / a.em_lim))); r2 = np.round(np.hypot(kk0[:, 0], kk0[:, 1])).astype(int); cols = np.where((r2 >= 2) & (r2 <= kb))[0]
kf = (np.hypot(fx2, fy2) / a.smpd).ravel()[cols]; ang = np.arctan2(fy2, fx2).ravel()[cols]
LAB = np.loadtxt(a.labels, dtype=int); L = [l.split() for l in open(a.leaves) if l.strip() and not l.startswith('#')]; K = len(L)
def crop_four(v):   # box-n real-space map -> the c^3 FFT-layout Fourier map of this script's frame (identity frame, as minor_em2)
    nb = v.shape[0]; sb = np.r_[0:c // 2, nb - c // 2:nb]   # the map's own box (same physical size as the particle box, e.g. refine3D's 88 at 2.66 A)
    F = np.fft.fftn(np.fft.ifftshift(v)); vc = np.fft.fftshift(np.real(np.fft.ifftn(F[np.ix_(sb, sb, sb)]))); return np.r_[np.fft.fftn(np.fft.ifftshift(vc)).ravel(), 0.].astype(np.complex64)
REF = [(crop_four(rd(l[1])), crop_four(rd(l[2]))) for l in L]; log('%d leaf maps cropped to %d^3 (band to %.1f A, %d pixels per image)' % (K, c, n * a.smpd / kb, len(cols)))
out_rows = []
for h, tag in enumerate('AB'):
    mr = np.loadtxt(SIMS[h] + '/meta_rows.txt').astype(int) - 1; lab = LAB[mr]; rows = np.where(lab > 0)[0]
    P = []; EO = []
    for l in open(a.poses.split(',')[h]):   # one line per half row (a blank line keeps its slot)
        d = kv(l); P.append([float(d.get(k, 0)) for k in ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast')]); EO.append(int(float(d.get('eo', 0))))
    assert len(P) == len(mr), 'poses %d lines vs %d half rows' % (len(P), len(mr))
    P = np.array(P)[rows]; EO = np.array(EO)[rows]; Fc = np.load('%s_%s_c%d.npy' % (a.cache, tag, c), mmap_mode='r'); N = len(rows)
    S = np.zeros((N, K), np.float32); own = lab[rows] - 1
    for b in range(0, N, 20000):
        rr = rows[b:b + 20000]; pp = P[b:b + 20000]; F = np.asarray(Fc[rr]).reshape(len(rr), -1)[:, cols]
        df = 0.5 * ((pp[:, 5] + pp[:, 6])[:, None] + (pp[:, 5] - pp[:, 6])[:, None] * np.cos(2. * (ang[None] - np.radians(pp[:, 7])[:, None]))) * 1e4
        chi = np.pi * lam * df * kf[None] ** 2 - 0.5 * np.pi * a.cs * 1e7 * lam ** 3 * kf[None] ** 4; C = (np.sqrt(1. - a.fraca ** 2) * np.sin(chi) + a.fraca * np.cos(chi)).astype(np.float32)
        if a.ctfabs: C = np.abs(C)
        ramp = np.exp(2j * np.pi * (fx2.ravel()[cols][None] * pp[:, 3, None] + fy2.ravel()[cols][None] * pp[:, 4, None])).astype(np.complex64)
        Wv = C * np.conj(ramp) * F; Wd = C ** 2; F2 = (np.abs(F) ** 2).astype(np.float32)
        Rm = rotm(pp[:, 0], pp[:, 1], pp[:, 2]); q = np.round(np.einsum('qj,nij->nqi', kk0[cols], Rm)).astype(np.int32)
        ok = (np.abs(q) <= c // 2 - 1).all(2); q = q % c; ix = q[..., 2] * c * c + q[..., 1] * c + q[..., 0]; ix[~ok] = c ** 3; shq = np.minimum(np.r_[sh, 0][ix], ns - 1)
        eo = EO[b:b + 20000]; ob = own[b:b + 20000]
        # per-shell reference scale (regression of the data on the own-leaf reference) and noise, per block (blocks of 20k are large enough)
        Pown = np.empty_like(Wv)
        for g in range(K):
            m = ob == g
            if not m.any(): continue
            for e_ in (0, 1):
                mm = m & (eo == e_)
                if mm.any(): Pown[mm] = REF[g][1 - e_][ix[mm]]   # eo 0 -> the odd map (index 1), eo 1 -> the even map
        bn = np.bincount(shq.ravel(), np.real(Pown * np.conj(Wv)).ravel(), ns); bd = np.bincount(shq.ravel(), (Wd * np.abs(Pown) ** 2).ravel(), ns)
        beta = np.clip(bn / np.maximum(bd, 1e-30), 0., 4.).astype(np.float32); bpx = beta[shq]
        res = F2 - 2. * np.real(bpx * Pown * np.conj(Wv)) + Wd * (bpx * np.abs(Pown)) ** 2
        wsh = (1. / np.maximum(np.bincount(shq.ravel(), res.ravel(), ns) / np.maximum(np.bincount(shq.ravel(), minlength=ns), 1), 1e-12)).astype(np.float32); wi = wsh[shq]
        for g in range(K):
            Pg = np.where((eo == 0)[:, None], REF[g][1][ix], REF[g][0][ix]) * bpx
            S[b:b + len(rr), g] = 0.5 * ((2. * np.real(Pg * np.conj(Wv)) - Wd * np.abs(Pg) ** 2) * wi).sum(1)
    full = np.zeros((len(mr), K), np.float16); full[rows] = S; np.save('%s/scores_%s.npy' % (a.out, tag), full)
    rank = (S > S[np.arange(N), own][:, None]).sum(1); log('half %s: %d particles scored against %d leaves; own leaf ranks first for %.1f %% (top 3 %.1f %%)' % (tag, N, K, 100 * (rank == 0).mean(), 100 * (rank < 3).mean()))
    out_rows.append((mr[rows], own, S, rank))
mrow = np.r_[out_rows[0][0], out_rows[1][0]]; own = np.r_[out_rows[0][1], out_rows[1][1]]; S = np.r_[out_rows[0][2], out_rows[1][2]]; rank = np.r_[out_rows[0][3], out_rows[1][3]]
Sc = S - S.mean(1, keepdims=True); Sc = (Sc - Sc.mean(0)) / np.maximum(Sc.std(0), 1e-12)
w, V = np.linalg.eigh(np.cov(Sc.T)); o = np.argsort(-w); w, V = w[o], V[:, o]; PC = Sc @ V[:, :a.npc]
rng = np.random.default_rng(0); sub = np.sort(rng.choice(len(S), min(a.nsub, len(S)), replace=False))
import umap
Y = umap.UMAP(n_neighbors=30, min_dist=0.3, random_state=0).fit_transform(PC[sub])
np.savez_compressed('%s/embed.npz' % a.out, xy=Y.astype(np.float32), leaf=own[sub] + 1, row=mrow[sub], pcs=PC[sub, :3].astype(np.float32), evr=(w / w.sum())[:12], rank=rank[sub], names=np.array([l[0] for l in L]))
log('score space: %d particles x %d leaves; PC shares %s; UMAP of %d' % (len(S), K, ' '.join('%.2f' % v for v in (w / w.sum())[:8]), len(sub)))
