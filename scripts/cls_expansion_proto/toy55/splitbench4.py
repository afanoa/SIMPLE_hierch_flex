#!/usr/bin/env python3
"""Fast in-memory split test of one tree node on a pair of independent half-sets (splitbench3 + --xem: the re-assignment
after the first decision is a CROSS-HALF EM, each half scored against the other half's filtered group maps).
Every half's particle images are Fourier-transformed once (cropped, CTF and slice positions from the refined poses);
a group's map is direct Fourier insertion sum(CTF conj(ramp) F) / (sum CTF^2 + eps) with ONE fixed per-shell regulariser
for every group, from the same number of particles (equal N), so neither size nor reconstruction software can make two
groups differ. Groups of the two halves are paired on the coarse shells (below the tested band), every pair difference
is tested by the best 3-shell window of its cross-half FSC (bandlo .. crop Nyquist), the threshold is a high quantile of
the same statistic over NSHUF label shuffles, then pairwise merging and archetype (blend) pruning give K*; the blend-residual
threshold is the same quantile of the blend-residual statistic over the shuffles (its own null).
usage: splitbench.py --sim SIM_A,SIM_B --poses P_A,P_B --cand C_A,C_B --out DIR [--gt GT_A,GT_B] [--truth NPZ]"""
import sys, os, argparse, numpy as np
from scipy.optimize import linear_sum_assignment, nnls
ap = argparse.ArgumentParser()
ap.add_argument('--sim', required=True); ap.add_argument('--poses', required=True); ap.add_argument('--cand', required=True); ap.add_argument('--out', required=True)
ap.add_argument('--gt', default=''); ap.add_argument('--truth', default=''); ap.add_argument('--crop', type=int, default=48); ap.add_argument('--nshuf', type=int, default=20)
ap.add_argument('--q', type=float, default=1.0); ap.add_argument('--bandlo', type=float, default=40.); ap.add_argument('--mskdiam', type=float, default=230.)
ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--kv', type=float, default=300.); ap.add_argument('--cs', type=float, default=2.7); ap.add_argument('--fraca', type=float, default=0.1)
ap.add_argument('--conv', default='1,1,1,1,1,1', help='transpose, e1 sign, y flip, e3 sign, shift sign, ctf sign (arm_llr.py calibration of these sets)')
ap.add_argument('--xem', type=int, default=0, help='> 0: cross-half EM iterations replace the per-half soft EM'); ap.add_argument('--em_lim', type=float, default=20., help='cross-half EM band limit (A)'); ap.add_argument('--em_min', type=int, default=100, help='a group below this many particles in either half is dropped'); ap.add_argument('--eps', type=float, default=0.05); ap.add_argument('--seed', type=int, default=0); ap.add_argument('--arch', type=int, default=1); ap.add_argument('--junkcorr', type=float, default=0.8, help='states whose coarse-shell map correlation with the node consensus is below this are flagged junk-like'); ap.add_argument('--softem', type=int, default=8, help='soft-EM iterations refining the accepted groups as a mixture before the final test (0 = off)'); ap.add_argument('--archmode', default='shell', help='global: one convex blend | shell: a convex blend per Fourier shell (absorbs the radially filtered imprint of hard classification along a real axis)'); ap.add_argument('--cache', default='')
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); rng = np.random.default_rng(a.seed)
cv = [float(v) for v in a.conv.split(',')]; convT, s1, yf, s3, ssh, sctf = bool(cv[0]), cv[1], cv[2], cv[3], cv[4], cv[5]
def log(*s): print(*s, flush=True)
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
# ---------------- geometry on the cropped grid (shell k = n*smpd/k A, as arm_llr.py / merge_cross_tree.py)
def setup(n):
    c = a.crop; x = np.arange(n) - n // 2; rr2 = np.hypot(x[:, None], x[None, :]); R = a.mskdiam / 2 / a.smpd
    msk2 = np.clip((R + 4 - rr2) / 8., 0, 1).astype(np.float32); sl = np.r_[0:c // 2, n - c // 2:n]
    fc = np.fft.fftfreq(c) * c / n; fy2, fx2 = np.meshgrid(fc, fc, indexing='ij'); kk0 = np.stack([fx2, fy2, np.zeros_like(fx2)], -1).reshape(-1, 3) * n
    f3 = np.fft.fftfreq(c) * c / n; kx, ky, kz = np.meshgrid(f3, f3, f3, indexing='ij'); sh3 = np.minimum(np.round(np.sqrt(kx**2 + ky**2 + kz**2) * n).astype(int), c // 2).ravel()
    g = np.arange(c) - c // 2; r3 = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2) * n / c; msk3 = np.clip((R + 4 - r3) / 8., 0, 1).astype(np.float32)
    return dict(n=n, c=c, msk2=msk2, sl=sl, fx2=fx2, fy2=fy2, kk0=kk0, sh3=sh3, msk3=msk3)
lam = 12.2643247 / np.sqrt(a.kv * 1e3 * (1. + a.kv * 1e3 * 0.978466e-6))
def rotm(e1, e2, e3):   # arm_llr.py rot(), vectorised
    a1, a2, a3 = np.radians(e1), np.radians(e2), np.radians(e3); c1, si1, c2, si2, c3, si3 = np.cos(a1), np.sin(a1), np.cos(a2), np.sin(a2), np.cos(a3), np.sin(a3)
    Z1 = np.zeros((len(e1), 3, 3)); Z1[:, 0, 0] = c1; Z1[:, 0, 1] = -si1; Z1[:, 1, 0] = si1; Z1[:, 1, 1] = c1; Z1[:, 2, 2] = 1
    Y2 = np.zeros((len(e1), 3, 3)); Y2[:, 0, 0] = c2; Y2[:, 0, 2] = si2; Y2[:, 2, 0] = -si2; Y2[:, 2, 2] = c2; Y2[:, 1, 1] = 1
    Z3 = np.zeros((len(e1), 3, 3)); Z3[:, 0, 0] = c3; Z3[:, 0, 1] = -si3; Z3[:, 1, 0] = si3; Z3[:, 1, 1] = c3; Z3[:, 2, 2] = 1
    return Z1 @ Y2 @ Z3
# ---------------- one half: images -> weighted Fourier coefficients and 3-D insertion indices for the node's particles
def load_half(simdir, posef, candf, tag):
    lab = np.loadtxt(candf).astype(int); rows = np.where(lab > 0)[0]
    P = np.array([[float(kv(l).get(k, 0)) for k in ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast')] for l in open(posef) if l.strip()])
    h = np.frombuffer(open(simdir + '/all.mrcs', 'rb').read(1024), np.int32); nx, ny, nz, nsym = [int(v) for v in (h[0], h[1], h[2], h[23])]
    I = np.memmap(simdir + '/all.mrcs', np.float32, 'r', offset=1024 + nsym, shape=(nz, ny, nx)); G = setup(nx); c = G['c']
    cache = '%s_%s_c%d.npy' % (a.cache, tag, c) if a.cache else ''   # every particle of the half, transformed once
    if cache and os.path.isfile(cache): F = np.load(cache, mmap_mode='r')[rows]
    else:
        Fall = np.zeros((nz, c, c), np.complex64)
        for b in range(0, nz, 4000):
            blk = np.array(I[b:b + 4000], np.float32) * G['msk2']; Fall[b:b + len(blk)] = np.fft.fft2(np.fft.ifftshift(blk, axes=(1, 2)))[:, G['sl']][:, :, G['sl']]
        if cache: np.save(cache, Fall)
        F = Fall[rows]; del Fall
    pp = P[rows]; kf = np.hypot(G['fx2'], G['fy2']) / a.smpd; ang = np.arctan2(G['fy2'], G['fx2'])
    df = 0.5 * ((pp[:, 5] + pp[:, 6])[:, None, None] + (pp[:, 5] - pp[:, 6])[:, None, None] * np.cos(2. * (ang[None] - np.radians(pp[:, 7])[:, None, None]))) * 1e4
    chi = np.pi * lam * df * kf[None] ** 2 - 0.5 * np.pi * a.cs * 1e7 * lam ** 3 * kf[None] ** 4; C = (sctf * (np.sqrt(1. - a.fraca ** 2) * np.sin(chi) + a.fraca * np.cos(chi))).astype(np.float32)
    ramp = np.exp(ssh * 2j * np.pi * (G['fx2'][None] * pp[:, 3, None, None] + G['fy2'][None] * pp[:, 4, None, None])).astype(np.complex64)
    Wv = (C * np.conj(ramp) * F).reshape(len(rows), -1); Wd = (C ** 2).reshape(len(rows), -1); F2 = (np.abs(F) ** 2).reshape(len(rows), -1).astype(np.float32)
    Rm = rotm(s1 * pp[:, 0], pp[:, 1], s3 * pp[:, 2]); kk = G['kk0'] * np.array([1., yf, 1.])
    p3 = np.einsum('qj,nij->nqi', kk, Rm) if convT else np.einsum('qj,nji->nqi', kk, Rm); q = np.round(p3).astype(np.int32)
    ok = (np.abs(q) <= c // 2 - 1).all(2); q = q % c; idx = (q[..., 2] * c * c + q[..., 1] * c + q[..., 0]).astype(np.int64); idx[~ok] = c ** 3
    log('%s: %d node particles (labels %s), crop %d (Nyquist %.1f A)' % (tag, len(rows), np.bincount(lab[rows]).tolist()[1:], c, nx * a.smpd / (c // 2)))
    return dict(lab=lab, rows=rows, Wv=Wv, Wd=Wd, F2=F2, idx=idx, G=G, nrow=len(lab))
def rec(H, sel, eps):   # Fourier map of the node particles 'sel' (indices into H['rows']), masked in real space
    c = H['G']['c']; ii = H['idx'][sel].ravel(); m = c ** 3 + 1
    num = np.bincount(ii, H['Wv'][sel].real.ravel(), m)[:-1] + 1j * np.bincount(ii, H['Wv'][sel].imag.ravel(), m)[:-1]; den = np.bincount(ii, H['Wd'][sel].ravel(), m)[:-1]
    V = num / (den + eps[H['G']['sh3']]); v = np.real(np.fft.ifftn(V.reshape(c, c, c))); v = np.fft.fftshift(v) * H['G']['msk3']
    return np.fft.fftn(np.fft.ifftshift(v)).ravel()
def shell_fsc(X, Y, sh, ns):
    num = np.bincount(sh, np.real(X * np.conj(Y)), ns); dx = np.bincount(sh, np.abs(X) ** 2, ns); dy = np.bincount(sh, np.abs(Y) ** 2, ns); return num / np.maximum(np.sqrt(dx * dy), 1e-30)
# ---------------- groups of one labelling -> equal-N maps of every group in both halves
def group_maps(HA, HB, labA, labB, eps, neq):
    out = []
    for H, lab in ((HA, labA), (HB, labB)):
        M = {}
        for g in np.unique(lab[lab > 0]):
            sel = np.where(lab == g)[0]; sel = rng.choice(sel, neq, replace=False); M[g] = rec(H, sel, eps)
        out.append(M)
    return out
def pair_and_test(MA, MB, cons_A, cons_B, sh, ns, klo, khi):
    ga, gb = sorted(MA), sorted(MB); low = (sh >= 2) & (sh < klo - 1)
    def lowcorr(X, Y): x, y = X[low], Y[low]; return float(np.real(np.vdot(x, y)) / max(np.linalg.norm(x) * np.linalg.norm(y), 1e-30))
    S = np.array([[lowcorr(MA[i] - cons_A, MB[j] - cons_B) for j in gb] for i in ga]); ia, ib = linear_sum_assignment(-S)
    pairs = [(ga[i], gb[j]) for i, j in zip(ia, ib)]
    extra = [{}, {}]   # unpaired groups (unequal counts): label -> the most similar paired label of the same half
    for side, M, cons, mine in ((0, MA, cons_A, [p[0] for p in pairs]), (1, MB, cons_B, [p[1] for p in pairs])):
        for u in [g for g in M if g not in mine]: extra[side][u] = max(mine, key=lambda m: lowcorr(M[u] - cons, M[m] - cons))
    pair_and_test.extra = extra
    def stat(DA, DB): f = shell_fsc(DA, DB, sh, ns); return max(float(f[k:k + 3].mean()) for k in range(klo, max(khi - 2, klo) + 1))
    R = {(p, q): stat(MA[pairs[p][0]] - MA[pairs[q][0]], MB[pairs[p][1]] - MB[pairs[q][1]]) for p in range(len(pairs)) for q in range(p + 1, len(pairs))}
    return pairs, R, stat, float(np.median([S[i, j] for i, j in zip(ia, ib)]))
HA = load_half(*[s.split(',')[0] for s in (a.sim, a.poses, a.cand)], 'A'); HB = load_half(*[s.split(',')[1] for s in (a.sim, a.poses, a.cand)], 'B')
G = HA['G']; c = G['c']; n = G['n']; sh = G['sh3']; ns = c // 2 + 1; klo = max(int(round(n * a.smpd / a.bandlo)), 3); khi = c // 2 - 1
labA, labB = HA['lab'][HA['rows']], HB['lab'][HB['rows']]
# one fixed per-shell regulariser: eps x the median per-shell weight of an N_eq reconstruction (same for every group map)
neq = int(min(np.bincount(labA)[1:][np.bincount(labA)[1:] > 0].min(), np.bincount(labB)[1:][np.bincount(labB)[1:] > 0].min()))
def shell_eps(H):
    c_ = H['G']['c']; ii = H['idx'].ravel(); den = np.bincount(ii, H['Wd'].ravel(), c_ ** 3 + 1)[:-1] * (neq / len(H['rows']))
    return np.array([np.median(den[(sh == k) & (den > 0)]) if ((sh == k) & (den > 0)).any() else 1. for k in range(ns)]) * a.eps
eps = 0.5 * (shell_eps(HA) + shell_eps(HB)); consA, consB = rec(HA, np.arange(len(HA['rows'])), eps * len(HA['rows']) / neq), rec(HB, np.arange(len(HB['rows'])), eps * len(HB['rows']) / neq)
log('equal N %d per group; band %d..%d (%.1f..%.1f A); %d shuffles, threshold quantile %.3f' % (neq, klo, khi, n * a.smpd / klo, n * a.smpd / khi, a.nshuf, a.q))
def _cblend(X, y):   # convex least squares: w >= 0, sum w = 1
    lamb = 10. * max(np.abs(X).max(), 1e-12) * np.sqrt(len(y)); return nnls(np.vstack([X, lamb * np.ones((1, X.shape[1]))]), np.r_[y, lamb])[0]
def blend_resid(Gd, p, others):
    if a.archmode == 'shell':   # per-shell convex weights
        r = Gd[p].copy(); wsum = np.zeros(len(others))
        for k in range(ns):
            m = sh == k
            if not m.any(): continue
            X = np.stack([np.concatenate([Gd[o][m].real, Gd[o][m].imag]) for o in others], 1); y = np.concatenate([Gd[p][m].real, Gd[p][m].imag]); w = _cblend(X, y)
            r[m] = Gd[p][m] - sum(wi * Gd[o][m] for wi, o in zip(w, others)); wsum += w
        return r, wsum
    X = np.stack([np.concatenate([Gd[o].real, Gd[o].imag]) for o in others], 1); y = np.concatenate([Gd[p].real, Gd[p].imag]); w = _cblend(X, y)
    return Gd[p] - sum(wi * Gd[o] for wi, o in zip(w, others)), w
# ---------------- shuffle null (pair differences and blend residuals)
NUL = []; NULB = []
for s in range(a.nshuf):
    pa, pb = labA.copy(), labB.copy(); rng.shuffle(pa); rng.shuffle(pb); MA, MB = group_maps(HA, HB, pa, pb, eps, neq)
    prs, R, st_, _ = pair_and_test(MA, MB, consA, consB, sh, ns, klo, khi); NUL += list(R.values())
    if a.arch and len(prs) > 2:
        gA = {p: MA[prs[p][0]] for p in range(len(prs))}; gB = {p: MB[prs[p][1]] for p in range(len(prs))}
        for p in gA: oth = [o for o in gA if o != p]; NULB.append(st_(blend_resid(gA, p, oth)[0], blend_resid(gB, p, oth)[0]))
NULB = np.array(NULB) if NULB else np.zeros(1); THRB = float(np.quantile(NULB, a.q)); NUL = np.array(NUL); THR = float(np.quantile(NUL, a.q)); log('blend-residual null (%d): median %.3f, 99 %% %.3f, max %.3f -> blend threshold %.3f' % (len(NULB), np.median(NULB), np.quantile(NULB, .99), NULB.max(), THRB)); log('shuffle null (%d pair statistics): median %.3f, 95 %% %.3f, 99 %% %.3f, max %.3f -> threshold %.3f' % (len(NUL), np.median(NUL), np.quantile(NUL, .95), np.quantile(NUL, .99), NUL.max(), THR))
def decide(labA, labB, tag):
  global pairs
  neq_ = int(min(np.bincount(labA)[1:][np.bincount(labA)[1:] > 0].min(), np.bincount(labB)[1:][np.bincount(labB)[1:] > 0].min()))
  MA, MB = group_maps(HA, HB, labA, labB, eps * neq_ / neq, neq_); pairs, R, stat, pc = pair_and_test(MA, MB, consA, consB, sh, ns, klo, khi); log('[%s] equal N %d' % (tag, neq_))
  log('paired %d groups (coarse-shell correlation median %.3f): %s' % (len(pairs), pc, ' '.join('%d<->%d' % p for p in pairs)))
  log('pair statistics at the start (sorted): %s' % ' '.join('%.2f' % v for v in sorted(R.values())))
  GA = {p: MA[pairs[p][0]] for p in range(len(pairs))}; GB = {p: MB[pairs[p][1]] for p in range(len(pairs))}; groups = {p: [p] for p in range(len(pairs))}
  def gstat(p, q): return stat(GA[p] - GA[q], GB[p] - GB[q])
  R = {k: gstat(*k) for k in R}; merged = []
  while R and min(R.values()) < THR:
      (p, q), v = min(R.items(), key=lambda kv: kv[1]); merged.append(v); np_, nq = len(groups[p]), len(groups[q])
      for D in (GA, GB): D[p] = (np_ * D[p] + nq * D[q]) / (np_ + nq); del D[q]
      groups[p] += groups.pop(q); R = {k: w for k, w in R.items() if q not in k}; R.update({tuple(sorted((p, o))): gstat(*sorted((p, o))) for o in groups if o != p})
  log('%d merges at %s -> %d groups' % (len(merged), ' '.join('%.2f' % v for v in merged), len(groups)))
  while a.arch and len(groups) > 2:
      best = None
      for p in groups:
          others = [o for o in groups if o != p]; rA, wA = blend_resid(GA, p, others); rB, wB = blend_resid(GB, p, others); v = stat(rA, rB)
          if best is None or v < best[1]: best = (p, v, others[int(np.argmax(wA + wB))])
      log('  archetype test: least-novel group %d blend residual %.2f (blend threshold %.3f)' % (best[0], best[1], THRB))
      if best[1] >= THRB: break
      groups[best[2]] += groups.pop(best[0]); GA.pop(best[0]); GB.pop(best[0])
  K = len(groups); log('K* = %d' % K)
  return groups, pairs, K
groups, pairs, K = decide(labA, labB, 'board groups')
def soft_em(H, lab, iters):   # mixture refinement at fixed poses: posteriors from the Gaussian likelihood of the component projections
    c_ = H['G']['c']; ks = sorted(set(lab[lab > 0])); n_ = len(lab); band_ = (sh[H['idx'] % (c_ ** 3)] >= 2) & (H['idx'] < c_ ** 3)
    W = np.zeros((n_, len(ks)), np.float32)
    for j, k in enumerate(ks): W[lab == k, j] = 1.
    for it in range(iters + 1):
        V = []
        for j in range(len(ks)):
            ii = H['idx'].ravel(); w = np.repeat(W[:, j], H['idx'].shape[1]); m = c_ ** 3 + 1
            num = np.bincount(ii, (H['Wv'].real.ravel() * w), m)[:-1] + 1j * np.bincount(ii, (H['Wv'].imag.ravel() * w), m)[:-1]; den = np.bincount(ii, H['Wd'].ravel() * w, m)[:-1]
            V.append(np.r_[num / (den + eps[sh] * max(W[:, j].sum(), 1.) / neq), 0.].astype(np.complex64))
        if it == iters: break
        S = np.zeros((n_, len(ks)), np.float32); R2 = np.zeros(ns); C2 = np.zeros(ns); shq = np.minimum(sh[H['idx'] % (c_ ** 3)], ns - 1)
        for b in range(0, n_, 4000):
            sl_ = slice(b, b + 4000); ix = H['idx'][sl_]; bm = band_[sl_]
            for j in range(len(ks)):
                P = V[j][ix]; t = 2. * np.real(P * np.conj(H['Wv'][sl_])) - H['Wd'][sl_] * np.abs(P) ** 2
                S[sl_, j] = np.where(bm, t, 0.).sum(1) if it == 0 else np.where(bm, t / sig2[shq[sl_]], 0.).sum(1)
                if it == 0: own = W[sl_, j] > 0; res = (H['F2'][sl_] - t); R2 += np.bincount(shq[sl_][own].ravel(), (res[own] * bm[own]).ravel(), ns); C2 += np.bincount(shq[sl_][own].ravel(), bm[own].ravel().astype(float), ns)
        if it == 0: sig2 = np.maximum(R2 / np.maximum(C2, 1), 1e-12); continue
        S = 0.5 * S + np.log(np.maximum(W.mean(0), 1e-6))[None]; S -= S.max(1, keepdims=True); W = np.exp(S); W /= W.sum(1, keepdims=True)
    return np.array(ks)[W.argmax(1)], W
def nmi_(tr, lab):
    m = tr > 0; t = np.unique(tr[m], return_inverse=True)[1]; g = np.unique(lab[m], return_inverse=True)[1]; J = np.zeros((t.max() + 1, g.max() + 1)); np.add.at(J, (t, g), 1); J /= J.sum()
    pt, pg = J.sum(1), J.sum(0); nz = J > 0; return (J[nz] * np.log(J[nz] / np.outer(pt, pg)[nz])).sum() / max(-(pt[pt > 0] * np.log(pt[pt > 0])).sum(), 1e-12)
def cross_em(stA, stB, K, iters):   # soft EM at fixed poses, each half scored against the OTHER half's maps: consensus + group deviation x cross-half FSC per shell
    kb = c // 2 - 1 if a.em_lim <= 0 else min(c // 2 - 1, int(round(n * a.smpd / a.em_lim))); m3 = c ** 3 + 1
    r2 = np.round(np.hypot(G['kk0'][:, 0], G['kk0'][:, 1])).astype(int); cols = np.where((r2 >= 2) & (r2 <= kb))[0]; P = {}
    for h, (H, st) in enumerate(((HA, stA), (HB, stB))):
        ix = np.ascontiguousarray(H['idx'][:, cols]); W = np.zeros((len(st), K), np.float32); W[st > 0, st[st > 0] - 1] = 1.
        P[h] = dict(Wv=np.ascontiguousarray(H['Wv'][:, cols]), Wd=np.ascontiguousarray(H['Wd'][:, cols]), F2=np.ascontiguousarray(H['F2'][:, cols]), ix=ix, shq=np.minimum(sh[ix % (c ** 3)], ns - 1), W=W)
        P[h]['conf'] = np.loadtxt(a.gt.split(',')[h], comments='#')[:, 1].astype(int)[H['rows']] if a.gt else None
    epsv = np.r_[eps[sh], eps[-1]]; alive = np.ones(K, bool)
    log('cross-half EM: %d groups, %d iterations, shells 2..%d (%.1f A), %d pixels per particle' % (K, iters, kb, n * a.smpd / kb, len(cols)))
    def maps(h):
        d = P[h]; ii = d['ix'].ravel(); cons = (np.bincount(ii, d['Wv'].real.ravel(), m3) + 1j * np.bincount(ii, d['Wv'].imag.ravel(), m3)) / (np.bincount(ii, d['Wd'].ravel(), m3) + epsv * len(d['W']) / neq)
        out = []
        for g in range(K):
            w = d['W'][:, g:g + 1]; num = np.bincount(ii, (d['Wv'].real * w).ravel(), m3) + 1j * np.bincount(ii, (d['Wv'].imag * w).ravel(), m3)
            out.append(num / (np.bincount(ii, (d['Wd'] * w).ravel(), m3) + epsv * max(float(w.sum()), 1.) / neq))
        return cons, out
    def estep(h, R):
        d = P[h]; nn = len(d['W']); mp = d['W'].argmax(1); R2 = np.zeros(ns); C2 = np.zeros(ns)
        for g in range(K):
            rows = np.where(mp == g)[0]
            for b in range(0, len(rows), 8000):
                r = rows[b:b + 8000]; Pg = R[g][d['ix'][r]]; res = d['F2'][r] - 2. * np.real(Pg * np.conj(d['Wv'][r])) + d['Wd'][r] * np.abs(Pg) ** 2
                R2 += np.bincount(d['shq'][r].ravel(), res.ravel(), ns); C2 += np.bincount(d['shq'][r].ravel(), minlength=ns)
        winv = (1. / np.maximum(R2 / np.maximum(C2, 1), 1e-12)).astype(np.float32)[d['shq']]; S = np.full((nn, K), -1e30, np.float32)
        for g in np.where(alive)[0]:
            for b in range(0, nn, 8000):
                s_ = slice(b, b + 8000); Pg = R[g][d['ix'][s_]]; S[s_, g] = 0.5 * ((2. * np.real(Pg * np.conj(d['Wv'][s_])) - d['Wd'][s_] * np.abs(Pg) ** 2) * winv[s_]).sum(1)
        S += np.log(np.maximum(d['W'].mean(0), 1e-6))[None]; S -= S.max(1, keepdims=True); W = np.exp(S); W[:, ~alive] = 0.; return W / W.sum(1, keepdims=True)
    for it in range(iters):
        (cA, MA_), (cB, MB_) = maps(0), maps(1); RA_, RB_, fl = [], [], []
        for g in range(K):
            DA, DB = MA_[g][:-1] - cA[:-1], MB_[g][:-1] - cB[:-1]
            f = np.bincount(sh, np.real(DA * np.conj(DB)), ns) / np.maximum(np.sqrt(np.bincount(sh, np.abs(DA) ** 2, ns) * np.bincount(sh, np.abs(DB) ** 2, ns)), 1e-30)
            f = np.clip(np.convolve(np.pad(f, 1, mode='edge'), np.ones(3) / 3., mode='valid'), 0., 1.); fl.append(f)
            RA_.append(np.r_[cA[:-1] + f[sh] * DA, 0.].astype(np.complex64)); RB_.append(np.r_[cB[:-1] + f[sh] * DB, 0.].astype(np.complex64))
        WA, WB = estep(0, RB_), estep(1, RA_)   # half A against half B's references and vice versa
        chA, chB = np.mean(WA.argmax(1) != P[0]['W'].argmax(1)), np.mean(WB.argmax(1) != P[1]['W'].argmax(1)); P[0]['W'], P[1]['W'] = WA, WB
        cA_, cB_ = np.bincount(WA.argmax(1), minlength=K), np.bincount(WB.argmax(1), minlength=K)   # MAP counts: the re-test's equal N comes from them
        for g in np.where(alive)[0]:
            if min(cA_[g], cB_[g]) < a.em_min: alive[g] = False; log('  group %d drained (%d / %d MAP particles) and dropped' % (g + 1, cA_[g], cB_[g]))
        s = '  iteration %d: MAP changed A %.3f B %.3f; weights A %s; deviation FSC per group (shell 4/8/12) %s' % (it + 1, chA, chB, ' '.join('%.3f' % v for v in WA.mean(0)), ' '.join('%.2f/%.2f/%.2f' % (f[4], f[min(8, ns - 1)], f[min(12, ns - 1)]) for f, al in zip(fl, alive) if al))
        if a.gt: s += '; NMI(conformer; group) A %.3f B %.3f' % (nmi_(P[0]['conf'], WA.argmax(1)), nmi_(P[1]['conf'], WB.argmax(1)))
        log(s)
    return P[0]['W'].argmax(1) + 1, P[1]['W'].argmax(1) + 1
if (a.softem or a.xem) and K >= 2:
    starts = []
    for H, lab, k in ((HA, labA, 0), (HB, labB, 1)):   # start from each half's merged groups (one numbering across the halves)
        start = np.zeros_like(lab)
        for gi, (p, mem) in enumerate(groups.items()):
            for m in mem: start[lab == pairs[m][k]] = gi + 1
        starts.append(start)
    if a.gt:
        for k in (0, 1): log('half %s start: NMI(conformer; group) %.3f' % ('AB'[k], nmi_(np.loadtxt(a.gt.split(',')[k], comments='#')[:, 1].astype(int)[(HA, HB)[k]['rows']], starts[k])))
    if a.xem:
        labA_em, labB_em = cross_em(starts[0], starts[1], len(groups), a.xem)
        for k, nl in ((0, labA_em), (1, labB_em)): log('cross-half EM half %s: MAP changed vs start %.3f' % ('AB'[k], np.mean(nl != starts[k])))
    else:
        for H, k in ((HA, 0), (HB, 1)):   # per half, independently
            newlab, W = soft_em(H, starts[k], a.softem); log('soft EM half %s: weights %s, MAP changed vs start %.3f' % ('AB'[k], ' '.join('%.3f' % v for v in W.mean(0)), np.mean(newlab != starts[k])))
            if k == 0: labA_em = newlab
            else: labB_em = newlab
    labA, labB = labA_em, labB_em; groups, pairs, K = decide(labA, labB, 'cross-half EM groups' if a.xem else 'soft-EM groups')
# ---------------- children labels (per project row) and truth readouts
for H, lab, side, k in ((HA, labA, 'A', 0), (HB, labB, 'B', 1)):
    ch = np.zeros(H['nrow'], int)
    for gi, (p, mem) in enumerate(sorted(groups.items(), key=lambda kv: -sum((labA == pairs[m][0]).sum() + (labB == pairs[m][1]).sum() for m in kv[1]))):   # one numbering for both halves (combined size)
        for m in mem:
            ch[H['rows'][lab == pairs[m][k]]] = gi + 1
            for u, t in pair_and_test.extra[k].items():
                if t == pairs[m][k]: ch[H['rows'][lab == u]] = gi + 1
    np.savetxt('%s/children_%s.txt' % (a.out, side), ch, fmt='%d')
    cons_ = consA if k == 0 else consB; lowm = (sh >= 2) & (sh < klo); jl = []   # junk-like flag: a state's map vs the node consensus on the coarse shells (real states share the particle)
    for g in range(1, K + 1):
        sel = np.where(ch[H['rows']] == g)[0]
        if len(sel) < 50: jl.append(float('nan')); continue
        Mg = rec(H, rng.choice(sel, min(len(sel), neq), replace=False), eps * min(len(sel), neq) / neq); x, y = Mg[lowm], cons_[lowm]
        jl.append(float(np.real(np.vdot(x, y)) / max(np.linalg.norm(x) * np.linalg.norm(y), 1e-30)))
    log(' half %s states: %s' % (side, ' '.join('%d:(n=%d, consensus corr %.2f%s)' % (g, (ch == g).sum(), jl[g - 1], ', junk-like' if jl[g - 1] < a.junkcorr else '') for g in range(1, K + 1))))
    if a.gt:
        conf = np.loadtxt(a.gt.split(',')[k], comments='#')[:, 1].astype(int); line = []
        for g in range(1, K + 1):
            cc = np.bincount(conf[ch == g]); nj = cc[0]; cc[0] = 0; o = np.argsort(-cc)[:3]; nn = max((ch == g).sum() - nj, 1)
            line.append('%d:(n=%d, %s)' % (g, (ch == g).sum(), ' '.join('%d:%.2f' % (j, cc[j] / nn) for j in o if cc[j] > 0)))
        log(' half %s children: %s' % (side, ' '.join(line)))
