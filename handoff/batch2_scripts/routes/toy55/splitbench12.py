#!/usr/bin/env python3
"""splitbench12 = splitbench11 + --dock (every group map, in the EM references and in the test, translated rigidly in 3-D to best match
its half's consensus over the EM band: groups that differ only by a translation, e.g. particles sorted by the direction of their shift
errors, are the same structure).
splitbench11 = splitbench10 + junk screen before the test (--junk_proj / --junk_ssnr: a unit whose signal projects on the node consensus
below JUNK_PROJ in the EM band or beyond it, or whose per-particle SSNR beyond the EM band is below JUNK_SSNR of the node's, is set
aside before any merge; junk_units.txt) + --specmatch 3 (pairwise, attenuation only: per shell the group with more signal is scaled
down to the other's amplitude, so noise is never amplified).
splitbench10 = splitbench9 + --specmatch (every group's half-maps scaled per shell to the node consensus signal spectrum from the
cross-half spectra before any difference: blur / contrast / CTF-envelope differences between groups cancel), per-group junk metrics
(signal power, projection on the consensus and per-particle SSNR relative to the node, per band) and --summary.
splitbench9 = splitbench8 + --strata: the label-shuffle null permutes labels only within strata of viewing direction (Fibonacci-sphere
bins of (e1, e2)) x defocus terciles, so every shuffled group keeps its group's view and defocus mix and a split explained by them
is not significant (view-biased subsets differ reproducibly in BOTH halves). Otherwise as splitbench8:
Fast in-memory split test of one tree node on a pair of independent half-sets (splitbench3 + --xem: the re-assignment
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
ap.add_argument('--strata', type=int, default=0, help='1: shuffle null stratified by view (NVIEW direction bins) x defocus terciles'); ap.add_argument('--nview', type=int, default=40)
ap.add_argument('--specmatch', type=int, default=0, help='1: group half-maps scaled per shell to the consensus signal power (cross-half spectra) before differencing; 2: per-shell regression on the consensus')
ap.add_argument('--rb_margin', type=int, default=2, help='shells between the EM band and the re-test band'); ap.add_argument('--dock', type=int, default=0, help='1: rigid 3-D translation of every group map to its half consensus (EM band) before use'); ap.add_argument('--shells', type=int, default=0, help='1: log the per-shell difference FSC of the four most different pairs'); ap.add_argument('--fulln', type=int, default=0, help='1: every group map from all its particles (regulariser scaled with the size) instead of the smallest group\'s N; merges weight by particles'); ap.add_argument('--junk_proj', type=float, default=0.); ap.add_argument('--junk_ssnr', type=float, default=0.)
ap.add_argument('--summary', default='', help='file: K*, child sizes, (with --gt) NMI(conformer; child) over both halves, junk metrics per child')
ap.add_argument('--gt', default=''); ap.add_argument('--truth', default=''); ap.add_argument('--crop', type=int, default=48); ap.add_argument('--nshuf', type=int, default=20)
ap.add_argument('--q', type=float, default=1.0); ap.add_argument('--bandlo', type=float, default=40.); ap.add_argument('--mskdiam', type=float, default=230.)
ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--kv', type=float, default=300.); ap.add_argument('--cs', type=float, default=2.7); ap.add_argument('--fraca', type=float, default=0.1)
ap.add_argument('--conv', default='1,1,1,1,1,1', help='transpose, e1 sign, y flip, e3 sign, shift sign, ctf sign (arm_llr.py calibration of these sets)')
ap.add_argument('--focus', default='', help='real-space mask (box of the stacks, e.g. a region board focus mask): every map is compared only inside it (local differences such as one arm)'); ap.add_argument('--force_em', type=int, default=0, help='1: the cross-half EM starts from every candidate group (no merges before it; a null test of the EM and re-test)'); ap.add_argument('--retest_beyond_em', type=int, default=0, help='1: the re-test after the cross-half EM uses only shells beyond the EM band (+1 margin), with its own shuffle null')
ap.add_argument('--pair', default='hungarian', help="'identity': candidate label k of half A is label k of half B (states defined on both halves by one earlier run)"); ap.add_argument('--protect', default='', help='colon list of candidate labels that never merge with or are pruned into unprotected groups and keep exactly their own particles in the cross-half EM'); ap.add_argument('--ctfabs', type=int, default=0, help='1: |CTF| for phase-flipped stacks (real RNase P)'); ap.add_argument('--xem', type=int, default=0, help='> 0: cross-half EM iterations replace the per-half soft EM'); ap.add_argument('--em_lim', type=float, default=20., help='cross-half EM band limit (A)'); ap.add_argument('--em_min', type=int, default=100, help='a group below this many particles in either half is dropped'); ap.add_argument('--eps', type=float, default=0.05); ap.add_argument('--seed', type=int, default=0); ap.add_argument('--arch', type=int, default=1); ap.add_argument('--junkcorr', type=float, default=0.8, help='states whose coarse-shell map correlation with the node consensus is below this are flagged junk-like'); ap.add_argument('--softem', type=int, default=8, help='soft-EM iterations refining the accepted groups as a mixture before the final test (0 = off)'); ap.add_argument('--archmode', default='shell', help='global: one convex blend | shell: a convex blend per Fourier shell (absorbs the radially filtered imprint of hard classification along a real axis)'); ap.add_argument('--cache', default='')
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
    if a.ctfabs: C = np.abs(C)
    ramp = np.exp(ssh * 2j * np.pi * (G['fx2'][None] * pp[:, 3, None, None] + G['fy2'][None] * pp[:, 4, None, None])).astype(np.complex64)
    Wv = (C * np.conj(ramp) * F).reshape(len(rows), -1); Wd = (C ** 2).reshape(len(rows), -1); F2 = (np.abs(F) ** 2).reshape(len(rows), -1).astype(np.float32)
    Rm = rotm(s1 * pp[:, 0], pp[:, 1], s3 * pp[:, 2]); kk = G['kk0'] * np.array([1., yf, 1.])
    p3 = np.einsum('qj,nij->nqi', kk, Rm) if convT else np.einsum('qj,nji->nqi', kk, Rm); q = np.round(p3).astype(np.int32)
    ok = (np.abs(q) <= c // 2 - 1).all(2); q = q % c; idx = (q[..., 2] * c * c + q[..., 1] * c + q[..., 0]).astype(np.int64); idx[~ok] = c ** 3
    log('%s: %d node particles (labels %s), crop %d (Nyquist %.1f A)' % (tag, len(rows), np.bincount(lab[rows]).tolist()[1:], c, nx * a.smpd / (c // 2)))
    return dict(lab=lab, rows=rows, Wv=Wv, Wd=Wd, F2=F2, idx=idx, G=G, nrow=len(lab), pp=pp)
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
            sel = np.where(lab == g)[0]
            if a.fulln: M[g] = rec(H, sel, eps * len(sel) / neq)   # every group from all its particles, regulariser scaled with its size
            else: sel = rng.choice(sel, neq, replace=False); M[g] = rec(H, sel, eps)
            if a.dock: M[g], group_maps.shifts[(H is HB, g)] = dock(M[g], consB if H is HB else consA)
        out.append(M)
    return out
group_maps.shifts = {}
def dock(M, C):   # rigid 3-D translation of M that best matches C over the EM band (Newton on the cross-correlation), applied to all of M
    m_ = DBAND; a_ = M[:c ** 3][m_] * np.conj(C[:c ** 3][m_]); k_ = KV[m_]; t = np.zeros(3)
    for _ in range(8):
        z = a_ * np.exp(-2j * np.pi * (k_ @ t)); g_ = 2 * np.pi * (k_ * z.imag[:, None]).sum(0); H_ = -(2 * np.pi) ** 2 * (k_.T * z.real) @ k_
        try: st_ = np.linalg.solve(H_, g_)
        except np.linalg.LinAlgError: break
        t = t - st_
        if np.abs(st_).max() < 1e-4: break
    if not np.all(np.isfinite(t)) or np.abs(t).max() > 3.: return M, np.zeros(3)
    out = M.copy(); out[:c ** 3] = M[:c ** 3] * np.exp(-2j * np.pi * (KV @ t)); return out, t
def pair_and_test(MA, MB, cons_A, cons_B, sh, ns, klo, khi):
    ga, gb = sorted(MA), sorted(MB); low = (sh >= 2) & (sh < klo - 1)
    def lowcorr(X, Y): x, y = X[low], Y[low]; return float(np.real(np.vdot(x, y)) / max(np.linalg.norm(x) * np.linalg.norm(y), 1e-30))
    S = np.array([[lowcorr(MA[i] - cons_A, MB[j] - cons_B) for j in gb] for i in ga])
    if a.pair == 'identity': common = [g for g in ga if g in gb]; ia = np.array([ga.index(g) for g in common], int); ib = np.array([gb.index(g) for g in common], int)
    else: ia, ib = linear_sum_assignment(-S)
    pairs = [(ga[i], gb[j]) for i, j in zip(ia, ib)]
    extra = [{}, {}]   # unpaired groups (unequal counts): label -> the most similar paired label of the same half
    for side, M, cons, mine in ((0, MA, cons_A, [p[0] for p in pairs]), (1, MB, cons_B, [p[1] for p in pairs])):
        for u in [g for g in M if g not in mine]: extra[side][u] = max(mine, key=lambda m: lowcorr(M[u] - cons, M[m] - cons))
    pair_and_test.extra = extra
    def stat(DA, DB): f = shell_fsc(DA, DB, sh, ns); return max(float(f[k:k + 3].mean()) for k in range(klo, max(khi - 2, klo) + 1))
    SA, SB = {}, {}
    for p, (ga_, gb_) in enumerate(pairs): SA[p], SB[p] = specm(MA[ga_], MB[gb_])
    pair_and_test.scaled = (SA, SB)
    def pstat(Ap, Aq, Bp, Bq):
        if a.specmatch != 3: return stat(Ap - Aq, Bp - Bq)
        fp, fq = pairmatch(Ap, Bp, Aq, Bq); return stat(Ap * fp - Aq * fq, Bp * fp - Bq * fq)
    pair_and_test.pstat = pstat
    R = {(p, q): pstat(SA[p], SA[q], SB[p], SB[q]) for p in range(len(pairs)) for q in range(p + 1, len(pairs))}
    return pairs, R, stat, float(np.median([S[i, j] for i, j in zip(ia, ib)]))
def pairmatch(Ap, Bp, Aq, Bq):   # attenuation only: per shell the group with more signal is scaled down to the other's amplitude (noise never amplified)
    fl = 1e-6 * max(np.abs(np.bincount(sh, np.real(consA * np.conj(consB)), ns)[2:]).max(), 1e-30); sm3 = lambda v: np.convolve(np.pad(v, 1, mode='edge'), np.ones(3) / 3., mode='valid')
    ap = sm3(np.sqrt(np.maximum(np.bincount(sh, np.real(Ap * np.conj(Bp)), ns), fl))); aq = sm3(np.sqrt(np.maximum(np.bincount(sh, np.real(Aq * np.conj(Bq)), ns), fl))); t = np.minimum(ap, aq)
    return (t / ap)[sh], (t / aq)[sh]
def specm(XA, XB):   # one group's half-maps scaled per shell to the consensus signal spectrum (cross-half estimates: the noise does not enter)
    if a.specmatch not in (1, 2): return XA, XB
    pc = np.bincount(sh, np.real(consA * np.conj(consB)), ns); fl = 1e-3 * np.abs(pc[2:]).max()
    if a.specmatch == 1: s = np.sqrt(np.maximum(np.bincount(sh, np.real(XA * np.conj(XB)), ns), fl) / np.maximum(pc, fl))
    else: s = 0.5 * (np.bincount(sh, np.real(XA * np.conj(consB)), ns) + np.bincount(sh, np.real(XB * np.conj(consA)), ns)) / np.maximum(pc, fl)
    s = np.clip(np.convolve(np.pad(s, 1, mode='edge'), np.ones(3) / 3., mode='valid'), 0.25, 4.); s[:2] = 1.
    return XA / s[sh], XB / s[sh]
HA = load_half(*[s.split(',')[0] for s in (a.sim, a.poses, a.cand)], 'A'); HB = load_half(*[s.split(',')[1] for s in (a.sim, a.poses, a.cand)], 'B')
if a.focus:   # region focus: the mask (stack box) resampled to the crop grid multiplies the sphere of every reconstruction
    from scipy.ndimage import zoom as _zoom
    _h = np.fromfile(a.focus, np.int32, 256); _m = np.fromfile(a.focus, np.float32, offset=1024 + int(_h[23]), count=int(_h[0] * _h[1] * _h[2])).reshape(int(_h[2]), int(_h[1]), int(_h[0]))
    _f = np.clip(_zoom(_m, HA['G']['c'] / _m.shape[0], order=1), 0, 1).astype(np.float32)[:HA['G']['c'], :HA['G']['c'], :HA['G']['c']]
    for _H in (HA, HB): _H['G']['msk3'] = _H['G']['msk3'] * _f
    log('focus mask %s: %.1f %% of the crop grid above 0.5' % (a.focus, 100. * np.mean(_f > 0.5)))
G = HA['G']; c = G['c']; n = G['n']; sh = G['sh3']; ns = c // 2 + 1; klo = max(int(round(n * a.smpd / a.bandlo)), 3); khi = c // 2 - 1
_q = np.arange(c ** 3); KV = np.stack([_q % c, (_q // c) % c, _q // (c * c)], 1); KV = np.where(KV >= c // 2, KV - c, KV).astype(np.float64) / n   # (x, y, z) frequency of every grid voxel, cycles per pixel
DBAND = (sh >= 2) & (sh <= (c // 2 - 1 if a.em_lim <= 0 else min(c // 2 - 1, int(round(n * a.smpd / a.em_lim)))))
labA, labB = HA['lab'][HA['rows']], HB['lab'][HB['rows']]
# one fixed per-shell regulariser: eps x the median per-shell weight of an N_eq reconstruction (same for every group map)
neq = int(min(np.bincount(labA)[1:][np.bincount(labA)[1:] > 0].min(), np.bincount(labB)[1:][np.bincount(labB)[1:] > 0].min()))
def shell_eps(H):
    c_ = H['G']['c']; ii = H['idx'].ravel(); den = np.bincount(ii, H['Wd'].ravel(), c_ ** 3 + 1)[:-1] * (neq / len(H['rows']))
    return np.array([np.median(den[(sh == k) & (den > 0)]) if ((sh == k) & (den > 0)).any() else 1. for k in range(ns)]) * a.eps
eps = 0.5 * (shell_eps(HA) + shell_eps(HB)); consA, consB = rec(HA, np.arange(len(HA['rows'])), eps * len(HA['rows']) / neq), rec(HB, np.arange(len(HB['rows'])), eps * len(HB['rows']) / neq)
log('equal N %d per group; band %d..%d (%.1f..%.1f A); %d shuffles, threshold quantile %.3f' % (neq, klo, khi, n * a.smpd / klo, n * a.smpd / khi, a.nshuf, a.q))
GTC = np.loadtxt(a.gt.split(',')[0], comments='#')[:, 1].astype(int) if a.gt else None
def gmetrics(selA, selB, cA=None, cB=None, nc=None):
    MAg = rec(HA, selA, eps * max(len(selA), 1) / neq); MBg = rec(HB, selB, eps * max(len(selB), 1) / neq); kb = min(c // 2 - 1, int(round(n * a.smpd / a.em_lim))); out = []
    cA, cB = (consA, consB) if cA is None else (cA, cB)
    for lo, hi in ((2, klo - 1), (klo, kb), (kb + 1, khi)):
        m = (sh >= lo) & (sh <= hi); xa, xb, ca, cb = MAg[m], MBg[m], cA[m], cB[m]
        pg, pc = np.real(np.vdot(xb, xa)), np.real(np.vdot(cb, ca)); xc = 0.5 * (np.real(np.vdot(cb, xa)) + np.real(np.vdot(ca, xb)))
        fg = pg / max(np.sqrt(np.real(np.vdot(xa, xa)) * np.real(np.vdot(xb, xb))), 1e-30); fc = pc / max(np.sqrt(np.real(np.vdot(ca, ca)) * np.real(np.vdot(cb, cb))), 1e-30)
        ng, nc_ = 0.5 * (len(selA) + len(selB)), (0.5 * (len(HA['rows']) + len(HB['rows'])) if nc is None else nc)
        out.append((pg / pc, xc / pc, (max(fg, 0.) / max(1 - fg, 1e-3) / ng) / max(fc / max(1 - fc, 1e-3) / nc_, 1e-12)))
    return out
def comp(rows):
    if GTC is None or len(rows) == 0: return ''
    cc = np.bincount(GTC[rows], minlength=1); nj = cc[0]; cc[0] = 0; o = np.argsort(-cc)[:3]
    return 'junk %.2f, %s' % (nj / len(rows), ' '.join('%d:%.2f' % (j, cc[j] / len(rows)) for j in o if cc[j] > 0))
def mline(g, selA, selB, rows):
    m = gmetrics(selA, selB); return '%s: n %d | power %.2f/%.2f/%.2f proj %.2f/%.2f/%.2f rel.SSNR %.2f/%.2f/%.2f | %s' % (g, len(selA) + len(selB), m[0][0], m[1][0], m[2][0], m[0][1], m[1][1], m[2][1], m[0][2], m[1][2], m[2][2], comp(rows))
log('unit metrics (bands >%.0f / %.0f-%.0f / %.0f-%.0f A: signal power and projection relative to the node consensus, per-particle SSNR relative to the node):' % (a.bandlo, a.bandlo, a.em_lim, a.em_lim, n * a.smpd / khi))
for g in sorted(set(labA[labA > 0]) & set(labB[labB > 0])):
    sA, sB = np.where(labA == g)[0], np.where(labB == g)[0]; log('  unit ' + mline(g, sA, sB, np.r_[HA['rows'][sA], HB['rows'][sB]]))
LAB0 = HA['lab'] + HB['lab'] if not np.intersect1d(HA['rows'], HB['rows']).size else HA['lab'].copy(); JUNK = []
if a.junk_proj > 0 or a.junk_ssnr > 0:   # junk screen before any merge: a unit lacking the node's signal (projection on the consensus or per-particle SSNR beyond the EM band)
    units_ = sorted(set(labA[labA > 0]) & set(labB[labB > 0]))
    def jflags(cA_, cB_, nc=None):
        out = {}
        for g in units_: m = gmetrics(np.where(labA == g)[0], np.where(labB == g)[0], cA_, cB_, nc); out[g] = (min(m[1][1], m[2][1]) < a.junk_proj or m[2][2] < a.junk_ssnr, m)
        return out
    fl_ = jflags(consA, consB); keep = [g for g in units_ if not fl_[g][0]]
    if keep and len(keep) < len(units_):
        sA_, sB_ = np.where(np.isin(labA, keep))[0], np.where(np.isin(labB, keep))[0]; fl_ = jflags(rec(HA, sA_, eps * len(sA_) / neq), rec(HB, sB_, eps * len(sB_) / neq), 0.5 * (len(sA_) + len(sB_)))
    JUNK = [g for g in units_ if fl_[g][0]]
    for g in units_: log('  junk screen unit %d: projection on the screened consensus %.2f / %.2f, rel. SSNR %.2f / %.2f -> %s' % (g, fl_[g][1][1][1], fl_[g][1][2][1], fl_[g][1][1][2], fl_[g][1][2][2], 'JUNK' if g in JUNK else 'ok'))
    open(a.out + '/junk_units.txt', 'w').write(''.join('%d\n' % g for g in JUNK))
    if JUNK and len(JUNK) < len(units_):   # the split test runs on the units that pass
        for H in (HA, HB):
            m = ~np.isin(H['lab'][H['rows']], JUNK); H['lab'] = np.where(np.isin(H['lab'], JUNK), 0, H['lab'])
            for k_ in ('Wv', 'Wd', 'F2', 'idx', 'pp'): H[k_] = H[k_][m]
            H['rows'] = H['rows'][m]
        labA, labB = HA['lab'][HA['rows']], HB['lab'][HB['rows']]; neq = int(min(np.bincount(labA)[1:][np.bincount(labA)[1:] > 0].min(), np.bincount(labB)[1:][np.bincount(labB)[1:] > 0].min()))
        eps = 0.5 * (shell_eps(HA) + shell_eps(HB)); consA, consB = rec(HA, np.arange(len(HA['rows'])), eps * len(HA['rows']) / neq), rec(HB, np.arange(len(HB['rows'])), eps * len(HB['rows']) / neq)
        log('junk screen: units %s set aside; the split test runs on %d / %d particles, equal N %d' % (' '.join(map(str, JUNK)), len(HA['rows']), len(HB['rows']), neq))
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
def strata(H):   # nuisance strata of the node particles: nearest of NVIEW Fibonacci directions to (e1, e2) x defocus tercile
    pp = H['pp']; ph, th = np.radians(pp[:, 0]), np.radians(pp[:, 1]); v = np.c_[np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)]
    i = np.arange(a.nview) + 0.5; z = 1 - 2 * i / a.nview; r = np.sqrt(np.maximum(1 - z * z, 0)); g = np.pi * (3 - np.sqrt(5)) * i
    vb = np.argmax(v @ np.c_[r * np.cos(g), r * np.sin(g), z].T, 1); df = 0.5 * (pp[:, 5] + pp[:, 6]); dq = np.searchsorted(np.quantile(df, [1 / 3., 2 / 3.]), df)
    return vb * 3 + dq
def sshuf(lab, st):
    out = lab.copy()
    for k in np.unique(st): ix = np.where(st == k)[0]; out[ix] = rng.permutation(lab[ix])
    return out
if a.strata: STA, STB = strata(HA), strata(HB); log('stratified null: %d / %d nuisance strata (view x defocus) in halves A / B' % (len(np.unique(STA)), len(np.unique(STB))))
NUL = []; NULB = []
kb_em = c // 2 - 1 if a.em_lim <= 0 else min(c // 2 - 1, int(round(n * a.smpd / a.em_lim))); klo2 = min(kb_em + a.rb_margin, khi - 2); NUL2, NULB2 = [], []   # beyond-EM band: shells the E-step never saw
for s in range(a.nshuf):
    if a.strata: pa, pb = sshuf(labA, STA), sshuf(labB, STB)
    else: pa, pb = labA.copy(), labB.copy(); rng.shuffle(pa); rng.shuffle(pb)
    MA, MB = group_maps(HA, HB, pa, pb, eps, neq)
    prs, R, st_, _ = pair_and_test(MA, MB, consA, consB, sh, ns, klo, khi); NUL += list(R.values())
    gA, gB = pair_and_test.scaled
    if a.retest_beyond_em:
        prs2, R2, st2_, _ = pair_and_test(MA, MB, consA, consB, sh, ns, klo2, khi); NUL2 += list(R2.values())
        if a.arch and len(prs2) > 2:
            gA2, gB2 = pair_and_test.scaled
            for p in gA2: oth = [o for o in gA2 if o != p]; NULB2.append(st2_(blend_resid(gA2, p, oth)[0], blend_resid(gB2, p, oth)[0]))
    if a.arch and len(prs) > 2:
        for p in gA: oth = [o for o in gA if o != p]; NULB.append(st_(blend_resid(gA, p, oth)[0], blend_resid(gB, p, oth)[0]))
NULB = np.array(NULB) if NULB else np.zeros(1); THRB = float(np.quantile(NULB, a.q)); NUL = np.array(NUL); THR = float(np.quantile(NUL, a.q)); log('blend-residual null (%d): median %.3f, 99 %% %.3f, max %.3f -> blend threshold %.3f' % (len(NULB), np.median(NULB), np.quantile(NULB, .99), NULB.max(), THRB)); log('shuffle null (%d pair statistics): median %.3f, 95 %% %.3f, 99 %% %.3f, max %.3f -> threshold %.3f' % (len(NUL), np.median(NUL), np.quantile(NUL, .95), np.quantile(NUL, .99), NUL.max(), THR))
if a.retest_beyond_em:
    THR2 = float(np.quantile(np.array(NUL2), a.q)); THRB2 = float(np.quantile(np.array(NULB2), a.q)) if NULB2 else THRB
    log('beyond-EM band %d..%d (%.1f..%.1f A): shuffle null median %.3f, max %.3f -> threshold %.3f; blend threshold %.3f' % (klo2, khi, n * a.smpd / klo2, n * a.smpd / khi, np.median(NUL2), np.max(NUL2), THR2, THRB2))
def decide(labA, labB, tag, protect=(), band=None):
  global pairs
  klo_, THR, THRB = (klo, globals()['THR'], globals()['THRB']) if band is None else band
  neq_ = int(min(np.bincount(labA)[1:][np.bincount(labA)[1:] > 0].min(), np.bincount(labB)[1:][np.bincount(labB)[1:] > 0].min()))
  MA, MB = group_maps(HA, HB, labA, labB, eps * neq_ / neq, neq_); pairs, R, stat, pc = pair_and_test(MA, MB, consA, consB, sh, ns, klo_, khi); log('[%s] equal N %d%s' % (tag, neq_, '' if klo_ == klo else ', re-test band %d..%d' % (klo_, khi)))
  GA, GB = [dict(d) for d in pair_and_test.scaled]; pst = pair_and_test.pstat; cntA = {p: int((labA == pairs[p][0]).sum()) for p in range(len(pairs))}; cntB = {p: int((labB == pairs[p][1]).sum()) for p in range(len(pairs))}
  if klo_ != klo: log('pair statistics at the start, full band %d..%d (sorted): %s' % (klo, khi, ' '.join('%.2f' % v for v in sorted(pair_and_test(MA, MB, consA, consB, sh, ns, klo, khi)[1].values()))))
  log('paired %d groups (coarse-shell correlation median %.3f): %s' % (len(pairs), pc, ' '.join('%d<->%d' % p for p in pairs)))
  log('pair statistics at the start (sorted): %s' % ' '.join('%.2f' % v for v in sorted(R.values())))
  log('pair statistics at the start: %s' % ' '.join('%d-%d:%.2f' % (pairs[p][0], pairs[q][0], v) for (p, q), v in sorted(R.items())))
  if a.dock: log('  docking shifts (px, half A): %s' % ' '.join('%d:%.2f/%.2f/%.2f' % ((g,) + tuple(t)) for (hb, g), t in sorted(group_maps.shifts.items()) if not hb))
  if a.shells:
      for (p, q), v in sorted(R.items(), key=lambda kv: -kv[1])[:4]: log('  shells 2..%d of pair %d-%d (stat %.2f): %s' % (khi, pairs[p][0], pairs[q][0], v, ' '.join('%.2f' % x for x in shell_fsc(GA[p] - GA[q], GB[p] - GB[q], sh, ns)[2:khi + 1])))
  groups = {p: [p] for p in range(len(pairs))}
  def gstat(p, q): return pst(GA[p], GA[q], GB[p], GB[q])   # the tested band's statistic (pair_and_test is called again for the full-band log)
  R = {k: gstat(*k) for k in R}; merged = []; protg = {p for p in range(len(pairs)) if pairs[p][0] in protect}
  for p in sorted(protg):   # what the data alone say about every protected group (it would merge with its nearest unprotected group if that statistic is below the threshold)
      un = [(gstat(*sorted((p, o))), o) for o in groups if o not in protg]
      if un: v, o = min(un); log('  protected label %s: nearest unprotected label %s at %.3f (threshold %.3f) -> %s without protection' % (pairs[p][0], pairs[o][0], v, THR, 'would MERGE' if v < THR else 'survives'))
  def adm(R): return {k: w for k, w in R.items() if (k[0] in protg) == (k[1] in protg)}   # protected groups merge only with each other
  while adm(R) and min(adm(R).values()) < THR:
      (p, q), v = min(adm(R).items(), key=lambda kv: kv[1]); merged.append(v)
      for D, cn in ((GA, cntA), (GB, cntB)):
          np_, nq = (cn[p], cn[q]) if a.fulln else (len(groups[p]), len(groups[q])); D[p] = (np_ * D[p] + nq * D[q]) / (np_ + nq); del D[q]
      cntA[p] += cntA.pop(q); cntB[p] += cntB.pop(q)
      groups[p] += groups.pop(q); protg.discard(q); R = {k: w for k, w in R.items() if q not in k}; R.update({tuple(sorted((p, o))): gstat(*sorted((p, o))) for o in groups if o != p})
  log('%d merges at %s -> %d groups' % (len(merged), ' '.join('%.2f' % v for v in merged), len(groups)))
  while a.arch and len(groups) > 2:
      best = None
      for p in groups:
          if p in protg: continue   # a protected group is never pruned
          others = [o for o in groups if o != p]; rA, wA = blend_resid(GA, p, others); rB, wB = blend_resid(GB, p, others); v = stat(rA, rB)
          un = [i for i, o in enumerate(others) if o not in protg]   # an unprotected group is never folded into a protected one
          if un and (best is None or v < best[1]): best = (p, v, others[max(un, key=lambda i: (wA + wB)[i])])
      if best is None: break
      log('  archetype test: least-novel group %d blend residual %.2f (blend threshold %.3f)' % (best[0], best[1], THRB))
      if best[1] >= THRB: break
      groups[best[2]] += groups.pop(best[0]); GA.pop(best[0]); GB.pop(best[0])
  K = len(groups); log('K* = %d%s' % (K, '' if not protg else ' (protected: %s)' % ' '.join(str(pairs[p][0]) for p in sorted(protg))))
  decide.protg = set(protg); return groups, pairs, K
PROT = set(int(v) for v in a.protect.replace(',', ':').split(':') if v)
groups, pairs, K = decide(labA, labB, 'board groups', PROT)
if a.force_em: groups = {p: [p] for p in range(len(pairs))}; K = len(groups); log('force_em: the cross-half EM starts from all %d candidate groups' % K)
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
def cross_em(stA, stB, K, iters, frozen=()):   # soft EM at fixed poses, each half scored against the OTHER half's maps: consensus + group deviation x cross-half FSC per shell
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
            if a.dock: out[-1] = dock(out[-1], cons)[0]
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
        for W_, st_ in ((WA, stA), (WB, stB)):   # protected groups keep exactly their own particles
            for g in frozen:
                m = st_ == g; W_[:, g - 1] = 0.; W_[m] = 0.; W_[m, g - 1] = 1.
            rs = W_.sum(1, keepdims=True); bad = rs[:, 0] <= 0
            if bad.any(): W_[bad] = np.where(alive & ~np.isin(np.arange(1, K + 1), list(frozen)), 1., 0.)[None]; rs = W_.sum(1, keepdims=True)
            W_ /= rs
        chA, chB = np.mean(WA.argmax(1) != P[0]['W'].argmax(1)), np.mean(WB.argmax(1) != P[1]['W'].argmax(1)); P[0]['W'], P[1]['W'] = WA, WB
        cA_, cB_ = np.bincount(WA.argmax(1), minlength=K), np.bincount(WB.argmax(1), minlength=K)   # MAP counts: the re-test's equal N comes from them
        for g in np.where(alive)[0]:
            if (g + 1) in frozen: continue
            if min(cA_[g], cB_[g]) < a.em_min: alive[g] = False; log('  group %d drained (%d / %d MAP particles) and dropped' % (g + 1, cA_[g], cB_[g]))
        s = '  iteration %d: MAP changed A %.3f B %.3f; weights A %s; deviation FSC per group (shell 4/8/12) %s' % (it + 1, chA, chB, ' '.join('%.3f' % v for v in WA.mean(0)), ' '.join('%.2f/%.2f/%.2f' % (f[4], f[min(8, ns - 1)], f[min(12, ns - 1)]) for f, al in zip(fl, alive) if al))
        if a.gt: s += '; NMI(conformer; group) A %.3f B %.3f' % (nmi_(P[0]['conf'], WA.argmax(1)), nmi_(P[1]['conf'], WB.argmax(1)))
        log(s)
    return [np.where(alive[None], P[h]['W'], -1.).argmax(1) + 1 for h in (0, 1)]   # a group drained in the last iteration keeps no particles
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
        FROZ = [gi + 1 for gi, p in enumerate(groups) if p in decide.protg]; log('cross-half EM: protected (frozen) groups %s' % FROZ) if FROZ else None
        labA_em, labB_em = cross_em(starts[0], starts[1], len(groups), a.xem, FROZ)
        for k, nl in ((0, labA_em), (1, labB_em)): log('cross-half EM half %s: MAP changed vs start %.3f' % ('AB'[k], np.mean(nl != starts[k])))
    else:
        for H, k in ((HA, 0), (HB, 1)):   # per half, independently
            newlab, W = soft_em(H, starts[k], a.softem); log('soft EM half %s: weights %s, MAP changed vs start %.3f' % ('AB'[k], ' '.join('%.3f' % v for v in W.mean(0)), np.mean(newlab != starts[k])))
            if k == 0: labA_em = newlab
            else: labB_em = newlab
    labA, labB = labA_em, labB_em; groups, pairs, K = decide(labA, labB, 'cross-half EM groups' if a.xem else 'soft-EM groups', set(FROZ) if a.xem else set(), (klo2, THR2, THRB2) if (a.xem and a.retest_beyond_em) else None)
# ---------------- children labels (per project row) and truth readouts
for H, lab, side, k in ((HA, labA, 'A', 0), (HB, labB, 'B', 1)):
    ch = np.zeros(H['nrow'], int)
    for gi, (p, mem) in enumerate(sorted(groups.items(), key=lambda kv: -sum((labA == pairs[m][0]).sum() + (labB == pairs[m][1]).sum() for m in kv[1]))):   # one numbering for both halves (combined size)
        if k == 0 and p in getattr(decide, 'protg', set()): log('protected group (EM label %s) -> final state %d' % (pairs[p][0], gi + 1))
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
def nmi_sc(x, y):   # scorecard.py's NMI (sqrt normalisation)
    _, xi = np.unique(x, return_inverse=True); _, yi = np.unique(y, return_inverse=True); J = np.zeros((xi.max() + 1, yi.max() + 1)); np.add.at(J, (xi, yi), 1); J /= J.sum()
    px, py = J.sum(1), J.sum(0); nz = J > 0; I = (J[nz] * np.log(J[nz] / np.outer(px, py)[nz])).sum(); return I / max(np.sqrt((px * np.log(px)).sum() * (py * np.log(py)).sum()), 1e-12)
chs = [np.loadtxt('%s/children_%s.txt' % (a.out, s_), dtype=int) for s_ in 'AB']; one = not np.intersect1d(HA['rows'], HB['rows']).size; chall = chs[0] + chs[1] if one else chs[0]
lines = ['  junk unit %d (set aside before the test)' % g for g in JUNK] + ['  child ' + mline(g, np.where(chs[0][HA['rows']] == g)[0], np.where(chs[1][HB['rows']] == g)[0], np.where(chall == g)[0] if one else np.array([], int)) for g in range(1, K + 1)]
log('child metrics:'); [log(l_) for l_ in lines]; nm = ''
if GTC is not None and one:
    u0 = LAB0; on = u0 > 0; cj = chall.copy()
    for i_, g in enumerate(JUNK): cj[u0 == g] = K + 1 + i_
    nm = ' NMI(conformer; unit) %.3f -> NMI(conformer; child%s) %.3f' % (nmi_sc(GTC[on], u0[on]), ' + junk units' if JUNK else '', nmi_sc(GTC[cj > 0], cj[cj > 0]))
jz = ' junk units %s (%d particles)' % ('+'.join(map(str, JUNK)), int(np.isin(LAB0, JUNK).sum())) if JUNK else ''
log('SUMMARY K* %d sizes %s%s%s' % (K, ' '.join(str((chall == g).sum()) for g in range(1, K + 1)), jz, nm))
if a.summary: open(a.summary, 'w').write('K* %d sizes %s%s%s\n' % (K, ' '.join(str((chall == g).sum()) for g in range(1, K + 1)), jz, nm) + '\n'.join(lines) + '\n')
