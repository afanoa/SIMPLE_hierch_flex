#!/usr/bin/env python3
"""Minority states: soft (responsibility-weighted, learned mixing weights) vs hard (argmax, classification EM) assignment of a node's
particles to its states at fixed poses, every particle scored against the OTHER half's references (consensus + the state's deviation
weighted by its cross-half FSC, the deviation restricted to the focus mask in real space), each reference scaled per shell by the data's
regression on it. With hard assignment a minority state's map is a blend of its members and its neighbours' tail; with soft assignment
every particle contributes to every state in proportion to its responsibility (the ML estimate). Writes each state's map (both halves
averaged, padded to box BOX, in SIMPLE's frame: the signed axis permutation that best matches --cons_ref) for the start labels and after
the EM: OUT/maps_start/state_NN.mrc, OUT/maps_<mode>/state_NN.mrc, plus OUT/labels_<mode>_{A,B}.txt (argmax per half row).
minor_em2.py (pipeline copy): --frame identity (default; r22 found SIMPLE's frame = identity, corr 0.907) skips the axis search, --nomaps 1
skips writing the state maps (the pipeline reconstructs the argmax labels with SIMPLE).
minor_em3.py: --post 1 also writes OUT/post_<mode>_{A,B}.txt, every component's responsibility per half row (0 outside the node), for
selecting a component's most confident particles (site-occupancy classification: 2 components that differ only inside a site mask).
minor_em4.py: --empty_init 1 = occupancy of a site FILLED in most particles: the labels only mark the node (all class 1); in the first
E-step class 1 = the consensus, class 2 = the consensus with the site's density removed (x (1 - mask)), mixing weights 0.9 / 0.1;
afterwards both classes are rebuilt from their responsibilities as usual (class 2 = site empty).
usage: minor_em4.py --sim SIM_A,SIM_B --poses P_A,P_B --labels L_A,L_B --cache C --cons_ref CONS.mrc --out DIR [--mask M.mrc --soft 1 --iters 8
       --em_lim 12 --ctfabs 1 --gt GT_A,GT_B --closed_col 2]"""
import os, itertools, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--sim', required=True); ap.add_argument('--poses', required=True); ap.add_argument('--labels', required=True)
ap.add_argument('--cache', default=''); ap.add_argument('--cons_ref', required=True); ap.add_argument('--out', required=True); ap.add_argument('--mask', default='')
ap.add_argument('--soft', type=int, default=1); ap.add_argument('--iters', type=int, default=8); ap.add_argument('--em_lim', type=float, default=12.); ap.add_argument('--box', type=int, default=160)
ap.add_argument('--gt', default=''); ap.add_argument('--closed_col', type=int, default=2); ap.add_argument('--crop', type=int, default=48); ap.add_argument('--mskdiam', type=float, default=230.)
ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--kv', type=float, default=300.); ap.add_argument('--cs', type=float, default=2.7); ap.add_argument('--fraca', type=float, default=0.1)
ap.add_argument('--conv', default='1,1,1,1,1,1'); ap.add_argument('--ctfabs', type=int, default=1); ap.add_argument('--eps', type=float, default=1.); ap.add_argument('--pifloor', type=float, default=1e-3); ap.add_argument('--frame', default='identity'); ap.add_argument('--nomaps', type=int, default=0); ap.add_argument('--post', type=int, default=0); ap.add_argument('--empty_init', type=int, default=0)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); mode = 'soft' if a.soft else 'hard'
cv = [float(v) for v in a.conv.split(',')]; convT, s1, yf, s3, ssh, sctf = bool(cv[0]), cv[1], cv[2], cv[3], cv[4], cv[5]
def log(*s): print(*s, flush=True)
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
def rd_mrc(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
def wr_mrc(p, x, smpd):
    h = np.zeros(256, np.int32); h[0:3] = x.shape; h[3] = 2; h[7:10] = x.shape; hf = h.view(np.float32); hf[10:13] = np.array(x.shape) * smpd; hf[13:16] = 90.
    h[16:19] = [1, 2, 3]; hf[19], hf[20], hf[21] = x.min(), x.max(), x.mean(); h[52] = int.from_bytes(b'MAP ', 'little'); h[53] = 0x4444
    with open(p, 'wb') as f: f.write(h.tobytes()); f.write(x.astype(np.float32).tobytes())
lam = 12.2643247 / np.sqrt(a.kv * 1e3 * (1. + a.kv * 1e3 * 0.978466e-6))
def rotm(e1, e2, e3):
    a1, a2, a3 = np.radians(e1), np.radians(e2), np.radians(e3); c1, si1, c2, si2, c3, si3 = np.cos(a1), np.sin(a1), np.cos(a2), np.sin(a2), np.cos(a3), np.sin(a3)
    Z1 = np.zeros((len(e1), 3, 3)); Z1[:, 0, 0] = c1; Z1[:, 0, 1] = -si1; Z1[:, 1, 0] = si1; Z1[:, 1, 1] = c1; Z1[:, 2, 2] = 1
    Y2 = np.zeros((len(e1), 3, 3)); Y2[:, 0, 0] = c2; Y2[:, 0, 2] = si2; Y2[:, 2, 0] = -si2; Y2[:, 2, 2] = c2; Y2[:, 1, 1] = 1
    Z3 = np.zeros((len(e1), 3, 3)); Z3[:, 0, 0] = c3; Z3[:, 0, 1] = -si3; Z3[:, 1, 0] = si3; Z3[:, 1, 1] = c3; Z3[:, 2, 2] = 1
    return Z1 @ Y2 @ Z3
# ---------------- geometry of the cropped grid (as splitbench13) and the scoring band
SIMS = a.sim.split(','); h0 = np.frombuffer(open(SIMS[0] + '/all.mrcs', 'rb').read(1024), np.int32); n = int(h0[0]); c = a.crop
x = np.arange(n) - n // 2; R0 = a.mskdiam / 2 / a.smpd; msk2 = np.clip((R0 + 4 - np.hypot(x[:, None], x[None, :])) / 8., 0, 1).astype(np.float32); sl = np.r_[0:c // 2, n - c // 2:n]
fc = np.fft.fftfreq(c) * c / n; fy2, fx2 = np.meshgrid(fc, fc, indexing='ij'); kk0 = np.stack([fx2, fy2, np.zeros_like(fx2)], -1).reshape(-1, 3) * n
f3 = np.fft.fftfreq(c) * c / n; kx, ky, kz = np.meshgrid(f3, f3, f3, indexing='ij'); sh = np.minimum(np.round(np.sqrt(kx ** 2 + ky ** 2 + kz ** 2) * n).astype(int), c // 2).ravel()
ns = c // 2 + 1; m3 = c ** 3 + 1; kb = min(c // 2 - 1, int(round(n * a.smpd / a.em_lim)))
r2 = np.round(np.hypot(kk0[:, 0], kk0[:, 1])).astype(int); cols = np.where((r2 >= 2) & (r2 <= kb))[0]
kf = (np.hypot(fx2, fy2) / a.smpd).ravel()[cols]; ang = np.arctan2(fy2, fx2).ravel()[cols]
g3 = np.arange(c) - c // 2; r3 = np.sqrt(g3[:, None, None] ** 2 + g3[None, :, None] ** 2 + g3[None, None, :] ** 2) * n / c; msk3 = np.clip((R0 + 4 - r3) / 8., 0, 1)
log('box %d, crop %d (Nyquist %.1f A), band shells 2..%d (to %.1f A), %d pixels per image' % (n, c, n * a.smpd / (c // 2), kb, n * a.smpd / kb, len(cols)))
def load_half(simdir, posef, lab, tag):   # only the band's pixels, in blocks (no full-crop arrays)
    rows = np.where(lab > 0)[0]
    Pz = np.array([[float(kv(l).get(k, 0)) for k in ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast')] for l in open(posef) if l.strip()])[rows]
    hh = np.frombuffer(open(simdir + '/all.mrcs', 'rb').read(1024), np.int32); nz, nsym = int(hh[2]), int(hh[23])
    cache = '%s_%s_c%d.npy' % (a.cache, tag, c) if a.cache else ''; Fc = np.load(cache, mmap_mode='r') if cache and os.path.isfile(cache) else None
    I = None if Fc is not None else np.memmap(simdir + '/all.mrcs', np.float32, 'r', offset=1024 + nsym, shape=(nz, n, n))
    N = len(rows); Wv = np.zeros((N, len(cols)), np.complex64); Wd = np.zeros((N, len(cols)), np.float32); F2 = np.zeros((N, len(cols)), np.float32); ix = np.zeros((N, len(cols)), np.int32)
    for b in range(0, N, 10000):
        rr = rows[b:b + 10000]; pp = Pz[b:b + 10000]
        if Fc is not None: F = np.asarray(Fc[rr]).reshape(len(rr), -1)[:, cols]
        else: blk = np.array(I[rr], np.float32) * msk2; F = np.fft.fft2(np.fft.ifftshift(blk, axes=(1, 2)))[:, sl][:, :, sl].reshape(len(rr), -1)[:, cols]
        df = 0.5 * ((pp[:, 5] + pp[:, 6])[:, None] + (pp[:, 5] - pp[:, 6])[:, None] * np.cos(2. * (ang[None] - np.radians(pp[:, 7])[:, None]))) * 1e4
        chi = np.pi * lam * df * kf[None] ** 2 - 0.5 * np.pi * a.cs * 1e7 * lam ** 3 * kf[None] ** 4; C = (sctf * (np.sqrt(1. - a.fraca ** 2) * np.sin(chi) + a.fraca * np.cos(chi))).astype(np.float32)
        if a.ctfabs: C = np.abs(C)
        ramp = np.exp(ssh * 2j * np.pi * (fx2.ravel()[cols][None] * pp[:, 3, None] + fy2.ravel()[cols][None] * pp[:, 4, None])).astype(np.complex64)
        Wv[b:b + len(rr)] = C * np.conj(ramp) * F; Wd[b:b + len(rr)] = C ** 2; F2[b:b + len(rr)] = np.abs(F) ** 2
        Rm = rotm(s1 * pp[:, 0], pp[:, 1], s3 * pp[:, 2]); kk = kk0[cols] * np.array([1., yf, 1.])
        p3 = np.einsum('qj,nij->nqi', kk, Rm) if convT else np.einsum('qj,nji->nqi', kk, Rm); q = np.round(p3).astype(np.int32)
        ok = (np.abs(q) <= c // 2 - 1).all(2); q = q % c; ixb = q[..., 2] * c * c + q[..., 1] * c + q[..., 0]; ixb[~ok] = c ** 3; ix[b:b + len(rr)] = ixb
    log('%s: %d node particles, labels %s' % (tag, N, np.bincount(lab[rows]).tolist()[1:]))
    return dict(rows=rows, Wv=Wv, Wd=Wd, F2=F2, ix=ix, shq=np.minimum(np.r_[sh, 0][ix], ns - 1).astype(np.int16), lab=lab[rows])
LABS = [np.loadtxt(f, dtype=int) for f in a.labels.split(',')]; K = int(max(l.max() for l in LABS))
if a.empty_init: LABS = [np.minimum(l, 1) for l in LABS]; K = 2
D = [load_half(SIMS[h], a.poses.split(',')[h], LABS[h], 'AB'[h]) for h in (0, 1)]
for h in (0, 1):
    d = D[h]; W = np.zeros((len(d['rows']), K), np.float32); W[np.arange(len(d['rows'])), d['lab'] - 1] = 1.; d['W'] = W
    if a.gt:
        G = np.loadtxt(a.gt.split(',')[h], comments='#'); cl = np.zeros(len(LABS[h]), bool); gi = G[:, 0].astype(int) - 1; okg = gi < len(cl)
        cl[gi[okg]] = G[okg, a.closed_col] > 0; d['closed'] = cl[d['rows']]
neq = int(np.median([max((D[h]['W'][:, g] > 0).sum(), 1) for h in (0, 1) for g in range(K)]))
def shell_eps(d):
    den = np.bincount(d['ix'].ravel(), d['Wd'].ravel(), m3)[:-1] * (neq / len(d['W'])); return np.array([np.median(den[(sh == k) & (den > 0)]) if ((sh == k) & (den > 0)).any() else 1. for k in range(ns)]) * a.eps
eps = 0.5 * (shell_eps(D[0]) + shell_eps(D[1])); epsv = np.r_[eps[sh], eps[-1]]
def wmap(d, w):   # Fourier map (FFT layout, c^3 + 1 with the out-of-band bin last) of the particles weighted by w, regulariser scaled with their number
    ii = d['ix'].ravel(); num = np.bincount(ii, (d['Wv'].real * w[:, None]).ravel(), m3) + 1j * np.bincount(ii, (d['Wv'].imag * w[:, None]).ravel(), m3)
    V = num / (np.bincount(ii, (d['Wd'] * w[:, None]).ravel(), m3) + epsv * max(float(w.sum()), 1.) / neq); V[-1] = 0.; return V
def to_real(V): return np.fft.fftshift(np.real(np.fft.ifftn(V[:-1].reshape(c, c, c))))
def to_four(v): return np.r_[np.fft.fftn(np.fft.ifftshift(v)).ravel(), 0.]
def crop_real(v):   # box-n real-space map -> the c^3 grid (Fourier crop), same physical frame
    F = np.fft.fftn(np.fft.ifftshift(v)); return np.fft.fftshift(np.real(np.fft.ifftn(F[np.ix_(sl, sl, sl)]))) * (c / v.shape[0]) ** 3
def pad_real(V, N):   # c^3 Fourier map -> box N real-space map
    F = np.zeros((N, N, N), complex); sN = np.r_[0:c // 2, N - c // 2:N]; F[np.ix_(sN, sN, sN)] = V[:-1].reshape(c, c, c); return np.fft.fftshift(np.real(np.fft.ifftn(F))) * (N / c) ** 3
def sperm(v, p, s):   # signed axis permutation about the centre (index N//2 stays fixed)
    y = np.transpose(v, p)
    for ax in range(3):
        if s[ax] < 0: y = np.roll(np.flip(y, ax), 1, ax)
    return y
def sperm_inv(v, p, s):
    y = v
    for ax in range(3):
        if s[ax] < 0: y = np.roll(np.flip(y, ax), 1, ax)
    return np.transpose(y, np.argsort(p))
# ---------------- frame: the signed axis permutation taking this script's reconstruction to SIMPLE's consensus
cons48 = 0.5 * (to_real(wmap(D[0], np.ones(len(D[0]['W']), np.float32))) + to_real(wmap(D[1], np.ones(len(D[1]['W']), np.float32))))
ref48 = crop_real(rd_mrc(a.cons_ref)); inm = msk3 > 0.5
def zc(u): u = u[inm] - u[inm].mean(); return u / max(np.linalg.norm(u), 1e-30)
cands = [((0, 1, 2), (1, 1, 1))] if a.frame == 'identity' else [(p, s) for p in itertools.permutations(range(3)) for s in itertools.product((1, -1), repeat=3)]
best = max(((float(zc(sperm(cons48, p, s)) @ zc(ref48)), p, s) for p, s in cands), key=lambda t: t[0])
CC, PERM, SGN = best; log('frame: axis permutation %s signs %s, correlation with %s %.3f' % (PERM, SGN, os.path.basename(a.cons_ref), CC))
MSK = None
if a.mask:
    MSK = np.clip(sperm_inv(crop_real(rd_mrc(a.mask)), PERM, SGN), 0, 1) * msk3; log('focus mask: %.1f %% of the c^3 sphere above 0.5' % (100 * (MSK > 0.5).sum() / inm.sum()))
def write_maps(tag):
    if a.nomaps: return
    od = '%s/maps_%s' % (a.out, tag); os.makedirs(od, exist_ok=True)
    for g in range(K):
        V = 0.5 * (wmap(D[0], D[0]['W'][:, g]) + wmap(D[1], D[1]['W'][:, g])); wr_mrc('%s/state_%02d.mrc' % (od, g + 1), sperm(pad_real(V, a.box), PERM, SGN), a.smpd)
def maps():   # consensus and state maps of both halves, once per iteration
    return [(wmap(D[hh], np.ones(len(D[hh]['W']), np.float32)), [wmap(D[hh], D[hh]['W'][:, g]) for g in range(K)]) for hh in (0, 1)]
def refs(h, out):   # references for half h (from the OTHER half): consensus + per-shell FSC-weighted deviation, the deviation inside the mask
    (cA, MA), (cB, MB) = out; R = []; fs = []
    for g in range(K):
        DA, DB = MA[g][:-1] - cA[:-1], MB[g][:-1] - cB[:-1]
        f = np.bincount(sh, np.real(DA * np.conj(DB)), ns) / np.maximum(np.sqrt(np.bincount(sh, np.abs(DA) ** 2, ns) * np.bincount(sh, np.abs(DB) ** 2, ns)), 1e-30); f = np.clip(f, 0., 1.); fs.append(f)
        Dg = f[sh] * (DB if h == 0 else DA)
        if MSK is not None: Dg = to_four(MSK * to_real(np.r_[Dg, 0.]))[:-1]
        R.append(np.r_[(cB if h == 0 else cA)[:-1] + Dg, 0.].astype(np.complex64))
    return R, np.array(fs)
def estep(h, R):
    d = D[h]; nn = len(d['W']); mp = d['W'].argmax(1); bn = np.zeros(ns); bd = np.zeros(ns); R2 = np.zeros(ns); C2 = np.zeros(ns)
    for g in range(K):
        rr = np.where(mp == g)[0]
        for b in range(0, len(rr), 8000):
            r = rr[b:b + 8000]; Pg = R[g][d['ix'][r]]; q = d['shq'][r].ravel()
            bn += np.bincount(q, np.real(Pg * np.conj(d['Wv'][r])).ravel(), ns); bd += np.bincount(q, (d['Wd'][r] * np.abs(Pg) ** 2).ravel(), ns)
    beta = np.clip(bn / np.maximum(bd, 1e-30), 0., 2.)
    for g in range(K):
        rr = np.where(mp == g)[0]
        for b in range(0, len(rr), 8000):
            r = rr[b:b + 8000]; Pg = R[g][d['ix'][r]] * beta[d['shq'][r]]; res = d['F2'][r] - 2. * np.real(Pg * np.conj(d['Wv'][r])) + d['Wd'][r] * np.abs(Pg) ** 2
            R2 += np.bincount(d['shq'][r].ravel(), res.ravel(), ns); C2 += np.bincount(d['shq'][r].ravel(), minlength=ns)
    wsh = (1. / np.maximum(R2 / np.maximum(C2, 1), 1e-12)).astype(np.float32); S = np.zeros((nn, K), np.float32)
    for b in range(0, nn, 8000):
        s_ = slice(b, b + 8000); q = d['shq'][s_]; bw = beta[q].astype(np.float32); wi = wsh[q]
        for g in range(K):
            Pg = R[g][d['ix'][s_]] * bw; S[s_, g] = 0.5 * ((2. * np.real(Pg * np.conj(d['Wv'][s_])) - d['Wd'][s_] * np.abs(Pg) ** 2) * wi).sum(1)
    return S, beta
def report(tag, fs=None):
    occ = [D[h]['W'].sum(0) for h in (0, 1)]; s = '%s: occupancy A %s | B %s' % (tag, ' '.join('%.0f' % v for v in occ[0]), ' '.join('%.0f' % v for v in occ[1]))
    if fs is not None: s += ' | deviation FSC 30-15 A per state %s' % ' '.join('%.2f' % v for v in fs[:, int(round(n * a.smpd / 30)):int(round(n * a.smpd / 15)) + 1].mean(1))
    if a.gt:
        pur = [sum((D[h]['W'][:, g] * D[h]['closed']).sum() for h in (0, 1)) / max(sum(D[h]['W'][:, g].sum() for h in (0, 1)), 1e-9) for g in range(K)]
        ncl = sum(D[h]['closed'].sum() for h in (0, 1)); g1 = int(np.argmax(pur)); rec = sum((D[h]['W'][:, g1] * D[h]['closed']).sum() for h in (0, 1)) / max(ncl, 1)
        s += ' | closed share per state %s; most closed state %d: purity %.2f, recall %.2f' % (' '.join('%.2f' % p for p in pur), g1 + 1, pur[g1], rec)
    log(s)
report('start'); write_maps('start')
for it in range(a.iters):
    Mi = maps(); RA, fs = refs(0, Mi); RB, _ = refs(1, Mi); ch = []
    if a.empty_init and it == 0:   # class 1 = the other half's consensus, class 2 = it with the site removed
        for R_, (cons_o, _) in ((RA, Mi[1]), (RB, Mi[0])):
            c_ = np.r_[cons_o[:-1], 0.]; R_[0] = c_.astype(np.complex64); R_[1] = (c_ - to_four(MSK * to_real(c_))).astype(np.complex64)
    for h, R in ((0, RA), (1, RB)):
        d = D[h]; S, beta = estep(h, R); pi = np.maximum(d['W'].mean(0), a.pifloor)
        if a.empty_init and it == 0: pi = np.array([0.9, 0.1])
        S += np.log(pi)[None]; old = d['W'].argmax(1)
        if a.soft: S -= S.max(1, keepdims=True); W = np.exp(S); d['W'] = (W / W.sum(1, keepdims=True)).astype(np.float32)
        else: W = np.zeros_like(d['W']); W[np.arange(len(W)), S.argmax(1)] = 1.; d['W'] = W
        ch.append((d['W'].argmax(1) != old).mean())
        if it == 0: log('  half %s reference scale per shell (2..%d): %s' % ('AB'[h], kb, np.array2string(beta[2:kb + 1], precision=2)))
    report('iteration %d (%s; argmax changed A %.3f B %.3f; median max responsibility A %.2f)' % (it + 1, mode, ch[0], ch[1], float(np.median(D[0]['W'].max(1)))), fs)
write_maps(mode)
for h in (0, 1):
    out = np.zeros(len(LABS[h]), int); out[D[h]['rows']] = D[h]['W'].argmax(1) + 1; np.savetxt('%s/labels_%s_%s.txt' % (a.out, mode, 'AB'[h]), out, fmt='%d')
    if a.post: P_ = np.zeros((len(LABS[h]), K), np.float32); P_[D[h]['rows']] = D[h]['W']; np.savetxt('%s/post_%s_%s.txt' % (a.out, mode, 'AB'[h]), P_, fmt='%.4f')
log('DONE')
