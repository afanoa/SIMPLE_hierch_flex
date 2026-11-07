#!/usr/bin/env python3
"""Junk component prototype (synthetic sets with truth): soft EM at fixed poses over a node's states plus junk components, each half
scored against the OTHER half's references (as splitbench13's cross_em). States: reference = half consensus + deviation weighted by the
deviation's cross-half FSC. Junk components: NOISE (no map, score 0: a particle no state explains better than noise) and FREE components
(klem's free components: their own reconstruction from their members, low-passed to LPJ A, never shrunk to the consensus), initialised
from the particles with the lowest own-state score. Per iteration: junk precision / recall against truth (conformer 0), NMI(conformer;
state) of the non-junk calls. --hide 1: truth-junk particles start in random real states (does the junk component pull them out?).
Geometry and image model copied from splitbench13.py (setup / rotm / load_half / rec / shell_eps).
v2: --free_init klem = klem's junk-like start (free 1 = the consensus low-passed to 35 A, free 2 = the consensus at half contrast, 2 % of the
mass), free maps shrunk toward the consensus (--free_shrink cons, klem's M-step) and scored over the wider band (--em_lim 12); the junk is
read out per type (blob / broken, typed by the CTF triple from SIM/oris_{blob,broken}.txt).
usage: junk_em2.py --sim S --poses P --eo EO --labels L --gt GT --cache C --out DIR [--free 2 --noise 0 --free_init klem --lpj 0 --iters 8 --hide 1]"""
import os, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--sim', required=True); ap.add_argument('--poses', required=True); ap.add_argument('--eo', required=True)
ap.add_argument('--labels', required=True); ap.add_argument('--gt', required=True); ap.add_argument('--cache', default=''); ap.add_argument('--out', required=True)
ap.add_argument('--free', type=int, default=2); ap.add_argument('--noise', type=int, default=0); ap.add_argument('--lpj', type=float, default=0.); ap.add_argument('--iters', type=int, default=8)
ap.add_argument('--hide', type=int, default=1); ap.add_argument('--init_frac', type=float, default=0.05); ap.add_argument('--em_lim', type=float, default=12.); ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--free_init', default='klem'); ap.add_argument('--free_shrink', default='cons'); ap.add_argument('--blob_lp', type=float, default=35.)
ap.add_argument('--crop', type=int, default=48); ap.add_argument('--mskdiam', type=float, default=230.); ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--kv', type=float, default=300.)
ap.add_argument('--cs', type=float, default=2.7); ap.add_argument('--fraca', type=float, default=0.1); ap.add_argument('--conv', default='1,1,1,1,1,1'); ap.add_argument('--ctfabs', type=int, default=0); ap.add_argument('--eps', type=float, default=1.)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); rng = np.random.default_rng(a.seed)
cv = [float(v) for v in a.conv.split(',')]; convT, s1, yf, s3, ssh, sctf = bool(cv[0]), cv[1], cv[2], cv[3], cv[4], cv[5]
def log(*s): print(*s, flush=True)
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
def setup(n):
    c = a.crop; x = np.arange(n) - n // 2; rr2 = np.hypot(x[:, None], x[None, :]); R = a.mskdiam / 2 / a.smpd
    msk2 = np.clip((R + 4 - rr2) / 8., 0, 1).astype(np.float32); sl = np.r_[0:c // 2, n - c // 2:n]
    fc = np.fft.fftfreq(c) * c / n; fy2, fx2 = np.meshgrid(fc, fc, indexing='ij'); kk0 = np.stack([fx2, fy2, np.zeros_like(fx2)], -1).reshape(-1, 3) * n
    f3 = np.fft.fftfreq(c) * c / n; kx, ky, kz = np.meshgrid(f3, f3, f3, indexing='ij'); sh3 = np.minimum(np.round(np.sqrt(kx**2 + ky**2 + kz**2) * n).astype(int), c // 2).ravel()
    g = np.arange(c) - c // 2; r3 = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2) * n / c; msk3 = np.clip((R + 4 - r3) / 8., 0, 1).astype(np.float32)
    return dict(n=n, c=c, msk2=msk2, sl=sl, fx2=fx2, fy2=fy2, kk0=kk0, sh3=sh3, msk3=msk3)
lam = 12.2643247 / np.sqrt(a.kv * 1e3 * (1. + a.kv * 1e3 * 0.978466e-6))
def rotm(e1, e2, e3):
    a1, a2, a3 = np.radians(e1), np.radians(e2), np.radians(e3); c1, si1, c2, si2, c3, si3 = np.cos(a1), np.sin(a1), np.cos(a2), np.sin(a2), np.cos(a3), np.sin(a3)
    Z1 = np.zeros((len(e1), 3, 3)); Z1[:, 0, 0] = c1; Z1[:, 0, 1] = -si1; Z1[:, 1, 0] = si1; Z1[:, 1, 1] = c1; Z1[:, 2, 2] = 1
    Y2 = np.zeros((len(e1), 3, 3)); Y2[:, 0, 0] = c2; Y2[:, 0, 2] = si2; Y2[:, 2, 0] = -si2; Y2[:, 2, 2] = c2; Y2[:, 1, 1] = 1
    Z3 = np.zeros((len(e1), 3, 3)); Z3[:, 0, 0] = c3; Z3[:, 0, 1] = -si3; Z3[:, 1, 0] = si3; Z3[:, 1, 1] = c3; Z3[:, 2, 2] = 1
    return Z1 @ Y2 @ Z3
def load_half(simdir, posef, lab, tag):
    rows = np.where(lab > 0)[0]
    P = np.array([[float(kv(l).get(k, 0)) for k in ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast')] for l in open(posef) if l.strip()])
    h = np.frombuffer(open(simdir + '/all.mrcs', 'rb').read(1024), np.int32); nx, ny, nz, nsym = [int(v) for v in (h[0], h[1], h[2], h[23])]
    I = np.memmap(simdir + '/all.mrcs', np.float32, 'r', offset=1024 + nsym, shape=(nz, ny, nx)); G = setup(nx); c = G['c']
    cache = '%s_%s_c%d.npy' % (a.cache, tag, c) if a.cache else ''
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
    log('%s: %d particles, crop %d (Nyquist %.1f A)' % (tag, len(rows), c, nx * a.smpd / (c // 2)))
    return dict(rows=rows, Wv=Wv, Wd=Wd, F2=F2, idx=idx, G=G)
lab0 = np.loadtxt(a.labels, dtype=int); eo = np.loadtxt(a.eo, dtype=int); conf = np.full(len(lab0), -1)
GT = np.loadtxt(a.gt, comments='#'); gi = GT[:, 0].astype(int) - 1; ok = gi < len(lab0); conf[gi[ok]] = GT[ok, 1].astype(int)
K = int(lab0.max()); start = lab0.copy(); tj = conf == 0
jtype = np.zeros(len(lab0), int)   # 0 = conformer, 1 = blob, 2 = broken, 3 = junk of unknown type
def ckey(l): d = kv(l); return (round(float(d['dfx']), 5), round(float(d['dfy']), 5), round(float(d['angast']), 3))
if os.path.isfile(a.sim + '/oris_blob.txt') and os.path.isfile(a.sim + '/all.txt'):
    tmap = {ckey(l): t for t, f in ((1, 'blob'), (2, 'broken')) for l in open('%s/oris_%s.txt' % (a.sim, f)) if l.strip()}
    for i, l in enumerate(l for l in open(a.sim + '/all.txt') if l.strip()):
        if i < len(jtype) and tj[i]: jtype[i] = tmap.get(ckey(l), 3)
    log('junk types: blob %d, broken %d, unknown %d' % ((jtype == 1).sum(), (jtype == 2).sum(), (jtype == 3).sum()))
if a.hide:   # truth-junk particles moved into random real states (states whose members are mostly not junk)
    real = [s for s in range(1, K + 1) if (lab0 == s).any() and tj[lab0 == s].mean() < 0.5]; m = tj & (lab0 > 0); start[m] = rng.choice(real, m.sum())
    log('hidden: %d truth-junk particles moved into the %d real states; %d states were junk-majority (now empty of their junk)' % (m.sum(), len(real), K - len(real)))
HA = load_half(a.sim, a.poses, np.where(eo == 0, start, 0), 'A'); HB = load_half(a.sim, a.poses, np.where(eo == 1, start, 0), 'B')
G = HA['G']; c = G['c']; n = G['n']; sh = G['sh3']; ns = c // 2 + 1; m3 = c ** 3 + 1
kb = min(c // 2 - 1, int(round(n * a.smpd / a.em_lim))); kj = min(c // 2 - 1, int(round(n * a.smpd / a.lpj))) if a.lpj > 0 else c // 2
r2 = np.round(np.hypot(G['kk0'][:, 0], G['kk0'][:, 1])).astype(int); cols = np.where((r2 >= 2) & (r2 <= kb))[0]
NJF, NJN = a.free, a.noise; NC = K + NJF + NJN   # components: states 1..K, free 1..NJF, noise
P = {}
for h, (H, lb) in enumerate(((HA, np.where(eo == 0, start, 0)), (HB, np.where(eo == 1, start, 0)))):
    st = lb[H['rows']]; ix = np.ascontiguousarray(H['idx'][:, cols]); W = np.zeros((len(st), NC), np.float32); W[np.arange(len(st)), st - 1] = 1.
    P[h] = dict(Wv=np.ascontiguousarray(H['Wv'][:, cols]), Wd=np.ascontiguousarray(H['Wd'][:, cols]), F2=np.ascontiguousarray(H['F2'][:, cols]), ix=ix,
                shq=np.minimum(sh[ix % (c ** 3)], ns - 1), W=W, conf=conf[H['rows']], tj=tj[H['rows']], jt=jtype[H['rows']])
neq = int(np.median([max((P[h]['W'][:, g] > 0).sum(), 1) for h in (0, 1) for g in range(K)]))
def shell_eps(d):
    den = np.bincount(d['ix'].ravel(), d['Wd'].ravel(), m3)[:-1] * (neq / len(d['W'])); return np.array([np.median(den[(sh == k) & (den > 0)]) if ((sh == k) & (den > 0)).any() else 1. for k in range(ns)]) * a.eps
eps = 0.5 * (shell_eps(P[0]) + shell_eps(P[1])); epsv = np.r_[eps[sh], eps[-1]]
def wmap(d, w, tgt=None):
    ii = d['ix'].ravel(); num = np.bincount(ii, (d['Wv'].real * w[:, None]).ravel(), m3) + 1j * np.bincount(ii, (d['Wv'].imag * w[:, None]).ravel(), m3)
    e = epsv * max(float(w.sum()), 1.) / neq; return (num + (0. if tgt is None else e * tgt)) / (np.bincount(ii, (d['Wd'] * w[:, None]).ravel(), m3) + e)
def lowpass(V, k): V = V.copy(); V[:-1][sh > k] = 0.; V[-1] = 0.; return V
def maps(first=False):   # consensus and component maps of both halves (once per iteration); free maps shrunk toward the consensus (klem) or zero
    out = []; kb35 = int(round(n * a.smpd / a.blob_lp))
    for hh in (0, 1):
        d = P[hh]; cons = wmap(d, np.ones(len(d['W']), np.float32))
        M = [wmap(d, d['W'][:, g], cons if (g >= K and a.free_shrink == 'cons') else None) for g in range(K + NJF)]
        if first and a.free_init == 'klem':   # klem's junk-like start: the consensus low-passed to 35 A, then the consensus at half contrast
            for j in range(NJF): M[K + j] = lowpass(cons, kb35) if j % 2 == 0 else 0.5 * cons
        out.append((cons, M))
    return out
def refs(h, out):   # references half h is scored against (the OTHER half's): states = consensus + deviation weighted by its cross-half FSC, free = own map low-passed
    (cA, MA), (cB, MB) = out; R = []
    for g in range(K):
        DA, DB = MA[g][:-1] - cA[:-1], MB[g][:-1] - cB[:-1]
        f = np.bincount(sh, np.real(DA * np.conj(DB)), ns) / np.maximum(np.sqrt(np.bincount(sh, np.abs(DA) ** 2, ns) * np.bincount(sh, np.abs(DB) ** 2, ns)), 1e-30)
        f = np.clip(f, 0., 1.); R.append(np.r_[(cB if h == 0 else cA)[:-1] + f[sh] * (DB if h == 0 else DA), 0.].astype(np.complex64))
    for g in range(K, K + NJF): R.append((lowpass((MB if h == 0 else MA)[g], kj) if a.lpj > 0 else np.r_[(MB if h == 0 else MA)[g][:-1], 0.]).astype(np.complex64))
    return R
def estep(h, R):
    d = P[h]; nn = len(d['W']); mp = d['W'].argmax(1); R2 = np.zeros(ns); C2 = np.zeros(ns)
    for g in range(K):   # per-shell noise from the residuals of the particles currently in real states
        rows = np.where(mp == g)[0]
        for b in range(0, len(rows), 8000):
            r = rows[b:b + 8000]; Pg = R[g][d['ix'][r]]; res = d['F2'][r] - 2. * np.real(Pg * np.conj(d['Wv'][r])) + d['Wd'][r] * np.abs(Pg) ** 2
            R2 += np.bincount(d['shq'][r].ravel(), res.ravel(), ns); C2 += np.bincount(d['shq'][r].ravel(), minlength=ns)
    winv = (1. / np.maximum(R2 / np.maximum(C2, 1), 1e-12)).astype(np.float32)[d['shq']]; S = np.zeros((nn, NC), np.float32)
    for g in range(K + NJF):
        for b in range(0, nn, 8000):
            s_ = slice(b, b + 8000); Pg = R[g][d['ix'][s_]]; S[s_, g] = 0.5 * ((2. * np.real(Pg * np.conj(d['Wv'][s_])) - d['Wd'][s_] * np.abs(Pg) ** 2) * winv[s_]).sum(1)
    # the noise component's score is 0 (column K + NJF stays zero)
    return S
def nmi(x, y):
    ux, ix = np.unique(x, return_inverse=True); uy, iy = np.unique(y, return_inverse=True); C = np.zeros((len(ux), len(uy))); np.add.at(C, (ix, iy), 1); Pm = C / C.sum()
    px, py = Pm.sum(1), Pm.sum(0); nz = Pm > 0; I = (Pm[nz] * np.log(Pm[nz] / np.outer(px, py)[nz])).sum(); return I / max(np.sqrt((px * np.log(px)).sum() * (py * np.log(py)).sum()), 1e-12)
def report(tag):
    s = []
    for h in (0, 1):
        d = P[h]; mp = d['W'].argmax(1); jc = mp >= K; t = d['tj']; tp = (jc & t).sum()
        prec = tp / max(jc.sum(), 1); rec_ = tp / max(t.sum(), 1); nj = ~jc & ~t
        j8 = d['W'][:, K:].sum(1) > 0.8; tp8 = (j8 & t).sum(); ty = ' (blob %.2f, broken %.2f)' % (jc[d['jt'] == 1].mean() if (d['jt'] == 1).any() else 0., jc[d['jt'] == 2].mean() if (d['jt'] == 2).any() else 0.)
        s.append('%s: junk calls %d (precision %.2f, recall %.2f%s; resp>0.8: %d, %.2f, %.2f), free %s noise %d, NMI(conformer; state) of the rest %.3f' % ('AB'[h], jc.sum(), prec, rec_, ty,
                 j8.sum(), tp8 / max(j8.sum(), 1), tp8 / max(t.sum(), 1), '/'.join(str((mp == K + j).sum()) for j in range(NJF)) or '-', (mp == K + NJF).sum() if NJN else 0,
                 nmi(d['conf'][nj], mp[nj]) if nj.sum() > 100 else float('nan')))
    log('  %s | %s' % (tag, ' | '.join(s)))
# free components start from the particles with the lowest own-state score (FRAC of each half), split at random among the free components (--free_init low),
# or as klem's junk-like maps holding 2 % of every particle's mass (--free_init klem)
M0 = maps(first=True); R0 = {0: refs(0, M0), 1: refs(1, M0)}
for h in (0, 1):
    if a.free_init == 'klem' and NJF:
        d = P[h]; d['W'][:, :K] *= 0.98; d['W'][:, K:K + NJF] = 0.02 / NJF; continue
    d = P[h]; S = estep(h, R0[h]); own = S[np.arange(len(d['W'])), d['W'].argmax(1)]; lo = np.argsort(own)[:int(a.init_frac * len(own))]
    if NJF: d['W'][lo] = 0.; d['W'][lo, K + rng.integers(NJF, size=len(lo))] = 1.
    log('half %s: %d lowest-scoring particles start the free components (truth junk among them %.2f; junk in the half %.2f)' % ('AB'[h], len(lo), d['tj'][lo].mean(), d['tj'].mean()))
report('start')
for it in range(a.iters):
    Mi = maps(); RA, RB = refs(0, Mi), refs(1, Mi)
    for h, R in ((0, RA), (1, RB)):
        d = P[h]; S = estep(h, R); pi = np.maximum(d['W'].mean(0), 1e-6); pi[K + NJF:] = max(pi[K + NJF:].max() if NJN else 1e-6, 0.01)   # the noise component keeps a floor weight
        S += np.log(pi)[None]; S -= S.max(1, keepdims=True); W = np.exp(S); d['W'] = (W / W.sum(1, keepdims=True)).astype(np.float32)
    report('iteration %d' % (it + 1))
for h, side in ((0, 'A'), (1, 'B')):
    out = np.zeros(len(lab0), int); H = (HA, HB)[h]; out[H['rows']] = P[h]['W'].argmax(1) + 1; np.savetxt('%s/calls_%s.txt' % (a.out, side), out, fmt='%d')
log('components: states 1..%d, free %s, noise %s' % (K, ','.join(str(K + 1 + j) for j in range(NJF)) or '-', str(K + NJF + 1) if NJN else '-'))
mpall = np.r_[P[0]['W'].argmax(1), P[1]['W'].argmax(1)]; tjall = np.r_[P[0]['tj'], P[1]['tj']]
log('final components (both halves; members, truth-junk fraction): ' + '  '.join('%s%d %d %.2f' % ('s' if g < K else ('f' if g < K + NJF else 'n'), g + 1, (mpall == g).sum(),
    tjall[mpall == g].mean() if (mpall == g).any() else 0.) for g in range(NC)))
