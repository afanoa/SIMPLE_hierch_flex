#!/usr/bin/env python3
"""latent_states.py: fixed-pose test of hetero_tree node states INSIDE the node's flex latent space, against the production
stage-B labels.

At the node's input poses each particle image is scored against the node mean mu (flex_pca's vol1) and a basis W = flex_pca's
modes (fflex/flex_pca_pcNNN.mrc), optionally plus SVD site factors of the node's recB0 half maps inside its saved site masks.
Over the band to --lp (the tree's flex lp, 15 A), with site_em's per-shell scale and noise and a per-particle amplitude a_i:
  b_i = (A_i W)' S^-1 (y_i - A_i mu),   G_i = (A_i W)' S^-1 (A_i W)    (A_i = a_i x scale x CTF x central slice at the input pose)
States are a Gaussian mixture in the latent space, z ~ sum_c pi_c N(xi_c, Omega) with Omega tied across states, fitted by EM with
z marginalised (extreme deconvolution on (b_i, G_i)): soft posteriors and their argmax.
Arms: P = production stage-B labels, R0 = the provisional states after stA0 (recB0), F = flex's delivered states (fflex), L = latent
states on the flex modes, LS = latent states on the flex modes + site factors; L and LS at the production K and at their own K
(cross-half held-out likelihood). Readouts: truth NMI / ARI, soft calibration (mean member posterior vs purity), quality association,
closed reference rows, and the latent held-out gain over one Gaussian of every arm's labels.
usage: latent_states.py --node NODE --dump DUMP --out OUT --cache planes.npy --truth KIND:ARGS [--lp 15 --nproc 8 --kmax 8]
  KIND: qtoy:<truth_qtoy.txt> | igmfix:<truth.txt> | igmnull:<truth.txt> | closed:<fixture_rows.txt>:<closed_rows.txt> | none
  --selftest: synthetic check of the latent EM only"""
import argparse, glob, json, os, re, sys, time
import numpy as np
from multiprocessing import Pool

PP = os.environ.get('PPINFO', '/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo')
sys.path.insert(0, PP)

MLP = 15.            # read_lp: Gaussian low-pass of the half maps for site discovery (half maximum at MLP)
CH = 400             # particles per sufficient-statistics chunk
G = {}               # arrays shared with the forked workers


def log(*a):
    print('[%s] %s' % (time.strftime('%H:%M:%S'), ' '.join(str(x) for x in a)), flush=True)


# ---------------------------------------------------------------- node, planes, geometry

def node_select(n0, node):
    """dump indices and project rows of the node: a local node project holds exactly rows.txt's rows, a full one is indexed by them."""
    rows = np.loadtxt(os.path.join(node, 'rows.txt'), dtype=int, ndmin=1)
    if n0 == len(rows):
        return np.arange(n0), rows
    if rows.max() <= n0:
        return rows - 1, rows
    raise SystemExit('rows.txt (%d rows, max %d) does not fit the dump (%d rows)' % (len(rows), rows.max(), n0))


def flex_keys(node):
    sh = open(os.path.join(node, 'fflex', 'hetero_step.sh')).read()
    vol1 = re.search(r'vol1=(\S+)', sh).group(1)
    md = re.search(r'mskdiam=(\S+)', sh)
    return vol1, float(md.group(1)) if md else 0.


def setup(a):
    from ppinfo_site import read_dump, read_stk, euler2m, prep_crop, ctf_eval
    nd = read_dump(os.path.join(a.dump, 'n_ptcl3D.txt'))
    recs, box, smpd, kv, cs, fraca = read_stk(os.path.join(a.dump, 'n_stk.txt'))
    n0 = len(nd['e1'])
    didx, prow = node_select(n0, a.node)
    ok = (nd['eo'][didx] == 0) | (nd['eo'][didx] == 1)
    didx, prow = didx[ok], prow[ok]
    if a.maxn and a.maxn < len(didx):
        s = np.sort(np.random.default_rng(3).choice(len(didx), a.maxn, replace=False))
        didx, prow = didx[s], prow[s]
    n = len(didx)
    vol1, mskdiam = flex_keys(a.node)
    mskdiam = a.mskdiam if a.mskdiam > 0 else mskdiam
    box_A = box * smpd
    kb = int(round(box_A / a.lp))
    B = min(box, max(26, 2 * kb + 4)); B += B % 2
    nyq, nh, khi = B // 2, B // 2 + 1, B // 2 - 1
    kb = min(kb, khi)
    log('node %s: %d particles (dump %d rows), box %d smpd %.3f, mskdiam %.1f -> box_crop %d, band shells 2..%d (%.1f A), mean %s'
        % (a.node, n, n0, box, smpd, mskdiam, B, kb, box_A / kb, vol1))
    e = np.c_[nd['e1'], nd['e2'], nd['e3'], nd['x'], nd['y']][didx]
    R = euler2m(e[:, :3]).astype(np.float32)
    indstk = nd['indstk'][didx].astype(int)
    stkind = np.maximum(nd['stkind'][didx].astype(int), 1)
    crop = prep_crop(recs, box, stkind, indstk, mskdiam / 2. / smpd, B, a.cache)
    kk, hh = np.meshgrid(np.r_[0:B // 2, -B // 2:0], np.arange(B // 2), indexing='ij')
    r2 = hh * hh + kk * kk
    sel = (r2 >= 3) & (r2 <= kb * (kb + 1))
    h, k = hh[sel].astype(np.int32), kk[sel].astype(np.int32)
    krow = np.where(k >= 0, k, k + B)
    sw = np.where(h == 0, 1., 2.)
    shell = np.rint(np.sqrt(h * h + k * k)).astype(int) - 1
    nsh = shell.max() + 1
    X = crop[:, krow, h]
    del crop
    X *= np.exp(-2j * np.pi * (h[None, :] * e[:, 3:4] + k[None, :] * e[:, 4:5]) / box).astype(np.complex64)
    ctf = ctf_eval(h, k, nd['dfx'][didx], nd['dfy'][didx], nd['angast'][didx], nd['phshift'][didx], kv, cs, fraca, box_A)
    f2 = np.zeros((n, nsh))
    for s in range(nsh):
        m = shell == s
        f2[:, s] = (sw[m][None, :] * (X[:, m].real ** 2 + X[:, m].imag ** 2)).sum(1)
    X *= ctf
    NC = B * B * nh
    idx = np.empty((n, len(h)), np.int32)
    cnj = np.empty((n, len(h)), bool)
    for b0 in range(0, n, 4000):
        b1 = min(n, b0 + 4000)
        p = np.rint(h[None, :, None] * R[b0:b1, None, 0, :] + k[None, :, None] * R[b0:b1, None, 1, :]).astype(np.int32)
        c = p[..., 0] < 0
        p[c] = -p[c]
        okk = (np.abs(p) < nyq).all(-1)
        idx[b0:b1] = np.where(okk, (p[..., 2] % B) * B * nh + (p[..., 1] % B) * nh + p[..., 0], NC)
        cnj[b0:b1] = c
    l = np.arange(B); ls = np.where(l < nyq, l, l - B); hc = np.arange(nh)
    cell_shell = np.minimum(np.rint(np.sqrt(ls[:, None, None] ** 2 + ls[None, :, None] ** 2 + hc[None, None, :] ** 2)), nyq).astype(np.int32).ravel()
    cell_fw = np.where((hc > 0) & (hc < nh - 1), 2., 1.)[None, None, :].repeat(B, 0).repeat(B, 1).ravel()
    eo = nd['eo'][didx].astype(int)
    G.update(dict(P=X, C=ctf * ctf, f2=f2, idx=idx, cnj=cnj, sw=sw, shell=shell, nsh=nsh, NC=NC, B=B, nh=nh, nyq=nyq, kb=kb, box=box,
                  smpd=smpd, box_A=box_A, mskdiam=mskdiam, cell_shell=cell_shell, cell_fw=cell_fw, eo=eo, n=n, n0=n0, didx=didx,
                  prow=prow, vol1=vol1))


def vol_ft(v):
    """real volume (any box) -> rfftn half volume at box_crop, raveled, with a trailing zero for out-of-band samples."""
    from mrcio import fcrop3
    if v.shape[0] != G['B']:
        v = fcrop3(v, G['B'])
    V = np.fft.rfftn(np.fft.ifftshift(v)).astype(np.complex128).ravel()
    return np.r_[V, 0.]


def slices(V, sl):
    s = V[G['idx'][sl]]
    return np.where(G['cnj'][sl], np.conj(s), s)


# ---------------------------------------------------------------- the basis

def lp_gauss(v, smpd, lp):
    n = v.shape[0]
    f = np.fft.fftfreq(n)
    f2 = f[:, None, None] ** 2 + f[None, :, None] ** 2 + f[None, None, :] ** 2
    return np.real(np.fft.ifftn(np.fft.fftn(v) * np.exp(-f2 / smpd ** 2 * lp ** 2 * np.log(2.))))


def site_factors(node, lab_r0):
    """SVD factors of the recB0 states' half-map deviations inside each saved site mask (the tree's site discovery input): per site
    the factors of the averaged halves' weighted deviations whose even and odd projections on the factor's state direction correlate
    >= 0.5 (at most two per site)."""
    from mrcio import read_mrc
    sites = sorted(glob.glob(os.path.join(node, 'sites', 'site_[0-9][0-9].mrc')))
    r0 = os.path.join(node, 'recB0')
    st = sorted(int(re.search(r'recvol_state(\d\d)_even_unfil', f).group(1))
                for f in glob.glob(os.path.join(r0, 'recvol_state??_even_unfil.mrc')))
    st = [s for s in st if (lab_r0 == s).sum() > 0 and os.path.exists(os.path.join(r0, 'recvol_state%02d_odd_unfil.mrc' % s))]
    if len(sites) == 0 or len(st) < 2:
        log('site factors: %d sites, %d recB0 states with maps -> none' % (len(sites), len(st)))
        return [], []
    w = np.array([(lab_r0 == s).sum() for s in st], float); w /= w.sum()
    halves = []
    for half in ('even', 'odd'):
        V = []
        for s in st:
            v, smpd = read_mrc(os.path.join(r0, 'recvol_state%02d_%s_unfil.mrc' % (s, half)))
            V.append(lp_gauss(v, smpd, MLP))
        V = np.array(V)
        n = V.shape[1]
        x = np.arange(n) - n // 2
        sph = (x[:, None, None] ** 2 + x[None, :, None] ** 2 + x[None, None, :] ** 2) <= (G['mskdiam'] / 2. / smpd) ** 2
        em = np.tensordot(w, V, 1)
        sc = np.array([(v[sph] * em[sph]).sum() / max((v[sph] ** 2).sum(), 1e-30) for v in V])
        halves.append(np.sqrt(w)[:, None, None, None] * (sc[:, None, None, None] * V - em[None]))
    de, do = halves
    maps, info = [], []
    for f in sites:
        ms, _ = read_mrc(f)
        ms = np.clip(ms, 0., 1.)
        vx = ms > 0.01
        Ae = (de[:, vx] * ms[vx][None]).T
        Ao = (do[:, vx] * ms[vx][None]).T
        U, S, Vt = np.linalg.svd(0.5 * (Ae + Ao), full_matrices=False)
        for j in range(min(2, len(S))):
            pe, po = Ae @ Vt[j], Ao @ Vt[j]
            r = float(np.corrcoef(pe, po)[0, 1])
            share = float(S[j] ** 2 / max((S ** 2).sum(), 1e-30))
            keep = r >= 0.5
            info.append(dict(site=os.path.basename(f), factor=j + 1, r_eo=r, share=share, kept=bool(keep)))
            if keep:
                m = np.zeros(ms.shape); m[vx] = U[:, j]
                maps.append(m)
        log('site %s: %d voxels, factors r_eo %s, shares %s' % (os.path.basename(f), int(vx.sum()),
            ' '.join('%.2f' % d['r_eo'] for d in info if d['site'] == os.path.basename(f)),
            ' '.join('%.2f' % d['share'] for d in info if d['site'] == os.path.basename(f))))
    return maps, info


def orthonormal(Vs, names):
    """orthonormalise the basis in the band (Fourier inner product over shells <= kb); the mixture is invariant to the basis's
    linear frame, so this only conditions G. Returns the new maps and, per input map, the share of its band energy outside the
    span of the maps before it (the flex modes come first)."""
    band = np.r_[(G['cell_shell'] <= G['kb']), False]
    fw = np.r_[G['cell_fw'], 0.] * band
    Vs = [V / np.sqrt(max(np.sum(fw * np.abs(V) ** 2), 1e-300)) for V in Vs]
    M = np.array([[np.sum(fw * (a.real * b.real + a.imag * b.imag)) for b in Vs] for a in Vs])
    resid = []
    for j in range(len(Vs)):
        if j == 0:
            resid.append(1.)
            continue
        A = M[:j, :j]; bj = M[:j, j]
        c = np.linalg.lstsq(A, bj, rcond=None)[0]
        resid.append(float(max(M[j, j] - bj @ c, 0.) / max(M[j, j], 1e-30)))
    lam, U = np.linalg.eigh(M)
    keep = lam > 1e-6 * lam.max()
    U, lam = U[:, keep], lam[keep]
    T = U / np.sqrt(lam)[None]
    W = [sum(T[j, l] * Vs[j] for j in range(len(Vs))) for l in range(T.shape[1])]
    log('basis: %d maps (%s) -> rank %d; band energy outside the preceding span: %s'
        % (len(Vs), ' '.join(names), len(W), ' '.join('%s %.2f' % (nm, r) for nm, r in zip(names, resid))))
    return W, resid


# ---------------------------------------------------------------- scale, noise, sufficient statistics

def ref_terms(V0, amp):
    n, nsh, shell, sw = G['n'], G['nsh'], G['shell'], G['sw']
    aa = np.zeros((n, nsh)); bb = np.zeros((n, nsh))
    for c0 in range(0, n, CH):
        sl = slice(c0, min(n, c0 + CH))
        s0 = slices(V0, sl)
        a = (s0.real * G['P'][sl].real + s0.imag * G['P'][sl].imag) * sw[None]
        b = (s0.real ** 2 + s0.imag ** 2) * G['C'][sl] * sw[None]
        nc = sl.stop - sl.start
        key = (np.arange(nc)[:, None] * nsh + shell[None, :]).ravel()
        aa[sl] = np.bincount(key, a.ravel(), nc * nsh).reshape(nc, nsh)
        bb[sl] = np.bincount(key, b.ravel(), nc * nsh).reshape(nc, nsh)
    return aa, bb


def scale_noise(V0):
    """per-shell model scale beta and noise sig2 (site_em's estimators) with a per-particle amplitude a_i (median 1)."""
    n = G['n']
    cw = np.bincount(G['shell'], G['sw'], G['nsh'])
    amp = np.ones(n)
    aa, bb = ref_terms(V0, amp)
    sh = np.arange(1, G['nsh'])                      # rings 2..kb
    for it in range(3):
        beta = np.maximum((amp[:, None] * aa).sum(0) / np.maximum((amp[:, None] ** 2 * bb).sum(0), 1e-30), 1e-12)
        res = G['f2'] - 2. * amp[:, None] * beta[None] * aa + (amp[:, None] * beta[None]) ** 2 * bb
        sig2 = np.maximum(res.sum(0) / np.maximum(n * cw, 1.), 1e-30)
        amp = ((beta[None, sh] * aa[:, sh] / sig2[None, sh]).sum(1) /
               np.maximum((beta[None, sh] ** 2 * bb[:, sh] / sig2[None, sh]).sum(1), 1e-30))
        amp = np.clip(amp, 0.05, 20.)
        amp /= np.median(amp)
    raw = (amp[:, None] * aa).sum(0) / np.maximum((amp[:, None] ** 2 * bb).sum(0), 1e-30)
    log('scale / noise: amplitude 5/50/95 %% %.2f %.2f %.2f; per shell (ring, A, beta raw, noise, beta^2 bb / noise):'
        % tuple(np.percentile(amp, [5, 50, 95])))
    snr = (beta ** 2 * (amp[:, None] ** 2 * bb).mean(0)) / np.maximum(sig2 * cw, 1e-30)
    for s in range(1, G['nsh']):
        log('   %2d %6.1f %10.3g %10.3g %8.3g' % (s + 1, G['box_A'] / (s + 1), raw[s], sig2[s], snr[s]))
    return beta, sig2, amp


def suff_stats(V0, Ws, beta, sig2, amp):
    n, q = G['n'], len(Ws)
    w1 = G['sw'] * (beta / sig2)[G['shell']]
    w2 = G['sw'] * (beta ** 2 / sig2)[G['shell']]
    b = np.zeros((n, q)); Gm = np.zeros((n, q, q))
    t0 = time.time()
    for c0 in range(0, n, CH):
        sl = slice(c0, min(n, c0 + CH))
        a = amp[sl][:, None]
        S = np.stack([slices(V, sl) for V in Ws], 1)
        s0 = slices(V0, sl)
        T = G['C'][sl] * w2[None] * a ** 2
        Rr = G['P'][sl] * (w1[None] * a) - s0 * T
        b[sl] = np.einsum('iqm,im->iq', S.real, Rr.real) + np.einsum('iqm,im->iq', S.imag, Rr.imag)
        Sr, Si = S.real, S.imag
        Gm[sl] = np.matmul(Sr * T[:, None], Sr.transpose(0, 2, 1)) + np.matmul(Si * T[:, None], Si.transpose(0, 2, 1))
    log('sufficient statistics: %d particles x %d modes in %.0f s; median trace G %.3g' % (n, q, time.time() - t0,
        np.median(np.trace(Gm, axis1=1, axis2=2))))
    return b, Gm


# ---------------------------------------------------------------- the latent mixture (extreme deconvolution, tied Omega)

def xd_estep(b, Gm, pi, xi, Om):
    Oi = np.linalg.inv(Om)
    H = Gm + Oi[None]
    Hi = np.linalg.inv(H)
    ldH = np.linalg.slogdet(H)[1]
    ldO = np.linalg.slogdet(Om)[1]
    h = b[:, :, None] + (Oi @ xi.T)[None]
    m = Hi @ h
    quad = (h * m).sum(1)
    c = -0.5 * np.einsum('kq,qp,kp->k', xi, Oi, xi)
    ll = np.log(np.maximum(pi, 1e-300))[None] + c[None] + 0.5 * quad - 0.5 * (ldO + ldH)[:, None]
    mx = ll.max(1, keepdims=True)
    lse = mx[:, 0] + np.log(np.exp(ll - mx).sum(1))
    return np.exp(ll - lse[:, None]), m, Hi, lse


def xd_mstep(r, m, Hi):
    N = len(r)
    Nk = r.sum(0) + 1e-12
    xi = np.einsum('ik,iqk->kq', r, m) / Nk[:, None]
    d = m - xi.T[None]
    Om = Hi.mean(0) + np.einsum('ik,iak,ibk->ab', r, d, d) / N
    return Nk / N, xi, 0.5 * (Om + Om.T)


def xd_fit(b, Gm, r0, Om0, mu0=None, maxit=300, tol=1e-6, fixed=False):
    """EM from initial responsibilities r0 (N, K) under the start Omega Om0 centred on mu0; fixed=True keeps r0 (a labels fit).
    Stops when the log-likelihood changes by less than tol nats per particle."""
    N, q = b.shape
    Oi0 = np.linalg.inv(Om0)
    Hi0 = np.linalg.inv(Gm + Oi0[None])
    mu0 = np.zeros(q) if mu0 is None else mu0
    m0 = Hi0 @ (b + Oi0 @ mu0)[:, :, None]
    pi, xi, Om = xd_mstep(r0, np.repeat(m0, r0.shape[1], 2), Hi0)
    if not fixed:
        Om = Om0.copy()
    prev = -np.inf
    for it in range(maxit):
        r, m, Hi, lse = xd_estep(b, Gm, pi, xi, Om)
        LL = float(lse.sum())
        if fixed:
            r = r0
        pi, xi, Om = xd_mstep(r, m, Hi)
        if it > 5 and abs(LL - prev) < tol * N:
            break
        prev = LL
    r, m, Hi, lse = xd_estep(b, Gm, pi, xi, Om)
    return dict(pi=pi, xi=xi, Om=Om, LL=float(lse.sum()), it=it + 1)


def kmeanspp(x, K, rng, its=25):
    c = [x[rng.integers(len(x))]]
    for _ in range(1, K):
        d = np.min(((x[:, None] - np.array(c)[None]) ** 2).sum(-1), 1)
        c.append(x[rng.choice(len(x), p=d / d.sum())])
    c = np.array(c)
    for _ in range(its):
        lab = np.argmin(((x[:, None] - c[None]) ** 2).sum(-1), 1)
        for k in range(K):
            if (lab == k).any():
                c[k] = x[lab == k].mean(0)
    return lab


def onehot(lab, K):
    r = np.zeros((len(lab), K)); r[np.arange(len(lab)), lab] = 1.
    return r


def fit_one(sel, Om0):
    """K = 1 from Om0 (converged harder: it is every gain's baseline)."""
    b, Gm = G['b'][sel], G['Gm'][sel]
    return xd_fit(b, Gm, np.ones((len(sel), 1)), Om0, maxit=2000)


def fit_job(job):
    """one EM start: (tag, sel, K, seed, Om1, evalsel) -> the fit (best kept by the caller) and the held-out per-particle lse."""
    tag, sel, K, seed, Om1, evalsel = job
    b, Gm, mu1 = G['b'][sel], G['Gm'][sel], G['mu1']
    if K == 1:
        f = xd_fit(b, Gm, np.ones((len(sel), 1)), Om1, mu1, maxit=2000)
    else:
        Oi = np.linalg.inv(Om1)
        Hi = np.linalg.inv(Gm + Oi[None])
        m0 = (Hi @ (b + Oi @ mu1)[:, :, None])[:, :, 0]
        lam, U = np.linalg.eigh(Om1)
        x = (m0 - mu1) @ U / np.sqrt(np.maximum(lam, 1e-30))
        lab = kmeanspp(x, K, np.random.default_rng(seed))
        f = xd_fit(b, Gm, onehot(lab, K) * 0.9 + 0.1 / K, Om1, mu1)
    out = dict(tag=tag, K=K, seed=seed, LL=f['LL'], it=f['it'], pi=f['pi'], xi=f['xi'], Om=f['Om'])
    if evalsel is not None:
        out['held'] = xd_estep(G['b'][evalsel], G['Gm'][evalsel], f['pi'], f['xi'], f['Om'])[3]
    return out


def run_jobs(jobs, nproc):
    if nproc <= 1:
        return [fit_job(j) for j in jobs]
    with Pool(nproc) as pool:
        return pool.map(fit_job, jobs, chunksize=1)


# ---------------------------------------------------------------- readouts

def entropy(lab):
    _, c = np.unique(lab, return_counts=True); p = c / c.sum()
    return float(-(p * np.log(p)).sum())


def nmi(a, b):
    ua, ia = np.unique(a, return_inverse=True); ub, ib = np.unique(b, return_inverse=True)
    t = np.zeros((len(ua), len(ub))); np.add.at(t, (ia, ib), 1.); t /= t.sum()
    pa, pb = t.sum(1), t.sum(0); nz = t > 0
    mi = float((t[nz] * np.log(t[nz] / (pa[:, None] * pb[None])[nz])).sum())
    ha, hb = entropy(a), entropy(b)
    return 0. if ha + hb == 0 else 2. * mi / (ha + hb)


def ari(a, b):
    ua, ia = np.unique(a, return_inverse=True); ub, ib = np.unique(b, return_inverse=True)
    t = np.zeros((len(ua), len(ub))); np.add.at(t, (ia, ib), 1.)
    c2 = lambda x: x * (x - 1) / 2.
    s, sa, sb, n2 = c2(t).sum(), c2(t.sum(1)).sum(), c2(t.sum(0)).sum(), c2(t.sum())
    e = sa * sb / max(n2, 1.)
    return float((s - e) / max(0.5 * (sa + sb) - e, 1e-12))


def eta2(lab, x):
    ok = np.isfinite(x)
    lab, x = lab[ok], x[ok]
    mu = x.mean(); tot = ((x - mu) ** 2).sum()
    btw = sum((lab == u).sum() * (x[lab == u].mean() - mu) ** 2 for u in np.unique(lab))
    return float(btw / max(tot, 1e-30))


def load_truth(spec, node, prow):
    kind, _, rest = spec.partition(':')
    if kind == 'none':
        return kind, {}
    if kind == 'qtoy':
        hdr = open(rest).readline().lstrip('#').split()
        t = np.loadtxt(rest, ndmin=2)
        assert (t[:, 0] == np.arange(1, len(t) + 1)).all(), 'truth rows are not 1..N'
        T = t[prow - 1]
        return kind, {c: T[:, hdr.index(c)] for c in ('combo', 'conformer', 'junk', 'snrq', 'bfac')}
    if kind in ('igmfix', 'igmnull'):
        t = np.loadtxt(rest, ndmin=2, comments='#')
        assert (t[:, 0] == np.arange(1, len(t) + 1)).all(), 'truth rows are not 1..N'
        T = t[prow - 1]
        if kind == 'igmfix':
            return kind, dict(combo=T[:, 1], **{'fa%d' % j: T[:, 1 + j] for j in range(1, 6)})
        return kind, dict(**{'na%d' % j: T[:, 6 + j] for j in range(1, 6)}, **{'ns%d' % j: T[:, 11 + j] for j in range(1, 6)})
    if kind == 'closed':
        fxr, clr = rest.split(':')
        ref_rows = np.loadtxt(clr, dtype=int, ndmin=1)
        root = np.loadtxt(os.path.join(os.path.dirname(node.rstrip('/')), 'node_r', 'rows.txt'), dtype=int, ndmin=1)
        rm = np.loadtxt(fxr, dtype=int, ndmin=1) if os.path.exists(fxr) else None
        ix = rm[prow - 1] if (rm is not None and root.max() <= len(rm)) else prow
        return kind, dict(closedref=np.isin(ix, ref_rows).astype(float))
    raise SystemExit('truth kind ' + kind)


def readout(kind, T, lab, post=None):
    """hard labels (and optional posteriors (N, K) in label order 0..K-1) against the truth."""
    out = dict(K=int(len(np.unique(lab))), sizes=[int((lab == u).sum()) for u in np.unique(lab)])
    if kind == 'qtoy':
        out.update(nmi_combo=nmi(lab, T['combo']), ari_combo=ari(lab, T['combo']),
                   eta2_snrq=eta2(lab, T['snrq']), eta2_bfac=eta2(lab, T['bfac']))
        tgt = T['combo']
    elif kind == 'igmfix':
        out.update({'nmi_' + c: nmi(lab, T[c]) for c in ('combo', 'fa1', 'fa2', 'fa3', 'fa4', 'fa5')})
        out['ari_fa1'] = ari(lab, T['fa1'])
        tgt = T['fa1']
    elif kind == 'igmnull':
        out.update({'eta2_' + c: eta2(lab, T[c]) for c in ('na1', 'na2', 'na3', 'na4', 'na5')})
        tgt = None
    elif kind == 'closed':
        r = T['closedref'] > 0.5
        nref, share = int(r.sum()), float(r.mean())
        gr = [(int((lab == u).sum()), int((r & (lab == u)).sum())) for u in np.unique(lab)]
        out.update(nref=nref, ref_share=share, groups=gr,
                   separates=bool(nref >= 500 and any(g[1] >= 0.5 * nref and g[1] / max(g[0], 1) >= 2 * share for g in gr)),
                   best_enrich=float(max(g[1] / max(g[0], 1) for g in gr) / max(share, 1e-12)))
        tgt = None
    else:
        tgt = None
    if post is not None and tgt is not None:
        cal = []
        for u in np.unique(lab):
            mm = lab == u
            vals, cnt = np.unique(tgt[mm], return_counts=True)
            cal.append(dict(size=int(mm.sum()), purity=float(cnt.max() / mm.sum()), mean_post=float(post[mm, u].mean())))
        out['calibration'] = cal
        out['mean_max_post'] = float(post.max(1).mean())
        out['mean_purity'] = float(sum(c['size'] * c['purity'] for c in cal) / len(lab))
    return out


def dump_states(fn, didx, n0):
    from ppinfo_site import read_dump
    if not os.path.exists(fn):
        return None
    d = read_dump(fn)
    if len(d['state']) != n0:
        log('WARNING %s has %d rows for %d' % (fn, len(d['state']), n0))
        return None
    return d['state'][didx].astype(int)


def relabel(lab):
    u, inv = np.unique(lab, return_inverse=True)
    return inv


# ---------------------------------------------------------------- selftest

def selftest(nproc):
    rng = np.random.default_rng(11)
    N, q, K = 9000, 4, 3
    xi = np.array([[0, 0, 0, 0], [2., 0, 0, 0], [0, 2.5, 0, 0]])
    lab = rng.choice(K, N, p=[.5, .3, .2])
    z = xi[lab] + rng.normal(0, .4, (N, q))
    A = rng.normal(0, 1, (N, q, q)) * np.exp(rng.normal(0, .7, N))[:, None, None]
    Gm = A @ A.transpose(0, 2, 1) * 0.3 + 0.05 * np.eye(q)[None]
    zhat = z + np.array([rng.multivariate_normal(np.zeros(q), np.linalg.inv(g)) for g in Gm])
    b = (Gm @ zhat[:, :, None])[:, :, 0]
    G.update(b=b, Gm=Gm, mu1=np.zeros(q))
    f1 = fit_one(np.arange(N), np.eye(q))
    G['mu1'] = f1['xi'][0]
    ia, ib = np.arange(0, N, 2), np.arange(1, N, 2)
    jobs = [('h', ia, k, s, f1['Om'], ib) for k in (1, 2, 3, 4, 5) for s in range(3)]
    res = run_jobs(jobs, nproc)
    gain = {}
    for k in (1, 2, 3, 4, 5):
        best = max((r for r in res if r['K'] == k), key=lambda r: r['LL'])
        gain[k] = float(best['held'].sum())
    base = gain[1]
    log('selftest held-out gain over K=1 (nats): ' + ' '.join('K%d %.1f' % (k, gain[k] - base) for k in gain))
    jobs = [('a', np.arange(N), K, s, f1['Om'], None) for s in range(3)]
    best = max(run_jobs(jobs, nproc), key=lambda r: r['LL'])
    r = xd_estep(b, Gm, best['pi'], best['xi'], best['Om'])[0]
    est = r.argmax(1)
    log('selftest K=3: NMI %.3f (truth), pi %s, mean max posterior %.3f, accuracy-proxy purity %.3f'
        % (nmi(est, lab), np.round(best['pi'], 3), r.max(1).mean(), readout('qtoy', dict(combo=lab, snrq=np.zeros(N), bfac=np.zeros(N)), est, r)['mean_purity']))


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--node'); ap.add_argument('--dump'); ap.add_argument('--out'); ap.add_argument('--cache', default='')
    ap.add_argument('--truth', default='none'); ap.add_argument('--lp', type=float, default=15.)
    ap.add_argument('--mskdiam', type=float, default=0.); ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--kmax', type=int, default=8); ap.add_argument('--restarts', type=int, default=5)
    ap.add_argument('--maxn', type=int, default=0); ap.add_argument('--selftest', action='store_true')
    a = ap.parse_args()
    if a.selftest:
        selftest(a.nproc)
        return
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    setup(a)
    n, n0, didx = G['n'], G['n0'], G['didx']
    from mrcio import read_mrc
    labs = {}
    for tag, fn in (('P', 'b_ptcl3D.txt'), ('R0', 'r0_ptcl3D.txt'), ('F', 'f_ptcl3D.txt')):
        s = dump_states(os.path.join(a.dump, fn), didx, n0)
        if s is not None:
            labs[tag] = s
            log('%s labels: %s' % (tag, ' '.join('%d:%d' % (u, (s == u).sum()) for u in np.unique(s))))
    kind, T = load_truth(a.truth, a.node, G['prow'])
    KP = int(len(np.unique(labs['P'][labs['P'] > 0])))
    # basis maps
    V0 = vol_ft(read_mrc(G['vol1'])[0])
    pcs = sorted(f for f in glob.glob(os.path.join(a.node, 'fflex', 'flex_pca_pc[0-9][0-9][0-9].mrc')))
    Vf = [vol_ft(read_mrc(f)[0]) for f in pcs]
    beta, sig2, amp = scale_noise(V0)
    sf_maps, sf_info = site_factors(a.node, labs.get('R0', np.zeros(n, int)))
    Vs = [vol_ft(m) for m in sf_maps]
    bases = {'L': (Vf, ['pc%d' % (j + 1) for j in range(len(Vf))])}
    if Vs:
        bases['LS'] = (Vf + Vs, ['pc%d' % (j + 1) for j in range(len(Vf))] + ['sf%d' % (j + 1) for j in range(len(Vs))])
    eo = G['eo']
    ia, ib = np.where(eo == 0)[0], np.where(eo == 1)[0]
    summary = dict(node=a.node, n=n, KP=KP, kind=kind, site_factors=sf_info, bases={}, arms={}, heldout={}, agreement={})
    allpost = {}
    for bname, (Vlist, names) in bases.items():
        W, resid = orthonormal(Vlist, names)
        b, Gm = suff_stats(V0, W, beta, sig2, amp)
        G.update(b=b, Gm=Gm)
        q = b.shape[1]
        Om0 = np.eye(q) / max(np.median(np.trace(Gm, axis1=1, axis2=2)) / q, 1e-30) * 10.
        G['mu1'] = np.zeros(q)
        f1 = fit_one(np.arange(n), Om0)
        Om1 = f1['Om']
        G['mu1'] = f1['xi'][0]
        log('%s: q %d, K=1 fit LL/particle %.4f (%d its), Omega eigenvalues %s' % (bname, q, f1['LL'] / n, f1['it'],
            np.array2string(np.sort(np.linalg.eigvalsh(Om1))[::-1], precision=3)))
        # cross-half held-out curve over K, both directions
        jobs = []
        for k in range(1, a.kmax + 1):
            for s in range(1 if k == 1 else 3):
                jobs += [('e', ia, k, 100 * k + s, Om1, ib), ('o', ib, k, 100 * k + s, Om1, ia)]
        res = run_jobs(jobs, a.nproc)
        held = {}
        for k in range(1, a.kmax + 1):
            tot = np.zeros(0)
            for d, ev in (('e', ib), ('o', ia)):
                best = max((r for r in res if r['K'] == k and r['tag'] == d), key=lambda r: r['LL'])
                tot = np.r_[tot, best['held']]
            held[k] = tot
        gains = {k: float((held[k] - held[1]).sum()) for k in held}
        ses = {k: float((held[k] - held[1]).std() * np.sqrt(len(held[k]))) for k in held}
        kstar = max(gains, key=lambda k: gains[k] - (2. * ses[k] if k > 1 else 0.))
        kstar = 1 if gains[kstar] - 2. * ses[kstar] <= 0 else kstar
        log('%s held-out gain over K=1 (nats, both halves; +-SE): %s -> own K %d'
            % (bname, ' '.join('K%d %.1f+-%.1f' % (k, gains[k], ses[k]) for k in gains), kstar))
        summary['heldout'][bname] = dict(gain=gains, se=ses, kstar=kstar, q=q, resid=resid, names=names)
        # full fits at the production K and at the own K
        for kk, tag in ((KP, '%s@KP' % bname), (kstar, '%s@own' % bname)):
            if kk < 1 or (tag.endswith('own') and kk == KP):
                continue
            jobs = [('a', np.arange(n), kk, 1000 + s, Om1, None) for s in range(a.restarts if kk > 1 else 1)]
            best = max(run_jobs(jobs, a.nproc), key=lambda r: r['LL'])
            post = xd_estep(b, Gm, best['pi'], best['xi'], best['Om'])[0]
            hard = post.argmax(1)
            allpost[tag] = post
            summary['arms'][tag] = readout(kind, T, hard, post)
            summary['arms'][tag].update(LL_per_ptcl=best['LL'] / n, pi=np.round(best['pi'], 4).tolist(),
                                        soft_occupancy=np.round(post.mean(0), 4).tolist())
            labs[tag] = hard + 1
            log('%s: K %d, LL/particle %.4f, pi %s, mean max posterior %.3f' % (tag, kk, best['LL'] / n, np.round(best['pi'], 3),
                post.max(1).mean()))
        # every arm's labels held out in this latent space: class Gaussians fitted on one half with the labels, scored on the other
        for tag in ('P', 'R0', 'F'):
            if tag not in labs:
                continue
            lab = relabel(labs[tag])
            K = lab.max() + 1
            tot = 0.
            per = []
            for fs, es in ((ia, ib), (ib, ia)):
                f = xd_fit(b[fs], Gm[fs], onehot(lab[fs], K), Om1, G['mu1'], fixed=True)
                lse = xd_estep(b[es], Gm[es], f['pi'], f['xi'], f['Om'])[3]
                per.append(lse)
            per = np.r_[per[0], per[1]]
            # held[1] concatenates (even->odd) then (odd->even) held-out K=1 scores: odd particles first, then even
            g = float((per - held[1]).sum())
            summary['heldout'][bname]['labels_' + tag] = dict(K=int(K), gain=g, se=float((per - held[1]).std() * np.sqrt(len(per))))
            log('%s: %s labels (K %d) held out in this latent space: gain over K=1 %.1f nats' % (bname, tag, K, g))
    # arm readouts and agreement
    for tag in ('P', 'R0', 'F'):
        if tag in labs:
            summary['arms'][tag] = readout(kind, T, relabel(labs[tag]))
    for tag in summary['arms']:
        if tag != 'P' and tag in labs:
            summary['agreement'][tag] = dict(nmi_P=nmi(labs[tag], labs['P']), ari_P=ari(labs[tag], labs['P']))
    summary['wall'] = time.time() - t0
    json.dump(summary, open(os.path.join(a.out, 'summary.json'), 'w'), indent=1, default=float)
    cols = sorted(labs)
    with open(os.path.join(a.out, 'labels.txt'), 'w') as fo:
        fo.write('# row ' + ' '.join(cols) + ' ' + ' '.join('post_%s_%d' % (t, j + 1) for t in allpost for j in range(allpost[t].shape[1])) + '\n')
        for i in range(n):
            fo.write('%d %s %s\n' % (G['prow'][i], ' '.join(str(int(labs[c][i])) for c in cols),
                     ' '.join('%.4f' % allpost[t][i, j] for t in allpost for j in range(allpost[t].shape[1]))))
    log('==== %s (%s, N %d, production K %d)' % (a.node, kind, n, KP))
    for tag, r in summary['arms'].items():
        keys = [k for k in r if k not in ('calibration', 'groups', 'sizes', 'pi', 'soft_occupancy')]
        log('  %-8s %s | agree P %s' % (tag, ' '.join('%s=%s' % (k, ('%.3f' % r[k]) if isinstance(r[k], float) else r[k]) for k in keys),
            ('%.3f' % summary['agreement'][tag]['nmi_P']) if tag in summary['agreement'] else '-'))
        if 'calibration' in r:
            log('           calibration (size purity mean_post): ' + '; '.join('%d %.2f %.2f' % (c['size'], c['purity'], c['mean_post'])
                for c in r['calibration']))
        if 'groups' in r:
            log('           groups (n, closed refs): ' + ' '.join('%d/%d' % g for g in r['groups']))
    log('done in %.0f s' % summary['wall'])


if __name__ == '__main__':
    main()
