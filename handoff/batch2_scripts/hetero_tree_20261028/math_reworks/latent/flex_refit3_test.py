#!/usr/bin/env python3
"""flex_refit_test.py: a Python prototype of two flex_pca changes, refitting the node's PPCA basis at the refinement's (stage-B)
poses, frozen, with per-particle reliability weights from the refinement's objective, judged by what the basis does downstream.
Bases: orig = flex_pca's delivered modes (fitted at the parent's poses); unw_in = the refit below at the input poses, unweighted
(the estimator control: unw_in vs unw is the poses alone); unw = refit at stage-B poses, unweighted; wt = refit at
stage-B poses, each particle weighted by its probability of fitting its stage-B state map as its noise predicts (a two-component
mixture on the log residual ratio, outliers = misaligned or junk). The refit is flex's EM in the Fourier lattice: E-step through
each particle's own precision (prior N(xi, Omega) from a one-component fit), M-step a q x q solve per voxel, band to 15 A,
orthonormal modes; it starts from flex's modes. Each basis, alone and with the node's site factors, is judged by the truth-trained
ceiling in its latent and by the calibrated stage-B-state posterior (calib_flexprior3's, two iterations).
usage: flex_refit_test.py --node NODE --dump DUMP --out OUT --cache planes.npy --truth KIND:ARGS [--maxn 15000 --refit_its 6]"""
import argparse, glob, json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import latent_states_v4 as LS
import tree_estep_v3 as TE
import flexprior_estep as FP
import calib_flexprior as CF

G = LS.G
log = LS.log
XD_SERIAL, SS_SERIAL = LS.xd_estep, CF.suff_stats_masked


def threaded(nthr):
    """route the latent mixture's E-step and the latent sufficient statistics (both independent per particle; numpy's batched
    linear algebra and einsum release the GIL) through nthr threads; the results are the serial ones"""
    from concurrent.futures import ThreadPoolExecutor
    pool = ThreadPoolExecutor(nthr)

    def xd_estep(b, Gm, pi, xi, Om):
        if len(b) < 64 * nthr:
            return XD_SERIAL(b, Gm, pi, xi, Om)
        e = np.linspace(0, len(b), nthr + 1).astype(int)
        parts = list(pool.map(lambda s: XD_SERIAL(b[s[0]:s[1]], Gm[s[0]:s[1]], pi, xi, Om), zip(e[:-1], e[1:])))
        return tuple(np.concatenate([p[j] for p in parts]) for j in range(4))

    def suff_stats_masked(V0, Ws, beta, sig2, amp, smask):
        n, q = G['n'], len(Ws)
        w1 = np.where(smask, G['sw'] * (beta / sig2)[G['shell']], 0.)
        w2 = np.where(smask, G['sw'] * (beta ** 2 / sig2)[G['shell']], 0.)
        b = np.zeros((n, q)); Gm = np.zeros((n, q, q))

        def chunk(c0):
            sl = slice(c0, min(n, c0 + 400))
            a = amp[sl][:, None]
            S = np.stack([LS.slices(V, sl) for V in Ws], 1)
            s0 = LS.slices(V0, sl)
            T = G['C'][sl] * w2[None] * a ** 2
            Rr = G['P'][sl] * (w1[None] * a) - s0 * T
            b[sl] = np.einsum('iqm,im->iq', S.real, Rr.real) + np.einsum('iqm,im->iq', S.imag, Rr.imag)
            Gm[sl] = np.matmul(S.real * T[:, None], S.real.transpose(0, 2, 1)) + np.matmul(S.imag * T[:, None], S.imag.transpose(0, 2, 1))
        list(pool.map(chunk, range(0, n, 400)))
        return b, Gm

    LS.xd_estep, CF.suff_stats_masked = xd_estep, suff_stats_masked
    log('threads: latent E-step and sufficient statistics on %d threads' % nthr)


def ppca_refit(W0, V0, beta, sig2, amp, wts, low, kf, nits):
    """flex's PPCA EM on the node's particles at the current geometry, weighted by wts, from the modes W0."""
    W = list(W0)
    q = len(W)
    NC, shell, sw = G['NC'], G['shell'], G['sw']
    c1 = np.where(low, sw * (beta / sig2)[shell], 0.)
    c2 = np.where(low, sw * (beta ** 2 / sig2)[shell], 0.)
    cs = np.r_[G['cell_shell'], 10 ** 6]
    for it in range(nits):
        b, Gm = CF.suff_stats_masked(V0, W, beta, sig2, amp, low)
        G.update(b=b, Gm=Gm, mu1=np.zeros(q))
        f1 = LS.fit_one(np.arange(G['n']), np.eye(q) / max(np.median(np.trace(Gm, axis1=1, axis2=2)) / q, 1e-30) * 10.)
        Oi = np.linalg.inv(f1['Om'])
        H = Gm + Oi[None]
        Hi = np.linalg.inv(H)
        m = (Hi @ (b + Oi @ f1['xi'][0])[:, :, None])[:, :, 0]
        M = Hi + m[:, :, None] * m[:, None, :]
        A = np.zeros((NC + 1, q, q)); Bv = np.zeros((NC + 1, q), complex)
        for c0 in range(0, G['n'], 500):
            sl = slice(c0, min(G['n'], c0 + 500))
            ii, cj, mm = G['idx'][sl], G['cnj'][sl], G['mate'][sl]
            a = amp[sl][:, None]; ww = wts[sl][:, None]
            s0 = LS.slices(V0, sl)
            R = G['P'][sl] * (c1[None] * a) - s0 * (G['C'][sl] * c2[None] * a ** 2)
            T = G['C'][sl] * c2[None] * a ** 2 * ww
            half = mm >= 0
            wf = np.where(half, 0.5, 1.)
            f, fm = ii.ravel(), mm[half]
            for j in range(q):
                v = R * (ww * m[sl, j][:, None])
                v = np.where(cj, np.conj(v), v)
                Bv[:, j] += np.bincount(f, (v.real * wf).ravel(), NC + 1) + 1j * np.bincount(f, (v.imag * wf).ravel(), NC + 1)
                if half.any():
                    Bv[:, j] += np.bincount(fm, 0.5 * v[half].real, NC + 1) - 1j * np.bincount(fm, 0.5 * v[half].imag, NC + 1)
                for k in range(j, q):
                    t = T * M[sl, j, k][:, None]
                    A[:, j, k] += np.bincount(f, (t * wf).ravel(), NC + 1)
                    if half.any():
                        A[:, j, k] += np.bincount(fm, 0.5 * t[half], NC + 1)
        iu = np.triu_indices(q, 1)
        A[:, iu[1], iu[0]] = A[:, iu[0], iu[1]]
        # a ridge per ring at a thousandth of the ring's median diagonal: unsampled cells stay zero; rings without data (beyond
        # flex's band, zeroed below) get a unit ridge
        dg = np.einsum('vjj->v', A) / q
        ridge = np.zeros(NC + 1)
        for s in range(G['nyq'] + 1):
            sel = (cs == s) & (dg > 0)
            if sel.any():
                ridge[cs == s] = 1e-3 * np.median(dg[sel])
        ridge[ridge <= 0.] = 1.
        A += ridge[:, None, None] * np.eye(q)[None]
        X = np.linalg.solve(A, Bv[:, :, None])[:, :, 0]
        X[cs > kf] = 0.
        W, _ = LS.orthonormal([X[:, j].copy() for j in range(q)], ['m%d' % (j + 1) for j in range(q)])
        log('   refit it %d: latent log-lik per particle %.4f, Omega trace %.3g, modes %d' % (it + 1, f1['LL'] / G['n'], np.trace(f1['Om']), len(W)))
    return W


def reliability(refs, beta, sig2, lab, ia, ib):
    """per particle the residual ratio against its own state's cross-half reference over the band, and its probability of
    being in-model under a two-component mixture on the log ratio (outliers above)."""
    shell, sw, nsh = G['shell'], G['sw'], G['nsh']
    cw = np.bincount(shell, sw, nsh)
    rr = np.zeros(G['n'])
    for h, other in ((ia, refs[1]), (ib, refs[0])):
        for c0 in range(0, len(h), 500):
            p = h[c0:c0 + 500]
            a = G['cn'][p][:, None]
            mref = np.zeros((len(p), len(shell)), complex)
            for c in np.unique(lab[p]):
                wm = lab[p] == c
                mref[wm] = LS.slices(other[c], p[wm])
            aa = (mref.real * G['P'][p].real + mref.imag * G['P'][p].imag) * sw[None]
            bb = np.abs(mref) ** 2 * G['C'][p] * sw[None]
            key = (np.arange(len(p))[:, None] * nsh + shell[None, :]).ravel()
            AA = np.bincount(key, aa.ravel(), len(p) * nsh).reshape(len(p), nsh)
            BB = np.bincount(key, bb.ravel(), len(p) * nsh).reshape(len(p), nsh)
            res = G['f2'][p] - 2. * a * beta[None] * AA + (a * beta[None]) ** 2 * BB
            rr[p] = (res / sig2[None]).sum(1) / cw.sum()
    x = np.log(np.maximum(rr, 1e-12))
    med = np.median(x); mad = 1.4826 * np.median(np.abs(x - med))
    mu = np.array([med, med + 4. * mad]); sd = np.array([mad, 3. * mad]); pi = np.array([0.95, 0.05])
    for _ in range(200):
        ll = np.log(pi)[None] - 0.5 * ((x[:, None] - mu[None]) / sd[None]) ** 2 - np.log(sd)[None]
        r = np.exp(ll - ll.max(1, keepdims=True)); r /= r.sum(1, keepdims=True)
        pi = r.mean(0); mu = (r * x[:, None]).sum(0) / r.sum(0)
        sd = np.sqrt((r * (x[:, None] - mu[None]) ** 2).sum(0) / r.sum(0)) + 1e-6
        if mu[1] < mu[0]:
            break
    w = r[:, 0] if mu[1] >= mu[0] else np.ones(G['n'])
    log('reliability: residual ratio median %.3f, mixture in-model %.3f at %.3f +- %.3f, outliers at %.3f; weights 5/50 %% %.2f %.2f'
        % (float(np.median(rr)), pi[0], np.exp(mu[0]), sd[0], np.exp(mu[1]), *np.percentile(w, [5, 50])))
    return w, rr


def calibrate(W, labf, K, kind, T, eps, ia, ib, vA, vB, nits, V0, beta_mu, sig2_mu, amp, kf, tag):
    """calib_flexprior3's calibrated stage-B-state posterior with the latent of W (no refit), two iterations; and the ceiling."""
    allm = np.ones(len(vA), bool)
    b, Gm = CF.suff_stats_masked(V0, W, beta_mu, sig2_mu, amp, allm)
    bA, GA = CF.suff_stats_masked(V0, W, beta_mu, sig2_mu, amp, vA)
    bB, GB = CF.suff_stats_masked(V0, W, beta_mu, sig2_mu, amp, vB)
    n, q = b.shape
    G.update(b=b, Gm=Gm, mu1=np.zeros(q))
    f1 = LS.fit_one(np.arange(n), np.eye(q) / max(np.median(np.trace(Gm, axis1=1, axis2=2)) / q, 1e-30) * 10.)
    Om1, G['mu1'] = f1['Om'], f1['xi'][0]
    out = dict(tag=tag, q=q)
    cl = LS.ceiling(kind, T, b, Gm, ia, ib, Om1, G['mu1'])
    if cl is not None:
        out['ceiling'] = cl
    w, wA, wB = [np.zeros((n, K)) for _ in range(3)]
    for fs, es in ((ia, ib), (ib, ia)):
        f = LS.xd_fit(b[fs], Gm[fs], LS.onehot(labf[fs], K), Om1, G['mu1'], fixed=True)
        for dst, bb, gg in ((w, b, Gm), (wA, bA, GA), (wB, bB, GB)):
            dst[es] = LS.xd_estep(bb[es], gg[es], f['pi'], f['xi'], f['Om'])[0]
    refs, _, _ = TE.references(LS.onehot(labf, K), ia, ib, eps)
    beta, sig2 = TE.ref_scale_noise([refs[1], refs[0]], [ia, ib], labf)
    tgt = {'qtoy': 'combo', 'closed': 'closedref'}.get(kind)
    lab = labf
    out['its'] = []
    for it in range(1, nits + 1):
        if it > 1:
            refs, _, _ = TE.references(LS.onehot(lab, K), ia, ib, eps)
            beta, sig2 = TE.ref_scale_noise([refs[1], refs[0]], [ia, ib], lab)
        lA, lB = np.zeros((n, K)), np.zeros((n, K))
        for hh_, other in ((ia, refs[1]), (ib, refs[0])):
            lA[hh_] = FP.score_band(hh_, other, beta, sig2, vA)
            lB[hh_] = FP.score_band(hh_, other, beta, sig2, vB)
        S = {(bt, al): CF.predictive(lA, lB, wA, wB, bt, al) for bt in CF.BETAS for al in CF.ALPHAS}
        bstar, astar = max(S, key=S.get)
        row = dict(it=it, beta=bstar, alpha=astar, gain=(S[(bstar, astar)] - S[(1.0, 0.)]) / n)
        for arm, (bt, al) in (('cal', (bstar, astar)), ('obj', (1.0, 0.))):
            r = FP.softmax(bt * (lA + lB) + al * np.log(np.maximum(w, 1e-300)))
            hl = r.argmax(1)
            if arm == 'cal':
                lab = hl
            if tgt is not None and len(np.unique(T[tgt])) > 1:
                rd = LS.readout(kind, T, hl, r)
                for kk in ('nmi_combo', 'best_enrich'):
                    if kk in rd:
                        row['%s_%s' % (arm, kk)] = rd[kk]
            row['%s_mmp' % arm] = float(r.max(1).mean())
        out['its'].append(row)
        log('   [%s] it %d: beta %.2f alpha %.2f gain %+.4f%s' % (tag, it, bstar, astar, row['gain'],
            ''.join(' %s %.3f' % (kk, row[kk]) for kk in sorted(row) if kk.startswith(('cal_', 'obj_')))))
    if cl is not None:
        log('   [%s] ceiling %s' % (tag, ' '.join('%s=%s' % (k, ('%.3f' % v) if isinstance(v, float) else v) for k, v in cl.items())))
    return out


def load_node(a, poses):
    """setup at the given poses; the node mean with its scale and noise, flex's modes in flex's band"""
    from mrcio import read_mrc
    a.poses = poses
    LS.setup(a)
    kf = int(round(G['box_A'] / FP.LPF))
    V0 = LS.vol_ft(read_mrc(G['vol1'])[0])
    beta_mu, sig2_mu, amp = LS.scale_noise(V0)
    G['cn'] = amp
    pcs = sorted(glob.glob(os.path.join(a.pcdir or os.path.join(a.node, 'fflex'), 'flex_pca_pc[0-9][0-9][0-9].mrc')))
    log('flex modes: %d from %s' % (len(pcs), a.pcdir or os.path.join(a.node, 'fflex')))
    Worig, _ = LS.orthonormal([FP.lowpass_cells(LS.vol_ft(read_mrc(f)[0]), kf) for f in pcs], ['pc%d' % (j + 1) for j in range(len(pcs))])
    return V0, beta_mu, sig2_mu, amp, Worig, kf, (G['shell'] + 1) <= kf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--node'); ap.add_argument('--dump'); ap.add_argument('--out'); ap.add_argument('--cache', default='')
    ap.add_argument('--truth', default='none'); ap.add_argument('--mskdiam', type=float, default=0.)
    ap.add_argument('--maxn', type=int, default=15000); ap.add_argument('--refit_its', type=int, default=6)
    ap.add_argument('--nit', type=int, default=2); ap.add_argument('--nproc', type=int, default=4)
    ap.add_argument('--arm', default='all', help='orig | unw_in | unw | wt | all: the basis (one per task in the wide array)')
    ap.add_argument('--pcdir', default='', help='flex modes from this folder instead of the node\'s fflex (flex_pca rerun at stage-B poses)')
    ap.add_argument('--sf', default='both', help='no | yes | both: the basis without / with the node\'s site factors')
    a = ap.parse_args()
    a.lp = FP.LPB
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    if a.nproc > 1:
        threaded(a.nproc)
    from ppinfo_site import soft_sphere
    Win = None
    if a.arm in ('unw_in', 'all'):
        V0, beta_mu, sig2_mu, amp, Worig, kf, low = load_node(a, 'input')
        log('refit unweighted at the input poses'); Win = ppca_refit(Worig, V0, beta_mu, sig2_mu, amp, np.ones(G['n']), low, kf, a.refit_its)
    V0, beta_mu, sig2_mu, amp, Worig, kf, low = load_node(a, 'stageB')
    n, n0, didx = G['n'], G['n0'], G['didx']
    G['sm_c'] = np.fft.ifftshift(soft_sphere(G['B'], G['msk_crop']))
    B, kb = G['B'], G['kb']
    kk, hh = np.meshgrid(np.r_[0:B // 2, -B // 2:0], np.arange(B // 2), indexing='ij')
    r2 = hh * hh + kk * kk
    sel = (r2 >= 3) & (r2 <= kb * (kb + 1))
    h, k = hh[sel], kk[sel]
    vA = ((h + k) % 2) == 0; vB = ~vA
    labf = LS.relabel(LS.dump_states(os.path.join(a.dump, 'b_ptcl3D.txt'), didx, n0))
    K = labf.max() + 1
    kind, T = LS.load_truth(a.truth, a.node, G['prow'])
    R0 = LS.dump_states(os.path.join(a.dump, 'r0_ptcl3D.txt'), didx, n0)
    sfm = LS.site_factors(a.node, R0 if R0 is not None else np.zeros(n, int))[0] if a.sf != 'no' else []
    Wsf = [FP.lowpass_cells(LS.vol_ft(m), kf) for m in sfm]
    eo = G['eo']
    ia, ib = np.where(eo == 0)[0], np.where(eo == 1)[0]
    eps = 0.5 * (TE.density(ia) + TE.density(ib))
    refs, _, _ = TE.references(LS.onehot(labf, K), ia, ib, eps)
    beta, sig2 = TE.ref_scale_noise([refs[1], refs[0]], [ia, ib], labf)
    wrel, rr = reliability(refs, beta, sig2, labf, ia, ib)
    summary = dict(node=a.node, n=n, K=int(K), kind=kind, rel_w_median=float(np.median(wrel)),
                   rel_w_frac_below_half=float((wrel < 0.5).mean()))
    if kind == 'qtoy':
        summary['rel_corr_snrq'] = float(np.corrcoef(wrel, T['snrq'])[0, 1])
        if T['junk'].any():
            summary['rel_junk_mean_w'] = float(wrel[T['junk'] > 0].mean())
        log('reliability vs truth: corr with snrq %.3f, junk mean weight %s' % (summary['rel_corr_snrq'], summary.get('rel_junk_mean_w', '-')))
    if kind == 'closed':
        cr = T['closedref'] > 0.5
        summary['rel_closed_mean_w'], summary['rel_other_mean_w'] = float(wrel[cr].mean()), float(wrel[~cr].mean())
        log('reliability: closed reference rows mean weight %.3f, other rows %.3f' % (summary['rel_closed_mean_w'], summary['rel_other_mean_w']))
    Wunw = Wwt = None
    if a.arm in ('unw', 'all'):
        log('refit unweighted at stage-B poses'); Wunw = ppca_refit(Worig, V0, beta_mu, sig2_mu, amp, np.ones(n), low, kf, a.refit_its)
    if a.arm in ('wt', 'all'):
        log('refit weighted at stage-B poses'); Wwt = ppca_refit(Worig, V0, beta_mu, sig2_mu, amp, wrel, low, kf, a.refit_its)
    summary['variants'] = []
    for tag, Wl in (('orig', Worig), ('unw_in', Win), ('unw', Wunw), ('wt', Wwt)):
        if Wl is None or a.arm not in (tag, 'all'):
            continue
        for sfx, extra in (('', []), ('+sf', Wsf)):
            if (sfx and not extra) or (a.sf == 'no' and sfx) or (a.sf == 'yes' and not sfx):
                continue
            Wv, _ = LS.orthonormal(list(Wl) + list(extra), ['b%d' % (j + 1) for j in range(len(Wl) + len(extra))])
            summary['variants'].append(calibrate(Wv, labf, K, kind, T, eps, ia, ib, vA, vB, a.nit, V0, beta_mu, sig2_mu, amp, kf, tag + sfx))
    summary['wall'] = time.time() - t0
    json.dump(summary, open(os.path.join(a.out, 'refit_summary.json'), 'w'), indent=1, default=float)
    log('done in %.0f s' % summary['wall'])


if __name__ == '__main__':
    main()
