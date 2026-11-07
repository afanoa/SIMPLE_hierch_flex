#!/usr/bin/env python3
"""f2r_calib.py: the flex -> refinement corner test (calib_flexprior2.py --states stageB --poses stageB --basis sites) repeated on
another subsample (--seed) or the whole node (--maxn 0), with independent chains. calib_flexprior2.py rebuilds the references of
every arm from the calibrated arm's labels, so its objective-only row after iteration 1 is one E-step at the calibrated chain's
references. Here three chains rebuild their own:
  C  q ~ exp(beta l + alpha log w), (beta, alpha) by the checkerboard held-out likelihood every iteration (calib_flexprior2's cal)
  O  the objective alone, argmax l
  S  stickiness, argmax l + kappa [s = input label], kappa set at iteration 1 so that S keeps as many input labels as C does:
     what the flex prior gives beyond pulling particles back to the labels its class Gaussians were fitted with
  P  C's iteration-1 (beta, alpha) kept for every iteration: the prior as refinement would hold it, calibrated once between stages
Flex weights (flex modes [+ site factors], class Gaussians fitted with the input labels on one half, evaluated on the other),
objective, references and readouts as in calib_flexprior2.py. The input labels and poses are the dump's b_ptcl3D.txt: stage B's,
or stage A's in a dump directory that links stage A's dump there (f2r_run.sbatch, STAGE A).
usage: f2r_calib.py --node NODE --dump DUMP --out OUT --cache planes.npy --truth KIND:ARGS [--seed 3 --maxn 15000 --nit 4
       --basis sites --cmp production_ptcl3D_dump]"""
import argparse, glob, json, os, sys, time
from concurrent.futures import ThreadPoolExecutor
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import latent_states_v4 as LS
import tree_estep_v3 as TE
import flexprior_estep as FP
import calib_flexprior2 as C2

G = LS.G
log = LS.log
KAPPAS = np.r_[0., np.geomspace(1e-3, 1e4, 400)]
_estep0 = LS.xd_estep
_POOL = [None, 1]


def xd_estep_par(b, Gm, pi, xi, Om):
    """LS.xd_estep over particle chunks in threads (numpy's batched inverse, slogdet and matmul release the GIL); per particle the
    same arithmetic."""
    ex, nth = _POOL
    n = len(b)
    if nth <= 1 or n < 100 * nth:
        return _estep0(b, Gm, pi, xi, Om)
    cuts = np.linspace(0, n, nth + 1).astype(int)
    parts = list(ex.map(lambda c: _estep0(b[c[0]:c[1]], Gm[c[0]:c[1]], pi, xi, Om), zip(cuts[:-1], cuts[1:])))
    return tuple(np.concatenate([p[j] for p in parts]) for j in range(4))


def xd_mstep_fast(r, m, Hi):
    """LS.xd_mstep with the scatter sum_ik r_ik d_ik d_ik' as one contraction (einsum's three-operand loop dominates at q ~ 50)."""
    N = len(r)
    Nk = r.sum(0) + 1e-12
    xi = np.einsum('ik,iqk->kq', r, m) / Nk[:, None]
    d = (m - xi.T[None]) * np.sqrt(r)[:, None, :]
    Om = Hi.mean(0) + np.tensordot(d, d, axes=([0, 2], [0, 2])) / N
    return Nk / N, xi, 0.5 * (Om + Om.T)


def setup_seeded(a):
    """LS.setup with its subsample drawn by default_rng(a.seed); the shared helper draws it with default_rng(3)."""
    rng0 = np.random.default_rng
    np.random.default_rng = lambda s=None: rng0(a.seed if s == 3 else s)
    try:
        LS.setup(a)
    finally:
        np.random.default_rng = rng0


def readout_row(r, lab0, kind, T, l_tgt):
    hl = r.argmax(1)
    row = dict(mean_max_post=float(r.max(1).mean()), agree_in=float((hl == lab0).mean()),
               sizes=[int((hl == c).sum()) for c in range(r.shape[1])])
    if l_tgt:
        rd = LS.readout(kind, T, hl, r)
        row.update({k_: v for k_, v in rd.items() if k_ in ('nmi_combo', 'ari_combo', 'best_enrich', 'mean_purity')})
    return hl, row


def fmt(row):
    return ' '.join('%s %.3f' % (k_, row[k_]) for k_ in ('nmi_combo', 'ari_combo', 'best_enrich', 'agree_in', 'mean_max_post')
                    if k_ in row)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--node'); ap.add_argument('--dump'); ap.add_argument('--out'); ap.add_argument('--cache', default='')
    ap.add_argument('--truth', default='none'); ap.add_argument('--mskdiam', type=float, default=0.)
    ap.add_argument('--maxn', type=int, default=15000); ap.add_argument('--seed', type=int, default=3)
    ap.add_argument('--nit', type=int, default=4); ap.add_argument('--nproc', type=int, default=4)
    ap.add_argument('--basis', default='sites', help='flex | sites: the latent of the prior, flex modes or flex modes + site factors')
    ap.add_argument('--cmp', default='', help='a ptcl3D dump whose labels are read out against the truth too (production stage B)')
    a = ap.parse_args()
    a.lp, a.poses = FP.LPB, 'stageB'
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    _POOL[:] = [ThreadPoolExecutor(a.nproc), a.nproc]
    LS.xd_estep, LS.xd_mstep = xd_estep_par, xd_mstep_fast
    setup_seeded(a)
    log('subsample: %s; %d threads' % ('the whole node' if not a.maxn else 'maxn %d, seed %d' % (a.maxn, a.seed), a.nproc))
    from mrcio import read_mrc
    from ppinfo_site import soft_sphere
    n, n0, didx = G['n'], G['n0'], G['didx']
    G['sm_c'] = np.fft.ifftshift(soft_sphere(G['B'], G['msk_crop']))
    # the sample order of setup, to split by checkerboard parity
    B, kb = G['B'], G['kb']
    kk, hh = np.meshgrid(np.r_[0:B // 2, -B // 2:0], np.arange(B // 2), indexing='ij')
    r2 = hh * hh + kk * kk
    sel = (r2 >= 3) & (r2 <= kb * (kb + 1))
    h, k = hh[sel], kk[sel]
    assert len(h) == len(G['shell'])
    vA = ((h + k) % 2) == 0
    vB = ~vA
    kf = int(round(G['box_A'] / FP.LPF))
    log('views: %d / %d samples; flex band to ring %d, objective to ring %d' % (int(vA.sum()), int(vB.sum()), kf, kb))
    fin = os.path.join(a.dump, 'b_ptcl3D.txt')
    lab0 = LS.relabel(LS.dump_states(fin, didx, n0))
    K = lab0.max() + 1
    log('input labels and poses %s, prior basis %s: %s' % (os.path.realpath(fin), a.basis,
        ' '.join('%d' % (lab0 == c).sum() for c in range(K))))
    kind, T = LS.load_truth(a.truth, a.node, G['prow'])
    tgt = {'qtoy': 'combo', 'igmfix': 'fa1', 'closed': 'closedref'}.get(kind)
    l_tgt = tgt is not None and len(np.unique(T[tgt])) > 1
    V0 = LS.vol_ft(read_mrc(G['vol1'])[0])
    beta_mu, sig2_mu, amp = LS.scale_noise(V0)
    G['cn'] = amp
    pcs = sorted(glob.glob(os.path.join(a.node, 'fflex', 'flex_pca_pc[0-9][0-9][0-9].mrc')))
    Vs = [FP.lowpass_cells(LS.vol_ft(read_mrc(f)[0]), kf) for f in pcs]
    names = ['pc%d' % (j + 1) for j in range(len(pcs))]
    if a.basis == 'sites':
        R0 = LS.dump_states(os.path.join(a.dump, 'r0_ptcl3D.txt'), didx, n0)
        sfm, _ = LS.site_factors(a.node, R0 if R0 is not None else np.zeros(n, int))
        Vs += [FP.lowpass_cells(LS.vol_ft(m), kf) for m in sfm]
        names += ['sf%d' % (j + 1) for j in range(len(sfm))]
    W, _ = LS.orthonormal(Vs, names)
    allm = np.ones(len(h), bool)
    b, Gm = C2.suff_stats_masked(V0, W, beta_mu, sig2_mu, amp, allm)
    bA, GA = C2.suff_stats_masked(V0, W, beta_mu, sig2_mu, amp, vA)
    bB, GB = C2.suff_stats_masked(V0, W, beta_mu, sig2_mu, amp, vB)
    G.update(b=b, Gm=Gm)
    q = b.shape[1]
    G['mu1'] = np.zeros(q)
    t1 = time.time()
    f1 = LS.fit_one(np.arange(n), np.eye(q) / max(np.median(np.trace(Gm, axis1=1, axis2=2)) / q, 1e-30) * 10.)
    Om1, G['mu1'] = f1['Om'], f1['xi'][0]
    log('one Gaussian (q %d): %d iterations in %.0f s' % (q, f1['it'], time.time() - t1))
    eo = G['eo']
    ia, ib = np.where(eo == 0)[0], np.where(eo == 1)[0]
    w, wA, wB = np.zeros((n, K)), np.zeros((n, K)), np.zeros((n, K))
    for fs, es in ((ia, ib), (ib, ia)):
        t1 = time.time()
        f = LS.xd_fit(b[fs], Gm[fs], LS.onehot(lab0[fs], K), Om1, G['mu1'], fixed=True)
        for dst, bb, gg in ((w, b, Gm), (wA, bA, GA), (wB, bB, GB)):
            dst[es] = LS.xd_estep(bb[es], gg[es], f['pi'], f['xi'], f['Om'])[0]
        log('class Gaussians on %d particles: %d iterations in %.0f s' % (len(fs), f['it'], time.time() - t1))
    log('flex weights: mean max all / A / B %.3f %.3f %.3f; argmax keeps the input label %.3f'
        % (w.max(1).mean(), wA.max(1).mean(), wB.max(1).mean(), (w.argmax(1) == lab0).mean()))
    np.savez(os.path.join(a.out, 'f2r_weights.npz'), w=w, wA=wA, wB=wB, lab0=lab0, prow=G['prow'])
    if l_tgt:
        rd = LS.readout(kind, T, lab0)
        log('input labels vs truth: %s' % {kk_: v for kk_, v in rd.items() if kk_ in ('nmi_combo', 'ari_combo', 'best_enrich')})
        lc = LS.dump_states(a.cmp, didx, n0) if a.cmp else None
        if lc is not None:
            rd = LS.readout(kind, T, LS.relabel(lc))
            log('--cmp labels (%s) vs truth: %s' % (a.cmp, {kk_: v for kk_, v in rd.items() if kk_ in ('nmi_combo', 'ari_combo', 'best_enrich')}))
    eps = 0.5 * (TE.density(ia) + TE.density(ib))

    def loglik(lab):
        refs, _, _ = TE.references(LS.onehot(lab, K), ia, ib, eps)
        beta, sig2 = TE.ref_scale_noise([refs[1], refs[0]], [ia, ib], lab)
        lA, lB = np.zeros((n, K)), np.zeros((n, K))
        for hh_, other in ((ia, refs[1]), (ib, refs[0])):
            lA[hh_] = FP.score_band(hh_, other, beta, sig2, vA)
            lB[hh_] = FP.score_band(hh_, other, beta, sig2, vB)
        return lA, lB

    lw = np.log(np.maximum(w, 1e-300))
    stick = LS.onehot(lab0, K)
    metric = {'qtoy': 'nmi_combo', 'closed': 'best_enrich'}.get(kind) if l_tgt else None
    summary = dict(node=a.node, n=n, K=int(K), kind=kind, seed=a.seed, maxn=a.maxn, basis=a.basis, q=q, its=[])
    labs = dict(C=lab0, O=lab0, S=lab0, P=lab0)
    kappa = None
    l1 = loglik(lab0)
    for it in range(1, a.nit + 1):
        rec = dict(it=it)
        for ch in ('C', 'O', 'S', 'P'):
            if ch == 'P' and it == 1:
                bP, aP = bstar, astar
                labs['P'], rec['P'], rec['gainP'] = labs['C'], rec['C'], rec['gain']
                continue
            lA, lB = l1 if it == 1 else loglik(labs[ch])
            ltot = lA + lB
            if ch == 'P':
                sobjP = max(C2.predictive(lA, lB, wA, wB, bt, 0.) for bt in C2.BETAS)
                rec['gainP'] = (C2.predictive(lA, lB, wA, wB, bP, aP) - sobjP) / n
                hl, rec['P'] = readout_row(FP.softmax(bP * ltot + aP * lw), lab0, kind, T, l_tgt)
            elif ch == 'C':
                S = {(bt, al): C2.predictive(lA, lB, wA, wB, bt, al) for bt in C2.BETAS for al in C2.ALPHAS}
                bstar, astar = max(S, key=S.get)
                sobj = max(S[(bt, 0.)] for bt in C2.BETAS)
                rec.update(beta=bstar, alpha=astar, gain=(S[(bstar, astar)] - sobj) / n, gain_b1=(S[(bstar, astar)] - S[(1., 0.)]) / n,
                           grid={'%g,%g' % kv: S[kv] for kv in S})
                hl, rec['C'] = readout_row(FP.softmax(bstar * ltot + astar * lw), lab0, kind, T, l_tgt)
                rec['obj@C'] = readout_row(FP.softmax(ltot), lab0, kind, T, l_tgt)[1]
                rec['prod@C'] = readout_row(FP.softmax(ltot + lw), lab0, kind, T, l_tgt)[1]
            elif ch == 'O':
                hl, rec['O'] = readout_row(FP.softmax(ltot), lab0, kind, T, l_tgt)
            else:
                if kappa is None:
                    agr = np.array([((ltot + kp * stick).argmax(1) == lab0).mean() for kp in KAPPAS])
                    kappa = float(KAPPAS[np.argmin(np.abs(agr - rec['C']['agree_in']))])
                hl, rec['S'] = readout_row(FP.softmax(ltot + kappa * stick), lab0, kind, T, l_tgt)
                rec['kappa'] = kappa
            labs[ch] = hl
        rec['C_vs_O_moved'] = float((labs['C'] != labs['O']).mean())
        log('it %d: C (beta %.2f, alpha %.2f): held-out gain per particle over the best objective-only %+.4f (over beta 1 %+.4f); '
            'C and O labels differ for %.3f; S kappa %.3g; P (beta %.2f, alpha %.2f) held-out gain %+.4f'
            % (it, bstar, astar, rec['gain'], rec['gain_b1'], rec['C_vs_O_moved'], kappa, bP, aP, rec['gainP']))
        for tag in ('C', 'obj@C', 'prod@C', 'O', 'S', 'P'):
            log('   %-6s %s' % (tag, fmt(rec[tag])))
        summary['its'].append(rec)
        json.dump(summary, open(os.path.join(a.out, 'f2r_summary.json'), 'w'), indent=1, default=float)
    log('table (%s; agree = share of input labels kept):' % (metric or 'no truth'))
    log('   it  beta alpha   gain    C      O      S    obj@C    P   gainP | agree C   O     S     P')
    for rec in summary['its']:
        mv = ['%6.3f' % rec[t][metric] if metric else '     -' for t in ('C', 'O', 'S', 'obj@C', 'P')]
        log('   %d   %.2f  %.2f %+.4f %s %+.4f | %.3f %.3f %.3f %.3f' % (rec['it'], rec['beta'], rec['alpha'], rec['gain'], ' '.join(mv),
            rec['gainP'], rec['C']['agree_in'], rec['O']['agree_in'], rec['S']['agree_in'], rec['P']['agree_in']))
    summary['wall'] = time.time() - t0
    json.dump(summary, open(os.path.join(a.out, 'f2r_summary.json'), 'w'), indent=1, default=float)
    log('done in %.0f s' % summary['wall'])
    _POOL[0].shutdown()


if __name__ == '__main__':
    main()
