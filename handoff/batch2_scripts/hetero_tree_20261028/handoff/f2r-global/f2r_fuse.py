#!/usr/bin/env python3
"""f2r_fuse.py: the two-distribution state posterior at frozen poses. Per state a learned distribution over BOTH indicators of the
state, the latent posterior (flex modes + site factors, class Gaussians fitted cross-half as in f2r_calib.py; K-1 log-ratios) and
the objective's score vector (K-1 gaps), fitted with the input labels on one particle half and applied to the other: per-state
means and a shared covariance, so each indicator is weighed by its learned within-state spread and the two are correlated through
the covariance, with no exponents. Judged as in f2r_calib.py: the posterior from the checkerboard view A predicting view B's
objective likelihood against the best tempered objective (the held-out gain), and the truth. The indicator model is refitted at
every iteration (the gaps change with the references) on the input labels; the latent weights are fitted once.
Chains (references rebuilt from each chain's own labels):
  C   calibrated product (f2r_calib's C)       O   objective only          F   fused, shared full covariance
  FI  fused, the indicators independent         OQ  the objective alone through the learned model (its recalibration)
  L   the latent alone (argmax w, no iteration)
usage: f2r_fuse.py --node NODE --dump DUMP --out OUT --cache planes.npy --truth KIND:ARGS [--seed 3 --maxn 15000 --nit 4]"""
import argparse, glob, json, os, sys, time
from concurrent.futures import ThreadPoolExecutor
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import latent_states_v4 as LS
import tree_estep_v3 as TE
import flexprior_estep as FP
import calib_flexprior2 as C2
import f2r_calib as FC

G = LS.G
log = LS.log
CLIP = 30.        # bound on the latent log-ratios (a floored weight gives -690)
MINC = 5          # training particles a state needs for its own mean
RIDGE = 1e-4


def feats(lw, l):
    K = l.shape[1]
    z = np.clip(lw[:, :K - 1] - lw[:, K - 1:K], -CLIP, CLIP)
    return np.c_[z, l[:, :K - 1] - l[:, K - 1:K]]


def cols(X, mode):
    d = X.shape[1] // 2
    return X[:, d:] if mode == 'OQ' else X


def lda_fit(X, lab, K, mode):
    """per-state means and the pooled within-state covariance of the indicator vector; FI drops the covariance between the
    latent and the objective blocks."""
    d = X.shape[1]
    mu = np.zeros((K, d)); pi = np.full(K, 1e-9)
    S = np.zeros((d, d)); n = 0
    gm = X.mean(0)
    for s in range(K):
        m = lab == s
        if m.sum() < MINC:
            mu[s] = gm
            continue
        mu[s] = X[m].mean(0); pi[s] = m.sum()
        D = X[m] - mu[s]
        S += D.T @ D; n += int(m.sum())
    pi /= pi.sum()
    S /= max(n - K, 1)
    S += RIDGE * np.trace(S) / d * np.eye(d)
    if mode == 'FI':
        h = d // 2
        S[:h, h:] = 0.; S[h:, :h] = 0.
    return mu, pi, S


def lda_post(X, mu, pi, S):
    Si = np.linalg.inv(S)
    D = X[:, None, :] - mu[None]
    return FP.softmax(np.log(pi)[None] - 0.5 * np.einsum('nkd,de,nke->nk', D, Si, D))


def latent_share(mu, S):
    """K = 2: the share of the standardised discriminant carried by the latent block."""
    c = np.linalg.solve(S, mu[1] - mu[0]) * np.sqrt(np.diag(S))
    h = len(c) // 2
    return float(np.abs(c[:h]).sum() / max(np.abs(c).sum(), 1e-30))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--node'); ap.add_argument('--dump'); ap.add_argument('--out'); ap.add_argument('--cache', default='')
    ap.add_argument('--truth', default='none'); ap.add_argument('--mskdiam', type=float, default=0.)
    ap.add_argument('--maxn', type=int, default=15000); ap.add_argument('--seed', type=int, default=3)
    ap.add_argument('--nit', type=int, default=4); ap.add_argument('--nproc', type=int, default=4)
    ap.add_argument('--basis', default='sites'); ap.add_argument('--cmp', default='')
    a = ap.parse_args()
    a.lp, a.poses = FP.LPB, 'stageB'
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    FC._POOL[:] = [ThreadPoolExecutor(a.nproc), a.nproc]
    LS.xd_estep, LS.xd_mstep = FC.xd_estep_par, FC.xd_mstep_fast
    FC.setup_seeded(a)
    log('subsample: %s; %d threads' % ('the whole node' if not a.maxn else 'maxn %d, seed %d' % (a.maxn, a.seed), a.nproc))
    from mrcio import read_mrc
    from ppinfo_site import soft_sphere
    n, n0, didx = G['n'], G['n0'], G['didx']
    G['sm_c'] = np.fft.ifftshift(soft_sphere(G['B'], G['msk_crop']))
    B, kb = G['B'], G['kb']
    kk, hh = np.meshgrid(np.r_[0:B // 2, -B // 2:0], np.arange(B // 2), indexing='ij')
    r2 = hh * hh + kk * kk
    sel = (r2 >= 3) & (r2 <= kb * (kb + 1))
    h, k = hh[sel], kk[sel]
    vA = ((h + k) % 2) == 0
    vB = ~vA
    kf = int(round(G['box_A'] / FP.LPF))
    lab0 = LS.relabel(LS.dump_states(os.path.join(a.dump, 'b_ptcl3D.txt'), didx, n0))
    K = lab0.max() + 1
    log('input labels, prior basis %s: %s' % (a.basis, ' '.join('%d' % (lab0 == c).sum() for c in range(K))))
    if K < 2:
        raise SystemExit('one state: nothing to fuse')
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
    f1 = LS.fit_one(np.arange(n), np.eye(q) / max(np.median(np.trace(Gm, axis1=1, axis2=2)) / q, 1e-30) * 10.)
    Om1, G['mu1'] = f1['Om'], f1['xi'][0]
    eo = G['eo']
    ia, ib = np.where(eo == 0)[0], np.where(eo == 1)[0]
    halves = ((ia, ib), (ib, ia))
    w, wA, wB = np.zeros((n, K)), np.zeros((n, K)), np.zeros((n, K))
    for fs, es in halves:
        f = LS.xd_fit(b[fs], Gm[fs], LS.onehot(lab0[fs], K), Om1, G['mu1'], fixed=True)
        for dst, bb, gg in ((w, b, Gm), (wA, bA, GA), (wB, bB, GB)):
            dst[es] = LS.xd_estep(bb[es], gg[es], f['pi'], f['xi'], f['Om'])[0]
    log('flex weights (q %d): mean max all / A / B %.3f %.3f %.3f; argmax keeps the input label %.3f'
        % (q, w.max(1).mean(), wA.max(1).mean(), wB.max(1).mean(), (w.argmax(1) == lab0).mean()))
    if l_tgt:
        rd = LS.readout(kind, T, lab0)
        log('input labels vs truth: %s' % {kk_: v for kk_, v in rd.items() if kk_ in ('nmi_combo', 'ari_combo', 'best_enrich')})
    lw, lwA, lwB = (np.log(np.maximum(x, 1e-300)) for x in (w, wA, wB))
    eps = 0.5 * (TE.density(ia) + TE.density(ib))

    def loglik(lab):
        refs, _, _ = TE.references(LS.onehot(lab, K), ia, ib, eps)
        beta, sig2 = TE.ref_scale_noise([refs[1], refs[0]], [ia, ib], lab)
        lA, lB = np.zeros((n, K)), np.zeros((n, K))
        for hh_, other in ((ia, refs[1]), (ib, refs[0])):
            lA[hh_] = FP.score_band(hh_, other, beta, sig2, vA)
            lB[hh_] = FP.score_band(hh_, other, beta, sig2, vB)
        return lA, lB

    def fused(X, mode):
        """the cross-half posterior from the indicator model and, for K = 2, the latent's share of its discriminant."""
        qq = np.zeros((n, K)); sh = []
        for fs, es in halves:
            mu, pi, S = lda_fit(cols(X[fs], mode), lab0[fs], K, mode)
            qq[es] = lda_post(cols(X[es], mode), mu, pi, S)
            if K == 2 and mode != 'OQ':
                sh.append(latent_share(mu, S))
        return qq, (float(np.mean(sh)) if sh else None)

    def pred_q(qq, l2):
        return float((C2.lse(np.log(np.maximum(qq, 1e-300)) + l2) - l2.max(1)).sum())

    metric = {'qtoy': 'nmi_combo', 'closed': 'best_enrich'}.get(kind) if l_tgt else None
    chains = ['C', 'O', 'F', 'FI', 'OQ']
    labs = {ch: lab0 for ch in chains}
    summary = dict(node=a.node, n=n, K=int(K), kind=kind, seed=a.seed, maxn=a.maxn, basis=a.basis, q=q, its=[])
    l1 = loglik(lab0)
    for it in range(1, a.nit + 1):
        rec = dict(it=it)
        for ch in chains:
            lA, lB = l1 if it == 1 else loglik(labs[ch])
            ltot = lA + lB
            sobj = max(C2.predictive(lA, lB, wA, wB, bt, 0.) for bt in C2.BETAS)
            share = None
            if ch == 'C':
                S = {(bt, al): C2.predictive(lA, lB, wA, wB, bt, al) for bt in C2.BETAS for al in C2.ALPHAS}
                bstar, astar = max(S, key=S.get)
                rec.update(beta=bstar, alpha=astar)
                gain = (S[(bstar, astar)] - sobj) / n
                qall = FP.softmax(bstar * ltot + astar * lw)
            elif ch == 'O':
                gain = 0.
                qall = FP.softmax(ltot)
            else:
                qA, _ = fused(feats(lwA, lA), ch)
                qB, _ = fused(feats(lwB, lB), ch)
                gain = (pred_q(qA, lB) + pred_q(qB, lA) - sobj) / n
                qall, share = fused(feats(lw, ltot), ch)
            if it == 1 and ch == 'O':
                gL = (pred_q(wA, lB) + pred_q(wB, lA) - sobj) / n
                _, rec['L'] = FC.readout_row(w, lab0, kind, T, l_tgt)
                rec['L']['gain'] = gL
            hl, row = FC.readout_row(qall, lab0, kind, T, l_tgt)
            row['gain'] = gain
            if share is not None:
                row['latent_share'] = share
            rec[ch] = row
            labs[ch] = hl
        log('it %d: held-out gain per particle over the best tempered objective: C (beta %.2f alpha %.2f) %+.4f, F %+.4f, FI %+.4f, OQ %+.4f%s%s'
            % (it, rec['beta'], rec['alpha'], rec['C']['gain'], rec['F']['gain'], rec['FI']['gain'], rec['OQ']['gain'],
               ', L %+.4f' % rec['L']['gain'] if it == 1 else '',
               '; latent share of the fused discriminant F %.2f FI %.2f' % (rec['F']['latent_share'], rec['FI']['latent_share'])
               if 'latent_share' in rec['F'] else ''))
        for tag in (['L'] if it == 1 else []) + chains:
            log('   %-3s %s' % (tag, FC.fmt(rec[tag])))
        summary['its'].append(rec)
        json.dump(summary, open(os.path.join(a.out, 'f2r_fuse_summary.json'), 'w'), indent=1, default=float)
    log('table (%s; gain in nats per particle over the best tempered objective; agree = share of input labels kept):' % (metric or 'no truth'))
    log('   it   gain C      F     FI     OQ |  C      O      F     FI     OQ   | agree  C     O     F     FI    OQ')
    for rec in summary['its']:
        mv = ['%6.3f' % rec[t][metric] if metric else '     -' for t in chains]
        log('   %d  %+.4f %+.4f %+.4f %+.4f | %s | %s' % (rec['it'], rec['C']['gain'], rec['F']['gain'], rec['FI']['gain'],
            rec['OQ']['gain'], ' '.join(mv), ' '.join('%.3f' % rec[t]['agree_in'] for t in chains)))
    r1 = summary['its'][0]
    log('latent alone: gain %+.4f%s' % (r1['L']['gain'], (', %s %.3f' % (metric, r1['L'][metric])) if metric else ''))
    summary['wall'] = time.time() - t0
    json.dump(summary, open(os.path.join(a.out, 'f2r_fuse_summary.json'), 'w'), indent=1, default=float)
    log('done in %.0f s' % summary['wall'])
    FC._POOL[0].shutdown()


if __name__ == '__main__':
    main()
