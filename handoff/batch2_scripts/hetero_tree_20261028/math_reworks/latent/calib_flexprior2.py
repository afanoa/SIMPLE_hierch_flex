#!/usr/bin/env python3
"""calib_flexprior.py: the flex prior's exponent (and the objective's temperature) set by held-out prediction within each
particle. Each image's Fourier samples are split by checkerboard parity of (h, k) into two views with independent noise. From
view A the state posterior q_A(s) ~ exp(beta l_A(s) + alpha log w_A(s)) (objective and flex weights from A alone) predicts
view B's likelihood, sum_i log sum_s q_A(s) p(y_B | s), and symmetrically; the (beta, alpha) maximising it assigns the states
from all samples, recalibrated every iteration (references rebuilt from the assignment). States are flex's own; poses frozen;
flex weights as in flexprior_estep.py (class Gaussians in flex's latent fitted with flex's labels on the other particle half).
usage: calib_flexprior.py --node NODE --dump DUMP --out OUT --cache planes.npy --truth KIND:ARGS [--maxn 15000 --nit 4]"""
import argparse, glob, json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import latent_states_v4 as LS
import tree_estep_v3 as TE
import flexprior_estep as FP

G = LS.G
log = LS.log
BETAS  = [0.25, 0.5, 0.75, 1.0]
ALPHAS = [0., 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]


def suff_stats_masked(V0, Ws, beta, sig2, amp, smask):
    n, q = G['n'], len(Ws)
    w1 = np.where(smask, G['sw'] * (beta / sig2)[G['shell']], 0.)
    w2 = np.where(smask, G['sw'] * (beta ** 2 / sig2)[G['shell']], 0.)
    b = np.zeros((n, q)); Gm = np.zeros((n, q, q))
    for c0 in range(0, n, 400):
        sl = slice(c0, min(n, c0 + 400))
        a = amp[sl][:, None]
        S = np.stack([LS.slices(V, sl) for V in Ws], 1)
        s0 = LS.slices(V0, sl)
        T = G['C'][sl] * w2[None] * a ** 2
        Rr = G['P'][sl] * (w1[None] * a) - s0 * T
        b[sl] = np.einsum('iqm,im->iq', S.real, Rr.real) + np.einsum('iqm,im->iq', S.imag, Rr.imag)
        Gm[sl] = np.matmul(S.real * T[:, None], S.real.transpose(0, 2, 1)) + np.matmul(S.imag * T[:, None], S.imag.transpose(0, 2, 1))
    return b, Gm


def lse(x, axis=1):
    m = x.max(axis, keepdims=True)
    return (m + np.log(np.exp(x - m).sum(axis, keepdims=True))).squeeze(axis)


def predictive(lA, lB, wA, wB, beta, alpha):
    """sum over particles of log sum_s q_A(s) exp(l_B(s)) and the swap (constants in s dropped)."""
    tot = 0.
    for l1, w1, l2 in ((lA, wA, lB), (lB, wB, lA)):
        z = beta * l1 + alpha * np.log(np.maximum(w1, 1e-300))
        lq = z - lse(z)[:, None]
        tot += float((lse(lq + l2) - l2.max(1)).sum())
    return tot


def tgt_kind_ready(kind, T):
    tgt = {'qtoy': 'combo', 'igmfix': 'fa1', 'closed': 'closedref'}.get(kind)
    return tgt is not None and len(np.unique(T[tgt])) > 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--node'); ap.add_argument('--dump'); ap.add_argument('--out'); ap.add_argument('--cache', default='')
    ap.add_argument('--truth', default='none'); ap.add_argument('--mskdiam', type=float, default=0.)
    ap.add_argument('--maxn', type=int, default=15000); ap.add_argument('--nit', type=int, default=4)
    ap.add_argument('--nproc', type=int, default=4)
    ap.add_argument('--states', default='flex', help='flex | stageB: whose states (and labels) the posterior is over')
    ap.add_argument('--poses', default='input', help='input | stageB')
    ap.add_argument('--basis', default='flex', help='flex | sites: the latent of the prior, flex modes or flex modes + site factors')
    a = ap.parse_args()
    a.lp = FP.LPB
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    LS.setup(a)
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
    F0 = LS.dump_states(os.path.join(a.dump, 'f_ptcl3D.txt' if a.states == 'flex' else 'b_ptcl3D.txt'), didx, n0)
    labf = LS.relabel(F0)
    log('states %s, poses %s, prior basis %s: %s' % (a.states, a.poses, a.basis, ' '.join('%d' % (labf == c).sum() for c in range(labf.max() + 1))))
    K = labf.max() + 1
    kind, T = LS.load_truth(a.truth, a.node, G['prow'])
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
    b, Gm = suff_stats_masked(V0, W, beta_mu, sig2_mu, amp, allm)
    bA, GA = suff_stats_masked(V0, W, beta_mu, sig2_mu, amp, vA)
    bB, GB = suff_stats_masked(V0, W, beta_mu, sig2_mu, amp, vB)
    G.update(b=b, Gm=Gm)
    q = b.shape[1]
    G['mu1'] = np.zeros(q)
    f1 = LS.fit_one(np.arange(n), np.eye(q) / max(np.median(np.trace(Gm, axis1=1, axis2=2)) / q, 1e-30) * 10.)
    Om1, G['mu1'] = f1['Om'], f1['xi'][0]
    eo = G['eo']
    ia, ib = np.where(eo == 0)[0], np.where(eo == 1)[0]
    w, wA, wB = np.zeros((n, K)), np.zeros((n, K)), np.zeros((n, K))
    for fs, es in ((ia, ib), (ib, ia)):
        f = LS.xd_fit(b[fs], Gm[fs], LS.onehot(labf[fs], K), Om1, G['mu1'], fixed=True)
        for dst, bb, gg in ((w, b, Gm), (wA, bA, GA), (wB, bB, GB)):
            dst[es] = LS.xd_estep(bb[es], gg[es], f['pi'], f['xi'], f['Om'])[0]
    log('flex weights: mean max all / A / B %.3f %.3f %.3f' % (w.max(1).mean(), wA.max(1).mean(), wB.max(1).mean()))
    eps = 0.5 * (TE.density(ia) + TE.density(ib))
    # stage-B states start from their own labels' references (the tree's), flex states from the weighted maps
    r_init = LS.onehot(labf, K) if a.states == 'stageB' else w
    refs, _, _ = TE.references(r_init, ia, ib, eps)
    beta, sig2 = TE.ref_scale_noise([refs[1], refs[0]], [ia, ib], r_init.argmax(1))
    if tgt_kind_ready(kind, T):
        rd = LS.readout(kind, T, labf)
        log('the states\' own labels: %s' % {kk_: v for kk_, v in rd.items() if kk_ in ('nmi_combo', 'ari_combo', 'best_enrich', 'separates')})
    tgt = {'qtoy': 'combo', 'igmfix': 'fa1', 'closed': 'closedref'}.get(kind)
    summary = dict(node=a.node, n=n, K=int(K), kind=kind, its=[])
    lab = None
    for it in range(1, a.nit + 1):
        if it > 1:
            refs, _, _ = TE.references(LS.onehot(lab, K), ia, ib, eps)
            beta, sig2 = TE.ref_scale_noise([refs[1], refs[0]], [ia, ib], lab)
        lA, lB = np.zeros((n, K)), np.zeros((n, K))
        for hh_, other in ((ia, refs[1]), (ib, refs[0])):
            lA[hh_] = FP.score_band(hh_, other, beta, sig2, vA)
            lB[hh_] = FP.score_band(hh_, other, beta, sig2, vB)
        S = {(bt, al): predictive(lA, lB, wA, wB, bt, al) for bt in BETAS for al in ALPHAS}
        bstar, astar = max(S, key=S.get)
        ref = S[(1.0, 0.)]
        log('it %d: held-out prediction per particle relative to objective-only at beta 1: best (beta %.2f, alpha %.2f) %+.4f; '
            'product (1, 1) %+.4f; objective (%.2f, 0) %+.4f' % (it, bstar, astar, (S[(bstar, astar)] - ref) / n, (S[(1.0, 1.0)] - ref) / n,
            max(BETAS, key=lambda x: S[(x, 0.)]), (max(S[(x, 0.)] for x in BETAS) - ref) / n))
        ltot = lA + lB
        rows = {}
        for tag, (bt, al) in (('cal', (bstar, astar)), ('obj', (1.0, 0.)), ('prod', (1.0, 1.0))):
            r = FP.softmax(bt * ltot + al * np.log(np.maximum(w, 1e-300)))
            hl = r.argmax(1)
            row = dict(beta=bt, alpha=al, mean_max_post=float(r.max(1).mean()), sizes=[int((hl == c).sum()) for c in range(K)])
            if tgt is not None and len(np.unique(T[tgt])) > 1:
                rd = LS.readout(kind, T, hl, r)
                row.update({kk_: v for kk_, v in rd.items() if kk_ in ('nmi_combo', 'ari_combo', 'best_enrich', 'mean_purity')})
            if kind == 'qtoy':
                row['eta2_snrq'] = LS.eta2(hl, T['snrq'])
            rows[tag] = row
            if tag == 'cal':
                lab = hl
            log('   %-4s (beta %.2f alpha %.2f): mean max post %.3f%s' % (tag, bt, al, row['mean_max_post'],
                ''.join(' %s %.3f' % (kk_, row[kk_]) for kk_ in ('nmi_combo', 'ari_combo', 'best_enrich', 'mean_purity', 'eta2_snrq') if kk_ in row)))
        summary['its'].append(dict(it=it, beta=bstar, alpha=astar, grid={'%g,%g' % kv: S[kv] for kv in S}, rows=rows))
    summary['wall'] = time.time() - t0
    json.dump(summary, open(os.path.join(a.out, 'calib2_summary.json'), 'w'), indent=1, default=float)
    log('done in %.0f s' % summary['wall'])


if __name__ == '__main__':
    main()
