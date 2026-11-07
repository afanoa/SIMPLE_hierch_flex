#!/usr/bin/env python3
"""flexprior_estep.py: flex's per-particle state weights as the prior of the refinement's state posterior, at frozen poses.
States are flex's own states at the node (fflex/f.simple). The flex weight w_is is the particle's posterior over those states
in flex's latent space (flex modes, band to flex's 15 A): class Gaussians with a tied covariance, z marginalised through each
particle's own precision, fitted with flex's labels on one half and evaluated on the other. The objective is the euclid
log-likelihood of the particle against cross-half state references (site_em's), started from the maps M-estimated with w.
Arms (states argmax of the posterior; references rebuilt from each arm's own labels after the first E-step):
  O15  objective, shells to 15 A          O10  objective, shells to 10 A
  F    flex weights alone                  P15  w x objective to 15 A (same band: the image counted twice)
  P10  w x objective to 10 A               PS   w x objective from 15 to 10 A only (band split: disjoint data)
usage: flexprior_estep.py --node NODE --dump DUMP --out OUT --cache planes.npy --truth KIND:ARGS [--maxn 15000 --nit 4]"""
import argparse, glob, json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import latent_states_v4 as LS
import tree_estep_v3 as TE

G = LS.G
log = LS.log
LPF = 15.          # flex's band
LPB = 10.          # stage B's band


def score_band(pi, refs, beta, sig2, smask):
    """per particle and state: the tree's log-likelihood (up to a constant) over the samples smask."""
    shell, sw = G['shell'], G['sw']
    K = len(refs)
    w1 = np.where(smask, sw * (beta / sig2)[shell], 0.)
    w2 = np.where(smask, sw * (beta ** 2 / sig2)[shell], 0.)
    s0 = np.zeros((len(pi), K))
    for c0 in range(0, len(pi), 400):
        p = pi[c0:c0 + 400]
        a = G['cn'][p][:, None]
        P, C = G['P'][p], G['C'][p]
        M = np.stack([LS.slices(rf, p) for rf in refs], 1)
        s0[c0:c0 + len(p)] = (np.einsum('ikm,im->ik', M.real, P.real * w1[None] * a) + np.einsum('ikm,im->ik', M.imag, P.imag * w1[None] * a)
                              - 0.5 * np.einsum('ikm,im->ik', np.abs(M) ** 2, C * w2[None] * a ** 2))
    return s0


def lowpass_cells(V, kmax):
    """zero the cells of a half-volume beyond ring kmax (flex's band)."""
    return np.r_[np.where(G['cell_shell'] <= kmax, V[:G['NC']], 0.), 0.]


def softmax(ll):
    ll = ll - ll.max(1, keepdims=True)
    r = np.exp(ll)
    return r / r.sum(1, keepdims=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--node'); ap.add_argument('--dump'); ap.add_argument('--out'); ap.add_argument('--cache', default='')
    ap.add_argument('--truth', default='none'); ap.add_argument('--mskdiam', type=float, default=0.)
    ap.add_argument('--maxn', type=int, default=15000); ap.add_argument('--nit', type=int, default=4)
    ap.add_argument('--nproc', type=int, default=4)
    a = ap.parse_args()
    a.lp, a.poses = LPB, 'input'
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    LS.setup(a)
    from mrcio import read_mrc
    from ppinfo_site import soft_sphere
    n, n0, didx = G['n'], G['n0'], G['didx']
    G['sm_c'] = np.fft.ifftshift(soft_sphere(G['B'], G['msk_crop']))
    kf = int(round(G['box_A'] / LPF))
    low = (G['shell'] + 1) <= kf                         # samples in flex's band
    log('bands: flex / low to ring %d (%.1f A), high rings %d..%d (to %.1f A); %d low, %d high samples'
        % (kf, G['box_A'] / kf, kf + 1, G['kb'], G['box_A'] / G['kb'], int(low.sum()), int((~low).sum())))
    F0 = LS.dump_states(os.path.join(a.dump, 'f_ptcl3D.txt'), didx, n0)
    labf = LS.relabel(F0)
    K = labf.max() + 1
    log('flex states: %s' % ' '.join('%d' % (labf == c).sum() for c in range(K)))
    kind, T = LS.load_truth(a.truth, a.node, G['prow'])
    # flex's latent sufficient statistics in its own band
    V0 = LS.vol_ft(read_mrc(G['vol1'])[0])
    beta_mu, sig2_mu, amp = LS.scale_noise(V0)
    G['cn'] = amp
    pcs = sorted(glob.glob(os.path.join(a.node, 'fflex', 'flex_pca_pc[0-9][0-9][0-9].mrc')))
    W, _ = LS.orthonormal([lowpass_cells(LS.vol_ft(read_mrc(f)[0]), kf) for f in pcs], ['pc%d' % (j + 1) for j in range(len(pcs))])
    b, Gm = LS.suff_stats(V0, W, beta_mu, sig2_mu, amp)
    G.update(b=b, Gm=Gm)
    q = b.shape[1]
    G['mu1'] = np.zeros(q)
    f1 = LS.fit_one(np.arange(n), np.eye(q) / max(np.median(np.trace(Gm, axis1=1, axis2=2)) / q, 1e-30) * 10.)
    Om1, G['mu1'] = f1['Om'], f1['xi'][0]
    # flex weights: class Gaussians fitted with flex's labels on one half, the other half's posterior
    eo = G['eo']
    ia, ib = np.where(eo == 0)[0], np.where(eo == 1)[0]
    w = np.zeros((n, K))
    for fs, es in ((ia, ib), (ib, ia)):
        f = LS.xd_fit(b[fs], Gm[fs], LS.onehot(labf[fs], K), Om1, G['mu1'], fixed=True)
        w[es] = LS.xd_estep(b[es], Gm[es], f['pi'], f['xi'], f['Om'])[0]
    log('flex weights: mean max %.3f, argmax agrees with flex labels %.3f' % (w.max(1).mean(), (w.argmax(1) == labf).mean()))
    # the weighted start: cross-half references M-estimated with w
    eps = 0.5 * (TE.density(ia) + TE.density(ib))
    refs0, alls0, _ = TE.references(w, ia, ib, eps)
    beta0, sig20 = TE.ref_scale_noise([refs0[1], refs0[0]], [ia, ib], w.argmax(1))
    tgt = {'qtoy': 'combo', 'igmfix': 'fa1', 'closed': 'closedref'}.get(kind)
    arms = {'O15': (False, low), 'O10': (False, np.ones_like(low)), 'F': (True, None), 'P15': (True, low),
            'P10': (True, np.ones_like(low)), 'PS': (True, ~low)}
    summary = dict(node=a.node, n=n, K=int(K), kind=kind, arms={})
    for arm, (l_prior, smask) in arms.items():
        refs, beta, sig2 = refs0, beta0, sig20
        rows = []
        for it in range(1, (1 if smask is None else a.nit) + 1):
            if it > 1:
                refs, _, _ = TE.references(LS.onehot(lab, K), ia, ib, eps)
                beta, sig2 = TE.ref_scale_noise([refs[1], refs[0]], [ia, ib], lab)
            ll = np.zeros((n, K))
            if smask is not None:
                for h, other in ((ia, refs[1]), (ib, refs[0])):
                    ll[h] = score_band(h, other, beta, sig2, smask)
            if l_prior:
                ll += np.log(np.maximum(w, 1e-300))
            r = softmax(ll)
            lab = r.argmax(1)
            row = dict(it=it, moved_vs_flex=float((lab != labf).mean()), mean_max_post=float(r.max(1).mean()),
                       sizes=[int((lab == c).sum()) for c in range(K)])
            if tgt is not None and len(np.unique(T[tgt])) > 1:
                rd = LS.readout(kind, T, lab, r)
                row.update({k: v for k, v in rd.items() if k in ('nmi_combo', 'ari_combo', 'nmi_fa1', 'best_enrich', 'separates', 'mean_purity')})
            if kind == 'qtoy':
                row['eta2_snrq'] = LS.eta2(lab, T['snrq'])
            rows.append(row)
            log('%-3s it %d: moved vs flex %.3f, mean max post %.3f, sizes %s%s' % (arm, it, row['moved_vs_flex'], row['mean_max_post'], row['sizes'],
                ''.join(' %s %s' % (k, ('%.3f' % row[k]) if isinstance(row[k], float) else row[k]) for k in
                        ('nmi_combo', 'ari_combo', 'best_enrich', 'mean_purity', 'eta2_snrq') if k in row)))
        summary['arms'][arm] = rows
    if tgt is not None and len(np.unique(T[tgt])) > 1:
        rd = LS.readout(kind, T, labf)
        summary['flex_labels'] = {k: v for k, v in rd.items() if k in ('nmi_combo', 'ari_combo', 'best_enrich', 'separates')}
        log('flex labels themselves: %s' % summary['flex_labels'])
    summary['wall'] = time.time() - t0
    json.dump(summary, open(os.path.join(a.out, 'flexprior_summary.json'), 'w'), indent=1, default=float)
    log('done in %.0f s' % summary['wall'])


if __name__ == '__main__':
    main()
