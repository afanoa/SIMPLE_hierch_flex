#!/usr/bin/env python3
"""f2r_seed_labels.py: K-state seed labels for a hetero_tree node from the flex + site-factor latent at the node's INPUT poses,
to seed stage A in place of site_em's labels. The latent sufficient statistics and the orthonormal band basis are
latent_states_v4's (flex modes + the node's site factors, band to flex's 15 A); the K-state mixture (tied Omega, z marginalised)
is fitted on all particles from --restarts k-means++ starts (best log-likelihood), then the labels are made cross-half: class
Gaussians refitted with the all-particle labels on one half (fixed responsibilities) assign the other half. Writes one integer
label per node row (1..K, the node project's row order) and reports the populations, the agreement with production's seeds and,
with a truth, NMI / closed enrichment of both.
usage: f2r_seed_labels.py --node NODE --dump DUMP --out OUT --cache planes.npy --K 5 [--truth KIND:ARGS --restarts 5 --nproc 8]"""
import argparse, glob, json, os, sys, time
from concurrent.futures import ThreadPoolExecutor
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import latent_states_v4 as LS
import flexprior_estep as FP
import f2r_calib as FC

G = LS.G
log = LS.log


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--node'); ap.add_argument('--dump'); ap.add_argument('--out'); ap.add_argument('--cache', default='')
    ap.add_argument('--truth', default='none'); ap.add_argument('--mskdiam', type=float, default=0.)
    ap.add_argument('--K', type=int, required=True); ap.add_argument('--restarts', type=int, default=5)
    ap.add_argument('--nproc', type=int, default=8); ap.add_argument('--maxn', type=int, default=0); ap.add_argument('--cmp', default='')
    a = ap.parse_args()
    a.lp, a.poses = FP.LPB, 'input'
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    FC._POOL[:] = [ThreadPoolExecutor(a.nproc), a.nproc]
    LS.xd_estep, LS.xd_mstep = FC.xd_estep_par, FC.xd_mstep_fast
    LS.setup(a)
    from mrcio import read_mrc
    n, n0, didx, K = G['n'], G['n0'], G['didx'], a.K
    if n != n0:
        raise SystemExit('the seed labels need every node row (%d of %d)' % (n, n0))
    kf = int(round(G['box_A'] / FP.LPF))
    kind, T = LS.load_truth(a.truth, a.node, G['prow'])
    tgt = {'qtoy': 'combo', 'closed': 'closedref'}.get(kind)
    l_tgt = tgt is not None and len(np.unique(T[tgt])) > 1
    V0 = LS.vol_ft(read_mrc(G['vol1'])[0])
    beta_mu, sig2_mu, amp = LS.scale_noise(V0)
    G['cn'] = amp
    pcs = sorted(glob.glob(os.path.join(a.node, 'fflex', 'flex_pca_pc[0-9][0-9][0-9].mrc')))
    Vs = [FP.lowpass_cells(LS.vol_ft(read_mrc(f)[0]), kf) for f in pcs]
    names = ['pc%d' % (j + 1) for j in range(len(pcs))]
    R0 = LS.dump_states(os.path.join(a.dump, 'r0_ptcl3D.txt'), didx, n0)
    sfm, _ = LS.site_factors(a.node, R0 if R0 is not None else np.zeros(n, int))
    Vs += [FP.lowpass_cells(LS.vol_ft(m), kf) for m in sfm]
    names += ['sf%d' % (j + 1) for j in range(len(sfm))]
    W, _ = LS.orthonormal(Vs, names)
    b, Gm = LS.suff_stats(V0, W, beta_mu, sig2_mu, amp)
    G.update(b=b, Gm=Gm)
    q = b.shape[1]
    G['mu1'] = np.zeros(q)
    f1 = LS.fit_one(np.arange(n), np.eye(q) / max(np.median(np.trace(Gm, axis1=1, axis2=2)) / q, 1e-30) * 10.)
    Om1, G['mu1'] = f1['Om'], f1['xi'][0]
    log('one Gaussian (q %d): LL/particle %.4f in %d iterations' % (q, f1['LL'] / n, f1['it']))
    jobs = [('a', np.arange(n), K, 1000 + s, Om1, None) for s in range(a.restarts)]
    res = LS.run_jobs(jobs, 1)
    best = max(res, key=lambda r: r['LL'])
    log('K %d mixture: LL/particle per start %s -> best seed %d, pi %s' % (K, ' '.join('%.4f' % (r['LL'] / n) for r in res), best['seed'],
        np.round(best['pi'], 3)))
    post = LS.xd_estep(b, Gm, best['pi'], best['xi'], best['Om'])[0]
    lab_all = post.argmax(1)
    eo = G['eo']
    ia, ib = np.where(eo == 0)[0], np.where(eo == 1)[0]
    lab = np.zeros(n, int); w = np.zeros((n, K))
    for fs, es in ((ia, ib), (ib, ia)):
        f = LS.xd_fit(b[fs], Gm[fs], LS.onehot(lab_all[fs], K), Om1, G['mu1'], fixed=True)
        w[es] = LS.xd_estep(b[es], Gm[es], f['pi'], f['xi'], f['Om'])[0]
        lab[es] = w[es].argmax(1)
    log('cross-half labels: mean max posterior %.3f, agree with the all-particle argmax %.3f' % (w.max(1).mean(), (lab == lab_all).mean()))
    # the seeds in production's order of population (largest first), 1-based, as the tree numbers its seeds
    order = np.argsort(-np.bincount(lab, minlength=K))
    rank = np.empty(K, int); rank[order] = np.arange(K)
    seeds = rank[lab] + 1
    np.savetxt(os.path.join(a.out, 'seed_labels.txt'), seeds, fmt='%d')
    prod = np.loadtxt(os.path.join(a.node, 'seeds', 'seeds.txt'), dtype=int, ndmin=1)
    out = dict(node=a.node, n=n, K=K, q=q, LL=best['LL'] / n, pops=[int((seeds == c).sum()) for c in range(1, K + 1)],
               prod_pops=[int((prod == c).sum()) for c in np.unique(prod)], nmi_vs_production_seeds=LS.nmi(seeds, prod))
    log('seed populations %s (production %s); NMI with production seeds %.3f' % (out['pops'], out['prod_pops'], out['nmi_vs_production_seeds']))
    if l_tgt:
        for tag, lb in (('latent', seeds), ('production', prod)):
            rd = LS.readout(kind, T, lb)
            out[tag] = {k_: v for k_, v in rd.items() if k_ in ('nmi_combo', 'ari_combo', 'best_enrich', 'separates')}
            log('%s seeds vs truth: %s' % (tag, out[tag]))
    out['wall'] = time.time() - t0
    json.dump(out, open(os.path.join(a.out, 'seed_labels.json'), 'w'), indent=1, default=float)
    log('done in %.0f s' % out['wall'])
    FC._POOL[0].shutdown()


if __name__ == '__main__':
    main()
