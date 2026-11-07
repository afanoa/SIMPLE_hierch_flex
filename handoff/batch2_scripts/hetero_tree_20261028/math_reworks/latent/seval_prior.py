#!/usr/bin/env python3
"""seval_prior.py: a node's flex prior for refine=state_eval. Per particle the log-posterior over the node's stage-B states in the
latent of flex modes plus site factors (flex band 15 A, stage-B poses), class Gaussians fitted with the stage-B labels on one particle
half and evaluated on the other (calib_flexprior2's weights), for every particle of the node. Writes one 'row lw_1 .. lw_K' line per
stage-B project row (natural log): the SIMPLE_STATE_PRIOR file.
usage: seval_prior.py --node NODE --dump DUMP --cache planes.npy --out PRIOR.txt [--nproc 8]"""
import argparse, glob, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import latent_states_v4 as LS
import flexprior_estep as FP
import calib_flexprior as CF
import flex_refit2_test as FR

G = LS.G
log = LS.log


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--node'); ap.add_argument('--dump'); ap.add_argument('--cache', default=''); ap.add_argument('--out')
    ap.add_argument('--mskdiam', type=float, default=0.); ap.add_argument('--nproc', type=int, default=8)
    a = ap.parse_args()
    a.lp, a.poses, a.maxn = FP.LPB, 'stageB', 0
    if a.nproc > 1:
        FR.threaded(a.nproc)
    LS.setup(a)
    from mrcio import read_mrc
    n, n0, didx = G['n'], G['n0'], G['didx']
    kf = int(round(G['box_A'] / FP.LPF))
    st = LS.dump_states(os.path.join(a.dump, 'b_ptcl3D.txt'), didx, n0)
    states = np.unique(st)
    assert (states == np.arange(1, len(states) + 1)).all(), 'stage-B states are not 1..K: %s' % states
    lab, K = st - 1, len(states)
    V0 = LS.vol_ft(read_mrc(G['vol1'])[0])
    beta_mu, sig2_mu, amp = LS.scale_noise(V0)
    G['cn'] = amp
    pcs = sorted(glob.glob(os.path.join(a.node, 'fflex', 'flex_pca_pc[0-9][0-9][0-9].mrc')))
    Vs = [FP.lowpass_cells(LS.vol_ft(read_mrc(f)[0]), kf) for f in pcs]
    R0 = LS.dump_states(os.path.join(a.dump, 'r0_ptcl3D.txt'), didx, n0)
    sfm, _ = LS.site_factors(a.node, R0 if R0 is not None else np.zeros(n, int))
    Vs += [FP.lowpass_cells(LS.vol_ft(m), kf) for m in sfm]
    W, _ = LS.orthonormal(Vs, ['b%d' % (j + 1) for j in range(len(Vs))])
    b, Gm = CF.suff_stats_masked(V0, W, beta_mu, sig2_mu, amp, np.ones(len(G['shell']), bool))
    q = b.shape[1]
    G.update(b=b, Gm=Gm, mu1=np.zeros(q))
    f1 = LS.fit_one(np.arange(n), np.eye(q) / max(np.median(np.trace(Gm, axis1=1, axis2=2)) / q, 1e-30) * 10.)
    Om1, G['mu1'] = f1['Om'], f1['xi'][0]
    eo = G['eo']
    ia, ib = np.where(eo == 0)[0], np.where(eo == 1)[0]
    w = np.zeros((n, K))
    for fs, es in ((ia, ib), (ib, ia)):
        f = LS.xd_fit(b[fs], Gm[fs], LS.onehot(lab[fs], K), Om1, G['mu1'], fixed=True)
        w[es] = LS.xd_estep(b[es], Gm[es], f['pi'], f['xi'], f['Om'])[0]
    np.savetxt(a.out, np.c_[didx + 1, np.log(np.maximum(w, 1e-300))], fmt=['%d'] + ['%.6f'] * K)
    log('prior: %d particles, %d states, latent rank %d; its argmax agrees with stage B %.3f; mean max w %.3f -> %s'
        % (n, K, q, (w.argmax(1) == lab).mean(), w.max(1).mean(), a.out))


if __name__ == '__main__':
    main()
