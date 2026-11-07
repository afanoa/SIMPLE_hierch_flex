#!/usr/bin/env python3
"""f2r_beta.py: the refinement's posterior temperature by held-out prediction over the joint (pose, state) neighbourhood that
stages A and B draw from. At each particle's stored stage-B rotation and --nperturb small random rotations about it (axis
uniform, angle |N(0, --sigma deg)|, in-plane included; shifts fixed so the cached samples stay valid) every state's cross-half
reference (tree_estep_v3.references from the stage-B labels) is scored on the two checkerboard views, for the stage-A band (to
15 A) and the stage-B band (to 10 A). For each beta the posterior from view A, q_A(r, s) ~ exp(beta l_A(r, s)) over the
neighbourhood x states, predicts view B, sum_i log sum_rs q_A(rs) exp(l_B(rs)), plus the swap; beta* maximises it. beta = 1 is
the Gaussian likelihood itself, master's prob_beta scale. Also the pose-only (within the stored state) and state-only (stored
pose) optima, the gains over 0.3 and 1, and the mean max posterior at each.
usage: f2r_beta.py --node NODE --dump DUMP --out OUT --cache planes.npy [--maxn 15000 --seed 3 --nperturb 64 --sigma 3.5]"""
import argparse, json, os, sys, time
from concurrent.futures import ThreadPoolExecutor
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import latent_states_v4 as LS
import tree_estep_v3 as TE
import flexprior_estep as FP
import f2r_calib as FC

G = LS.G
log = LS.log
BETAS = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 1.0, 1.25, 1.5, 2.0]


def lse(x):
    m = x.max(1, keepdims=True)
    return (m + np.log(np.exp(x - m).sum(1, keepdims=True)))[:, 0]


def pose_indices(R, h, k):
    """setup's slice indices for rotations R (n, 3, 3): the cell of each sample (NC when out of the box) and its conjugation."""
    B, nh, nyq, NC = G['B'], G['nh'], G['nyq'], G['NC']
    n = len(R)
    idx = np.empty((n, len(h)), np.int32); cnj = np.empty((n, len(h)), bool)
    for b0 in range(0, n, 4000):
        b1 = min(n, b0 + 4000)
        p = np.rint(h[None, :, None] * R[b0:b1, None, 0, :] + k[None, :, None] * R[b0:b1, None, 1, :]).astype(np.int32)
        c = p[..., 0] < 0
        p[c] = -p[c]
        okk = (np.abs(p) < nyq).all(-1)
        idx[b0:b1] = np.where(okk, (p[..., 2] % B) * B * nh + (p[..., 1] % B) * nh + p[..., 0], NC)
        cnj[b0:b1] = c
    return idx, cnj


def random_rotations(n, sigma_deg, rng):
    ax = rng.standard_normal((n, 3)); ax /= np.linalg.norm(ax, axis=1, keepdims=True)
    ang = np.abs(rng.standard_normal(n)) * np.radians(sigma_deg)
    Kx = np.zeros((n, 3, 3))
    Kx[:, 0, 1], Kx[:, 0, 2], Kx[:, 1, 0], Kx[:, 1, 2], Kx[:, 2, 0], Kx[:, 2, 1] = -ax[:, 2], ax[:, 1], ax[:, 2], -ax[:, 0], -ax[:, 1], ax[:, 0]
    s, c = np.sin(ang)[:, None, None], np.cos(ang)[:, None, None]
    return np.eye(3)[None] + s * Kx + (1. - c) * (Kx @ Kx), np.degrees(ang)


def score_masks(pi, refs, beta, sig2, masks):
    """per particle and state the log-likelihood (up to a constant) summed under each sample mask (masks: (4, m))."""
    shell, sw = G['shell'], G['sw']
    K = len(refs)
    w1 = sw * (beta / sig2)[shell]; w2 = sw * (beta ** 2 / sig2)[shell]
    out = np.zeros((len(pi), K, masks.shape[0]))
    for c0 in range(0, len(pi), 400):
        p = pi[c0:c0 + 400]
        a = G['cn'][p][:, None]
        P, C = G['P'][p], G['C'][p]
        M = np.stack([LS.slices(rf, p) for rf in refs], 1)
        contrib = (M.real * (P.real * w1[None] * a)[:, None] + M.imag * (P.imag * w1[None] * a)[:, None]
                   - 0.5 * (M.real ** 2 + M.imag ** 2) * (C * w2[None] * a ** 2)[:, None])
        out[c0:c0 + len(p)] = contrib @ masks.T
    return out


def predictive(LA, LB, beta):
    tot = 0.
    for l1, l2 in ((LA, LB), (LB, LA)):
        z = beta * l1
        lq = z - lse(z)[:, None]
        tot += float((lse(lq + l2) - l2.max(1)).sum())
    return tot


def mean_max_post(LA, LB, beta):
    return float(0.5 * (np.exp(beta * LA - lse(beta * LA)[:, None]).max(1).mean() + np.exp(beta * LB - lse(beta * LB)[:, None]).max(1).mean()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--node'); ap.add_argument('--dump'); ap.add_argument('--out'); ap.add_argument('--cache', default='')
    ap.add_argument('--truth', default='none'); ap.add_argument('--mskdiam', type=float, default=0.); ap.add_argument('--cmp', default='')
    ap.add_argument('--maxn', type=int, default=15000); ap.add_argument('--seed', type=int, default=3)
    ap.add_argument('--nperturb', type=int, default=64); ap.add_argument('--sigma', type=float, default=3.5)
    ap.add_argument('--nproc', type=int, default=4)
    a = ap.parse_args()
    a.lp, a.poses = FP.LPB, 'stageB'
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    FC._POOL[:] = [ThreadPoolExecutor(a.nproc), a.nproc]
    FC.setup_seeded(a)
    from ppinfo_site import read_dump, euler2m, soft_sphere
    from mrcio import read_mrc
    n, n0, didx = G['n'], G['n0'], G['didx']
    G['sm_c'] = np.fft.ifftshift(soft_sphere(G['B'], G['msk_crop']))
    B, kb = G['B'], G['kb']
    kk, hh = np.meshgrid(np.r_[0:B // 2, -B // 2:0], np.arange(B // 2), indexing='ij')
    r2 = hh * hh + kk * kk
    sel = (r2 >= 3) & (r2 <= kb * (kb + 1))
    h, k = hh[sel].astype(np.int32), kk[sel].astype(np.int32)
    vA = ((h + k) % 2) == 0
    vB = ~vA
    kf = int(round(G['box_A'] / FP.LPF))
    lo = (G['shell'] + 1) <= kf
    masks = np.array([vA & lo, vB & lo, vA, vB], float)
    bands = {'A (15 A)': (0, 1), 'B (10 A)': (2, 3)}
    lab0 = LS.relabel(LS.dump_states(os.path.join(a.dump, 'b_ptcl3D.txt'), didx, n0))
    K = lab0.max() + 1
    bd = read_dump(os.path.join(a.dump, 'b_ptcl3D.txt'))
    R0 = euler2m(np.c_[bd['e1'], bd['e2'], bd['e3']][didx])
    log('subsample %s; states %s; bands A to ring %d, B to ring %d; %d perturbations of sigma %.1f deg'
        % ('maxn %d seed %d' % (a.maxn, a.seed) if a.maxn else 'the whole node', ' '.join('%d' % (lab0 == c).sum() for c in range(K)),
           kf, kb, a.nperturb, a.sigma))
    V0 = LS.vol_ft(read_mrc(G['vol1'])[0])
    _, _, amp = LS.scale_noise(V0)
    G['cn'] = amp
    eo = G['eo']
    ia, ib = np.where(eo == 0)[0], np.where(eo == 1)[0]
    eps = 0.5 * (TE.density(ia) + TE.density(ib))
    refs, _, _ = TE.references(LS.onehot(lab0, K), ia, ib, eps)
    bsc, sig2 = TE.ref_scale_noise([refs[1], refs[0]], [ia, ib], lab0)
    P = a.nperturb + 1
    S, ang = random_rotations(P - 1, a.sigma, np.random.default_rng(a.seed))
    log('perturbation angles: median %.2f, 90%% %.2f deg' % (np.median(ang), np.percentile(ang, 90)))
    idx0, cnj0 = G['idx'], G['cnj']
    l = np.zeros((n, P, K, 4))
    t1 = time.time()
    for p in range(P):
        R = R0 if p == 0 else R0 @ S[p - 1][None]
        G['idx'], G['cnj'] = pose_indices(R, h, k)
        for hh_, other in ((ia, refs[1]), (ib, refs[0])):
            l[hh_, p] = score_masks(hh_, other, bsc, sig2, masks)
        if p == 0:
            log('stored pose scored in %.1f s; %d more poses' % (time.time() - t1, P - 1))
    G['idx'], G['cnj'] = idx0, cnj0
    log('neighbourhood scored in %.0f s' % (time.time() - t1))
    summary = dict(node=a.node, n=n, K=int(K), seed=a.seed, nperturb=a.nperturb, sigma=a.sigma, bands={})
    for bname, (ca, cb) in bands.items():
        LA = l[:, :, :, ca].reshape(n, P * K); LB = l[:, :, :, cb].reshape(n, P * K)
        at_stored = float(((l[:, :, :, ca] + l[:, :, :, cb]).reshape(n, P * K).argmax(1) // K == 0).mean())
        forms = {'joint': (LA, LB),
                 'pose only': (l[np.arange(n), :, lab0, ca], l[np.arange(n), :, lab0, cb]),
                 'state only': (l[:, 0, :, ca], l[:, 0, :, cb])}
        rec = dict(argmax_at_stored_pose=at_stored)
        log('band %s: the full-band argmax is the stored pose for %.3f of the particles' % (bname, at_stored))
        for fname, (A_, B_) in forms.items():
            curve = {bt: predictive(A_, B_, bt) / n for bt in BETAS}
            bstar = max(curve, key=curve.get)
            rec[fname] = dict(beta=bstar, curve=curve, gain_over_03=curve[bstar] - curve[0.3], gain_over_1=curve[bstar] - curve[1.0],
                              mmp={bt: mean_max_post(A_, B_, bt) for bt in (bstar, 0.3, 1.0)})
            log('   %-10s beta* %.2f; held-out gain per particle over 0.3 %+.4f, over 1 %+.4f; mean max posterior at beta* %.3f, 0.3 %.3f, 1 %.3f'
                % (fname, bstar, rec[fname]['gain_over_03'], rec[fname]['gain_over_1'], rec[fname]['mmp'][bstar], rec[fname]['mmp'][0.3], rec[fname]['mmp'][1.0]))
        summary['bands'][bname] = rec
    log('table: band | joint beta* (gain over 0.3 / over 1) | pose-only beta* | state-only beta*')
    for bname, rec in summary['bands'].items():
        log('   %-9s | %.2f (%+.4f / %+.4f) | %.2f | %.2f' % (bname, rec['joint']['beta'], rec['joint']['gain_over_03'], rec['joint']['gain_over_1'],
            rec['pose only']['beta'], rec['state only']['beta']))
    summary['wall'] = time.time() - t0
    json.dump(summary, open(os.path.join(a.out, 'f2r_beta_summary.json'), 'w'), indent=1, default=float)
    log('done in %.0f s' % summary['wall'])
    FC._POOL[0].shutdown()


if __name__ == '__main__':
    main()
