#!/usr/bin/env python3
"""tree_estep.py: does adding the node's latent model to the tree's own state E-step help?
At the node's stage-B poses, starting from the production stage-B labels, cross-half nearest-neighbour state references are
built as site_em builds them (all-map + FSC-weighted state deviation inside the mask sphere; each half scored against the other
half's references), and every particle is reassigned (one E-step, then soft EM iterations that rebuild the references) by
  E0  the tree's model: state references only (per-shell scale and noise, per-particle amplitude)
  EL  E0 plus within-state motion along the node's latent basis W (flex modes + site factors), marginalised over z ~ N(xi, Omega)
      (Omega from the residuals at the production labels); W deflated of the states' difference directions so it cannot absorb
      the split
  ELn EL without the deflation
Readouts per arm and iteration: moved fraction vs production, truth NMI / closed-reference enrichment, mean max posterior vs
purity (calibration), and on the nulls the posteriors' confidence.
usage: tree_estep.py --node NODE --dump DUMP --out OUT --cache planes.npy --truth KIND:ARGS [--maxn 15000 --nit 4 --nproc 8]"""
import argparse, json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import latent_states_v4 as LS

G = LS.G
log = LS.log


# ---------------------------------------------------------------- nearest-neighbour maps (ppinfo / site_em conventions)

def nn_acc(pidx, pw):
    P, C, cn, idx, cnj, mate, sw, NC = G['P'], G['C'], G['cn'], G['idx'], G['cnj'], G['mate'], G['sw'], G['NC']
    num = np.zeros(NC + 1, np.complex128); den = np.zeros(NC + 1)
    for b0 in range(0, len(pidx), 2000):
        pi = pidx[b0:b0 + 2000]; w = pw[b0:b0 + 2000]
        if not np.any(w > 0):
            continue
        ii = idx[pi]
        v = P[pi] * (w[:, None] * cn[pi][:, None] * sw[None, :])
        v = np.where(cnj[pi], np.conj(v), v)
        d = C[pi] * (w[:, None] * cn[pi][:, None] ** 2 * sw[None, :])
        mm = mate[pi]; half = mm >= 0; wf = np.where(half, 0.5, 1.)
        f = ii.ravel()
        num.real += np.bincount(f, (v.real * wf).ravel(), NC + 1)
        num.imag += np.bincount(f, (v.imag * wf).ravel(), NC + 1)
        den += np.bincount(f, (d * wf).ravel(), NC + 1)
        if half.any():
            fm = mm[half]; vm = v[half]
            num.real += np.bincount(fm, 0.5 * vm.real, NC + 1)
            num.imag += np.bincount(fm, -0.5 * vm.imag, NC + 1)
            den += np.bincount(fm, 0.5 * d[half], NC + 1)
    num[NC] = 0.; den[NC] = 0.
    return num, den


def density(pidx):
    _, den = nn_acc(pidx, np.ones(len(pidx)))
    d, cs = den[:G['NC']], G['cell_shell']
    eps = np.ones(G['nyq'] + 1)
    for s in range(G['nyq'] + 1):
        v = d[(cs == s) & (d > 0)]
        if len(v):
            eps[s] = np.sort(v)[(len(v) + 1) // 2 - 1] / len(pidx)
    return eps


def nn_norm(num, den, eps, wsum):
    NC = G['NC']
    return np.r_[num[:NC] / (den[:NC] + eps[G['cell_shell']] * max(wsum, 1.)), 0.]


def shell_fsc(x, y):
    fw, cs, NC = G['cell_fw'], G['cell_shell'], G['NC']
    x, y = x[:NC], y[:NC]
    sxy = np.bincount(cs, fw * (x.real * y.real + x.imag * y.imag), G['nyq'] + 1)
    sxx = np.bincount(cs, fw * np.abs(x) ** 2, G['nyq'] + 1)
    syy = np.bincount(cs, fw * np.abs(y) ** 2, G['nyq'] + 1)
    return sxy / np.maximum(np.sqrt(sxx * syy), 1e-30)


def masked(V):
    B, nh, NC = G['B'], G['nh'], G['NC']
    v = np.fft.irfftn(V[:NC].reshape(B, B, nh), s=(B, B, B)) * G['sm_c']
    return np.r_[np.fft.rfftn(v).ravel(), 0.]


def references(r, ia, ib, eps):
    """site_em's references: per half the all-map plus each state's FSC-weighted deviation inside the mask sphere."""
    K = r.shape[1]
    out = []
    alls = [nn_norm(*nn_acc(h, np.ones(len(h))), eps, len(h)) for h in (ia, ib)]
    st = [[nn_norm(*nn_acc(h, r[h, c]), eps, r[h, c].sum()) for c in range(K)] for h in (ia, ib)]
    refs = [[], []]
    fscs = []
    for c in range(K):
        f = np.clip(shell_fsc(st[0][c] - alls[0], st[1][c] - alls[1]), 0., 1.)
        fscs.append(f)
        fc = np.r_[f[G['cell_shell']], 0.]
        for j in range(2):
            refs[j].append(alls[j] + masked((st[j][c] - alls[j]) * fc))
    return refs, alls, fscs


# ---------------------------------------------------------------- scoring

def ref_scale_noise(refs_other, half_idx, lab):
    """per-shell scale of the references against the data (hard labels) and the residual noise, as nn_estep."""
    nsh, shell, sw = G['nsh'], G['shell'], G['sw']
    sa = np.zeros(nsh); sb = np.zeros(nsh); sf = np.zeros(nsh); n = 0
    for h, refs in zip(half_idx, refs_other):
        for c0 in range(0, len(h), 400):
            pi = h[c0:c0 + 400]
            a = G['cn'][pi][:, None]
            m = np.zeros((len(pi), len(shell)), np.complex128)
            for c in np.unique(lab[pi]):
                w = lab[pi] == c
                m[w] = LS.slices(refs[c], pi[w])
            aa = (m.real * G['P'][pi].real + m.imag * G['P'][pi].imag) * sw[None] * a
            bb = np.abs(m) ** 2 * G['C'][pi] * sw[None] * a ** 2
            sa += np.bincount(np.tile(shell, len(pi)), aa.ravel(), nsh)
            sb += np.bincount(np.tile(shell, len(pi)), bb.ravel(), nsh)
            sf += G['f2'][pi].sum(0)
            n += len(pi)
    beta = np.maximum(sa / np.maximum(sb, 1e-30), 1e-12)
    cw = np.bincount(shell, sw, nsh)
    sig2 = np.maximum((sf - 2. * beta * sa + beta ** 2 * sb) / np.maximum(n * cw, 1.), 1e-30)
    return beta, sig2


def score_half(pi, refs, Wd, beta, sig2, lat):
    """per particle and state: the tree's log-likelihood (up to a constant) and, with lat = (Oi, xi), the marginalised latent term."""
    shell, sw = G['shell'], G['sw']
    K = len(refs)
    w1 = sw * (beta / sig2)[shell]; w2 = sw * (beta ** 2 / sig2)[shell]
    s0 = np.zeros((len(pi), K)); sl = np.zeros((len(pi), K))
    for c0 in range(0, len(pi), 300):
        p = pi[c0:c0 + 300]
        a = G['cn'][p][:, None]
        P, C = G['P'][p], G['C'][p]
        M = np.stack([LS.slices(rf, p) for rf in refs], 1)                       # (n, K, m)
        s0[c0:c0 + len(p)] = (np.einsum('ikm,im->ik', M.real, (P.real * w1[None] * a)) + np.einsum('ikm,im->ik', M.imag, (P.imag * w1[None] * a))
                              - 0.5 * np.einsum('ikm,im->ik', np.abs(M) ** 2, C * w2[None] * a ** 2))
        if lat is not None and Wd:
            Oi, xi = lat
            S = np.stack([LS.slices(V, p) for V in Wd], 1)                       # (n, q, m), data units
            T = C * (sw / sig2[shell])[None] * a ** 2                              # model already in data units
            b0 = (np.einsum('iqm,im->iq', S.real, P.real * (sw / sig2[shell])[None] * a) +
                  np.einsum('iqm,im->iq', S.imag, P.imag * (sw / sig2[shell])[None] * a))
            cross = (np.einsum('iqm,ikm->iqk', S.real * (T * beta[shell][None])[:, None], M.real) +
                     np.einsum('iqm,ikm->iqk', S.imag * (T * beta[shell][None])[:, None], M.imag))
            Gm = np.matmul(S.real * T[:, None], S.real.transpose(0, 2, 1)) + np.matmul(S.imag * T[:, None], S.imag.transpose(0, 2, 1))
            H = Gm + Oi[None]
            Hi = np.linalg.inv(H)
            h = b0[:, :, None] - cross + (Oi @ xi)[None, :, None]                 # (n, q, K)
            sl[c0:c0 + len(p)] = 0.5 * np.einsum('iqk,iqp,ipk->ik', h, Hi, h)
    return s0, sl


def residual_stats(pi, refs, Wd, beta, sig2, lab):
    """b and G of the residual after each particle's own (production) state reference: the within-state latent statistics."""
    shell, sw = G['shell'], G['sw']
    q = len(Wd)
    b = np.zeros((len(pi), q)); Gm = np.zeros((len(pi), q, q))
    for c0 in range(0, len(pi), 300):
        p = pi[c0:c0 + 300]
        a = G['cn'][p][:, None]
        P, C = G['P'][p], G['C'][p]
        m = np.zeros((len(p), len(shell)), np.complex128)
        for c in np.unique(lab[p]):
            w = lab[p] == c
            m[w] = LS.slices(refs[c], p[w])
        S = np.stack([LS.slices(V, p) for V in Wd], 1)
        T = C * (sw / sig2[shell])[None] * a ** 2
        R = P * ((sw / sig2[shell])[None] * a) - m * (T * beta[shell][None])
        b[c0:c0 + len(p)] = np.einsum('iqm,im->iq', S.real, R.real) + np.einsum('iqm,im->iq', S.imag, R.imag)
        Gm[c0:c0 + len(p)] = np.matmul(S.real * T[:, None], S.real.transpose(0, 2, 1)) + np.matmul(S.imag * T[:, None], S.imag.transpose(0, 2, 1))
    return b, Gm


def band_ip(x, y):
    fw = np.r_[G['cell_fw'], 0.] * np.r_[(G['cell_shell'] <= G['kb']), False]
    return float(np.sum(fw * (x.real * y.real + x.imag * y.imag)))


def deflate(W, D):
    """W minus its projection on span(D) (band inner product), re-orthonormalised; returns the kept maps."""
    if D:
        M = np.array([[band_ip(a, b) for b in D] for a in D])
        Mi = np.linalg.pinv(M)
        W = [w - sum(c * d for c, d in zip(Mi @ np.array([band_ip(d, w) for d in D]), D)) for w in W]
    M = np.array([[band_ip(a, b) for b in W] for a in W])
    lam, U = np.linalg.eigh(M)
    keep = lam > 1e-6 * lam.max()
    T = U[:, keep] / np.sqrt(lam[keep])[None]
    return [sum(T[j, l] * W[j] for j in range(len(W))) for l in range(T.shape[1])]


def to_data_units(V, beta_mu):
    """a map in the mean map's units -> data units: per-cell scale of the cell's ring."""
    cs = np.clip(G['cell_shell'] - 1, 0, len(beta_mu) - 1)
    f = np.where(G['cell_shell'] <= G['kb'], beta_mu[cs], 0.)
    return np.r_[V[:G['NC']] * f, 0.]


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--node'); ap.add_argument('--dump'); ap.add_argument('--out'); ap.add_argument('--cache', default='')
    ap.add_argument('--truth', default='none'); ap.add_argument('--lp', type=float, default=15.)
    ap.add_argument('--mskdiam', type=float, default=0.); ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--maxn', type=int, default=15000); ap.add_argument('--nit', type=int, default=4)
    ap.add_argument('--sites', type=int, default=1)
    a = ap.parse_args()
    a.poses = 'stageB'
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    LS.setup(a)
    from mrcio import read_mrc
    from ppinfo_site import soft_sphere
    n, n0, didx = G['n'], G['n0'], G['didx']
    G['sm_c'] = np.fft.ifftshift(soft_sphere(G['B'], G['msk_crop']))
    P0 = LS.dump_states(os.path.join(a.dump, 'b_ptcl3D.txt'), didx, n0)
    R0 = LS.dump_states(os.path.join(a.dump, 'r0_ptcl3D.txt'), didx, n0)
    kind, T = LS.load_truth(a.truth, a.node, G['prow'])
    lab0 = LS.relabel(P0)
    K = lab0.max() + 1
    log('production stage-B labels: %s' % ' '.join('%d' % (lab0 == c).sum() for c in range(K)))
    V0 = LS.vol_ft(read_mrc(G['vol1'])[0])
    beta_mu, _, amp = LS.scale_noise(V0)
    G['cn'] = amp
    import glob
    pcs = sorted(glob.glob(os.path.join(a.node, 'fflex', 'flex_pca_pc[0-9][0-9][0-9].mrc')))
    Wmaps = [LS.vol_ft(read_mrc(f)[0]) for f in pcs]
    if a.sites:
        sfm, _ = LS.site_factors(a.node, R0 if R0 is not None else np.zeros(n, int))
        Wmaps += [LS.vol_ft(m) for m in sfm]
    Wdata = [to_data_units(V, beta_mu) for V in Wmaps]
    eo = G['eo']
    ia, ib = np.where(eo == 0)[0], np.where(eo == 1)[0]
    eps = 0.5 * (density(ia) + density(ib))
    # production references and the within-state latent model (Omega) at the production labels
    r = LS.onehot(lab0, K)
    refs, alls, fscs = references(r, ia, ib, eps)
    beta, sig2 = ref_scale_noise([refs[1], refs[0]], [ia, ib], lab0)
    D = [0.5 * ((refs[0][c] - alls[0]) + (refs[1][c] - alls[1])) for c in range(K)]
    Wd_def = deflate(Wdata, D)
    Wd_raw = deflate(Wdata, [])
    lat = {}
    for tag, Wd in (('EL', Wd_def), ('ELn', Wd_raw)):
        b = np.zeros((n, len(Wd))); Gm = np.zeros((n, len(Wd), len(Wd)))
        for h, other in ((ia, refs[1]), (ib, refs[0])):
            bh, gh = residual_stats(h, other, Wd, beta, sig2, lab0)
            b[h] = bh; Gm[h] = gh
        q = b.shape[1]
        Om0 = np.eye(q) / max(np.median(np.trace(Gm, axis1=1, axis2=2)) / q, 1e-30) * 10.
        f = LS.xd_fit(b, Gm, np.ones((n, 1)), Om0, maxit=150)
        lat[tag] = (Wd, np.linalg.inv(f['Om']), f['xi'][0])
        log('%s: q %d (of %d maps), within-state Omega eigenvalues %s, mean |xi| %.3g' % (tag, q, len(Wdata),
            np.array2string(np.sort(np.linalg.eigvalsh(f['Om']))[::-1][:6], precision=3), np.linalg.norm(f['xi'][0])))
    summary = dict(node=a.node, n=n, K=int(K), kind=kind, arms={})
    tgt = {'qtoy': 'combo', 'igmfix': 'fa1', 'closed': 'closedref'}.get(kind)
    for arm in ('E0', 'EL', 'ELn'):
        r = LS.onehot(lab0, K)
        pi_ = r.mean(0)
        rows = []
        for it in range(1, a.nit + 1):
            if it > 1:
                refs, alls, fscs = references(r, ia, ib, eps)
                hard = r.argmax(1)
                beta, sig2 = ref_scale_noise([refs[1], refs[0]], [ia, ib], hard)
            ll = np.zeros((n, K))
            for h, other in ((ia, refs[1]), (ib, refs[0])):
                L = None if arm == 'E0' else (lat[arm][1], lat[arm][2])
                Wd = [] if arm == 'E0' else lat[arm][0]
                s0, sl = score_half(h, other, Wd, beta, sig2, L)
                ll[h] = s0 + sl
            ll += np.log(np.maximum(pi_, 1e-12))[None]
            ll -= ll.max(1, keepdims=True)
            r = np.exp(ll); r /= r.sum(1, keepdims=True)
            pi_ = r.mean(0)
            hard = r.argmax(1)
            row = dict(it=it, moved=float((hard != lab0).mean()), mean_max_post=float(r.max(1).mean()),
                       sizes=[int((hard == c).sum()) for c in range(K)])
            if tgt is not None and len(np.unique(T[tgt])) > 1:
                rd = LS.readout(kind, T, hard, r)
                row.update({k: v for k, v in rd.items() if k in ('nmi_combo', 'ari_combo', 'nmi_fa1', 'ari_fa1', 'best_enrich', 'separates', 'groups', 'mean_purity', 'mean_max_post')})
            if kind == 'qtoy':
                row['eta2_snrq'] = LS.eta2(hard, T['snrq'])
            if kind == 'igmnull':
                row['eta2_na1'] = LS.eta2(hard, T['na1'])
            rows.append(row)
            log('%-3s it %d: moved %.3f, mean max posterior %.3f, sizes %s%s' % (arm, it, row['moved'], row['mean_max_post'], row['sizes'],
                ''.join(' %s %s' % (k, ('%.3f' % row[k]) if isinstance(row[k], float) else row[k]) for k in
                        ('nmi_combo', 'ari_combo', 'nmi_fa1', 'best_enrich', 'separates', 'mean_purity', 'eta2_snrq', 'eta2_na1') if k in row)))
        summary['arms'][arm] = rows
    # the production labels themselves
    base = dict(sizes=[int((lab0 == c).sum()) for c in range(K)])
    if tgt is not None and len(np.unique(T[tgt])) > 1:
        rd = LS.readout(kind, T, lab0)
        base.update({k: v for k, v in rd.items() if k in ('nmi_combo', 'ari_combo', 'nmi_fa1', 'ari_fa1', 'best_enrich', 'separates', 'groups')})
    summary['production'] = base
    log('production: %s' % ' '.join('%s %s' % (k, ('%.3f' % v) if isinstance(v, float) else v) for k, v in base.items()))
    summary['wall'] = time.time() - t0
    json.dump(summary, open(os.path.join(a.out, 'estep_summary.json'), 'w'), indent=1, default=float)
    log('done in %.0f s' % summary['wall'])


if __name__ == '__main__':
    main()
