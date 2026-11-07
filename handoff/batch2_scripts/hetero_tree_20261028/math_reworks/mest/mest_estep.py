#!/usr/bin/env python3
"""mest_estep.py: label-level M-estimation formulations at a node's frozen stage-B poses, from the production stage-B labels, on
tree_estep_v3's cross-half machinery (site_em references: all-map + FSC-weighted state deviation inside the mask sphere, each half
scored against the other half's references; uniform state prior, beta 1, one hard draw per particle per iteration). Arms:
  B0  the tree's Gaussian state model (tree_estep_v3's D0b1): the baseline every other arm is read against
  BV  predictive E-step: each reference's own half-map noise inflates the per-shell noise of its state
      (sigma2_k = sigma2 + beta^2 nv_k <a^2 CTF^2>) with the matching log-determinant: p(y | k) marginalised over the map error
  BT  Student-t noise: nu by maximum likelihood on a grid (own-state residuals, scale refitted by IRLS), IRLS particle weights in
      the reference M-step and the Student-t log-likelihood in the E-step (the principled M-estimator for junk and outliers)
Readouts per arm and iteration as tree_estep_v3 (moved share, mean max posterior, state n_eff, sizes, truth NMI / closed
enrichment), BV's per-shell inflation ratio, BT's nu-hat, weight percentiles and the junk AUROC of its weights (qtoy truth).
usage: mest_estep.py --node NODE --dump DUMP --out OUT --cache planes.npy --truth KIND:ARGS [--maxn 15000 --nit 4 --nproc 8]"""
import argparse, json, os, sys, time
import numpy as np
from scipy.special import gammaln
sys.path.insert(0, os.environ.get('MEST_LATENT', '/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent'))
import latent_states_v4 as LS
import tree_estep_v3 as TE

G = LS.G
log = LS.log
NU_GRID = [1., 2., 3., 5., 8., 12., 20., 30., 50., 100., 200., 500., 1e4]


def references(r, ia, ib, eps, pw=None):
    """tree_estep_v3.references with per-particle weights and the per-half raw state maps returned."""
    K = r.shape[1]
    w = np.ones(G['n']) if pw is None else pw
    alls = [TE.nn_norm(*TE.nn_acc(h, w[h]), eps, w[h].sum()) for h in (ia, ib)]
    st = [[TE.nn_norm(*TE.nn_acc(h, r[h, c] * w[h]), eps, (r[h, c] * w[h]).sum()) for c in range(K)] for h in (ia, ib)]
    refs = [[], []]
    for c in range(K):
        f = np.clip(TE.shell_fsc(st[0][c] - alls[0], st[1][c] - alls[1]), 0., 1.)
        fc = np.r_[f[G['cell_shell']], 0.]
        for j in range(2):
            refs[j].append(alls[j] + TE.masked((st[j][c] - alls[j]) * fc))
    return refs, alls


def ref_noise(refs):
    """per state and 2D band shell: the noise variance of one half's reference (map units) from the halves' difference."""
    cs, fw, NC, nsh = G['cell_shell'], G['cell_fw'], G['NC'], G['nsh']
    sw3 = np.bincount(cs, fw, G['nyq'] + 1)
    out = []
    for e, o in zip(refs[0], refs[1]):
        d = np.bincount(cs, fw * np.abs(e[:NC] - o[:NC]) ** 2, G['nyq'] + 1) / (4. * np.maximum(sw3, 1.))
        out.append(np.array([d[min(s + 1, G['nyq'])] for s in range(nsh)]))
    return np.array(out)


def own_terms(pi, refs, lab):
    """per particle and shell against the particle's own state reference: aa = sum sw a Re(M P*), bb = sum sw a^2 C |M|^2."""
    shell, sw, nsh = G['shell'], G['sw'], G['nsh']
    aa = np.zeros((len(pi), nsh)); bb = np.zeros((len(pi), nsh))
    for c0 in range(0, len(pi), 400):
        p = pi[c0:c0 + 400]; a = G['cn'][p][:, None]
        m = np.zeros((len(p), len(shell)), np.complex128)
        for c in np.unique(lab[p]):
            w = lab[p] == c
            m[w] = LS.slices(refs[c], p[w])
        x = (m.real * G['P'][p].real + m.imag * G['P'][p].imag) * sw[None] * a
        y = np.abs(m) ** 2 * G['C'][p] * sw[None] * a ** 2
        nc = len(p); key = (np.arange(nc)[:, None] * nsh + shell[None, :]).ravel()
        aa[c0:c0 + nc] = np.bincount(key, x.ravel(), nc * nsh).reshape(nc, nsh); bb[c0:c0 + nc] = np.bincount(key, y.ravel(), nc * nsh).reshape(nc, nsh)
    return aa, bb


def student_fit(aa, bb, f2, cw):
    """nu, beta, sigma2 and the IRLS weights: for each nu of the grid three IRLS rounds (weighted beta / sigma2 refit), the
    marginal log-likelihood of the own-state residuals decides nu (1e4 = Gaussian)."""
    n = len(aa); m = cw.sum(); best = None
    for nu in NU_GRID:
        w = np.ones(n)
        for it in range(3):
            beta = np.maximum((w[:, None] * aa).sum(0) / np.maximum((w[:, None] * bb).sum(0), 1e-30), 1e-12)
            res = np.maximum(f2 - 2. * beta[None] * aa + beta[None] ** 2 * bb, 0.)
            sig2 = np.maximum((w[:, None] * res).sum(0) / np.maximum(n * cw, 1.), 1e-30)
            Q = (res / sig2[None]).sum(1)
            w = (nu + m) / (nu + Q)
        LL = (gammaln((nu + m) / 2.) - gammaln(nu / 2.) - 0.5 * m * np.log(nu) - 0.5 * (nu + m) * np.log1p(Q / nu)).sum() - 0.5 * n * (cw * np.log(sig2)).sum()
        if best is None or LL > best[0]:
            best = (LL, nu, beta, sig2, w, Q)
    return best


def auroc(score, pos):
    o = np.argsort(score, kind='stable'); r = np.empty(len(score)); r[o] = np.arange(1, len(score) + 1)
    npos = pos.sum(); nneg = len(pos) - npos
    return float((r[pos].sum() - npos * (npos + 1) / 2.) / max(npos * nneg, 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--node'); ap.add_argument('--dump'); ap.add_argument('--out'); ap.add_argument('--cache', default='')
    ap.add_argument('--truth', default='none'); ap.add_argument('--lp', type=float, default=15.)
    ap.add_argument('--mskdiam', type=float, default=0.); ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--maxn', type=int, default=15000); ap.add_argument('--nit', type=int, default=4)
    ap.add_argument('--arms', default='B0,BV,BT')
    a = ap.parse_args(); a.poses = 'stageB'
    os.makedirs(a.out, exist_ok=True); t0 = time.time()
    LS.setup(a)
    from mrcio import read_mrc
    from ppinfo_site import soft_sphere
    n, n0, didx = G['n'], G['n0'], G['didx']
    G['sm_c'] = np.fft.ifftshift(soft_sphere(G['B'], G['msk_crop']))
    P0 = LS.dump_states(os.path.join(a.dump, 'b_ptcl3D.txt'), didx, n0)
    kind, T = LS.load_truth(a.truth, a.node, G['prow'])
    lab0 = LS.relabel(P0); K = lab0.max() + 1
    log('production stage-B labels: %s' % ' '.join('%d' % (lab0 == c).sum() for c in range(K)))
    V0 = LS.vol_ft(read_mrc(G['vol1'])[0])
    _, _, amp = LS.scale_noise(V0); G['cn'] = amp
    eo = G['eo']; ia, ib = np.where(eo == 0)[0], np.where(eo == 1)[0]
    eps = 0.5 * (TE.density(ia) + TE.density(ib))
    shell, sw, nsh = G['shell'], G['sw'], G['nsh']; cw = np.bincount(shell, sw, nsh); f2 = G['f2']
    # <a^2 CTF^2> per shell for the predictive inflation
    m2 = np.zeros(nsh)   # the first band ring (shell 0) holds no plane component
    for s in range(nsh):
        mm = shell == s
        if mm.any(): m2[s] = (G['C'][:, mm] * G['cn'][:, None] ** 2).mean()
    tgt = {'qtoy': 'combo', 'igmfix': 'fa1', 'closed': 'closedref'}.get(kind)
    summary = dict(node=a.node, n=n, K=int(K), kind=kind, arms={})
    base = dict(sizes=[int((lab0 == c).sum()) for c in range(K)])
    if tgt is not None and len(np.unique(T[tgt])) > 1:
        rd = LS.readout(kind, T, lab0); base.update({k: v for k, v in rd.items() if k in ('nmi_combo', 'ari_combo', 'best_enrich', 'separates')})
    summary['production'] = base
    log('production: %s' % ' '.join('%s %s' % (k, ('%.3f' % v) if isinstance(v, float) else v) for k, v in base.items()))

    def state_neff(r):
        return float(np.exp(-(r * np.log(np.maximum(r, 1e-300))).sum(1)).mean())

    for arm in a.arms.split(','):
        lab = lab0.copy(); pw = None; rng = np.random.default_rng(17); rows = []
        for it in range(1, a.nit + 1):
            refs, alls = references(LS.onehot(lab, K), ia, ib, eps, pw)
            other = [refs[1], refs[0]]
            beta, sig2 = TE.ref_scale_noise(other, [ia, ib], lab)
            ll = np.zeros((n, K)); extra = {}
            if arm == 'B0':
                for h, rf in ((ia, other[0]), (ib, other[1])):
                    s0, _ = TE.score_half(h, rf, [], beta, sig2, None); ll[h] = s0
            elif arm == 'BV':
                nvk = ref_noise(refs)                                   # (K, nsh), map units of each half's reference
                for k in range(K):
                    infl = beta ** 2 * nvk[k] * m2
                    sig2k = sig2 + infl
                    if k == 0 or it == 1:
                        extra.setdefault('inflation_ratio', []).append(np.round(infl / sig2, 4).tolist())
                    for h, rf in ((ia, other[0]), (ib, other[1])):
                        s0, _ = TE.score_half(h, [rf[k]], [], beta, sig2k, None)
                        ll[h, k] = s0[:, 0] - 0.5 * (f2[h] / sig2k[None]).sum(1) - 0.5 * (cw * np.log(sig2k)).sum()
                extra['inflation_ratio_median_over_states'] = float(np.median([np.median(x) for x in extra['inflation_ratio']]))
            elif arm == 'BT':
                aa = np.zeros((n, nsh)); bb = np.zeros((n, nsh))
                for h, rf in ((ia, other[0]), (ib, other[1])):
                    aa[h], bb[h] = own_terms(h, rf, lab)
                LLt, nu, beta_t, sig2_t, w, Q = student_fit(aa, bb, f2, cw)
                m = cw.sum()
                for h, rf in ((ia, other[0]), (ib, other[1])):
                    s0, _ = TE.score_half(h, rf, [], beta_t, sig2_t, None)
                    Qk = (f2[h] / sig2_t[None]).sum(1)[:, None] - 2. * s0
                    ll[h] = -0.5 * (nu + m) * np.log1p(np.maximum(Qk, 0.) / nu)
                pw = w / np.median(w)
                extra.update(nu=nu, m=float(m), w_pct=np.round(np.percentile(pw, [1, 5, 25, 50, 75, 95, 99]), 3).tolist(), Q_over_m_pct=np.round(np.percentile(Q / m, [1, 5, 50, 95, 99]), 3).tolist())
                if kind == 'qtoy':
                    extra['junk_auroc_of_Q'] = auroc(Q, T['junk'] > 0.5); extra['junk_share'] = float((T['junk'] > 0.5).mean())
                    extra['junk_mean_w'] = float(pw[T['junk'] > 0.5].mean()); extra['nonjunk_mean_w'] = float(pw[T['junk'] <= 0.5].mean())
            lt = ll - ll.max(1, keepdims=True); r = np.exp(lt); r /= r.sum(1, keepdims=True)
            lab = (r.cumsum(1) > rng.random(n)[:, None]).argmax(1); amax = r.argmax(1)
            row = dict(it=it, moved=float((lab != lab0).mean()), mean_max_post=float(r.max(1).mean()), state_neff=state_neff(r), sizes=[int((lab == c).sum()) for c in range(K)], **extra)
            for tag, hl in (('draw', lab), ('amax', amax)):
                if tgt is not None and len(np.unique(T[tgt])) > 1:
                    rd = LS.readout(kind, T, hl, r)
                    for k in ('nmi_combo', 'ari_combo', 'best_enrich', 'mean_purity'):
                        if k in rd: row['%s_%s' % (tag, k)] = rd[k]
                if kind == 'qtoy': row['%s_eta2_snrq' % tag] = LS.eta2(hl, T['snrq'])
            rows.append(row)
            log('%-3s it %d: moved %.3f, mean max post %.3f, n_eff %.2f, sizes %s%s%s' % (arm, it, row['moved'], row['mean_max_post'], row['state_neff'], row['sizes'],
                ''.join(' %s %.3f' % (k, row[k]) for k in sorted(row) if k.startswith(('draw_', 'amax_'))),
                ''.join(' %s=%s' % (k, extra[k]) for k in ('nu', 'w_pct', 'junk_auroc_of_Q', 'junk_mean_w', 'nonjunk_mean_w', 'inflation_ratio_median_over_states') if k in extra)))
        summary['arms'][arm] = rows
    summary['wall'] = time.time() - t0
    json.dump(summary, open(os.path.join(a.out, 'mest_estep_summary.json'), 'w'), indent=1, default=float)
    log('done in %.0f s' % summary['wall'])


if __name__ == '__main__':
    main()
