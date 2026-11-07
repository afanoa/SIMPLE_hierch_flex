#!/usr/bin/env python3
"""f2r_regref.py: the latent regularising the state REFERENCES, not the particle posterior. Each state's deviation from the
node's all-particle map is projected onto the fixed latent span (flex modes + site factors, orthonormal in the band, fitted once at
frozen poses) before the cross-half FSC weighting and masking of tree_estep_v3.references; the objective then scores against
references that differ only along directions the latent says vary. Chains (references rebuilt from each chain's own labels):
  O   objective, raw references (the drifting baseline)      C   calibrated product on raw references (f2r_calib's C)
  R   objective, every deviation projected onto the span (differences only inside flex's band)
  RH  objective, projected inside flex's band, the raw FSC-weighted deviation kept beyond it
Per chain and iteration: the held-out (beta, alpha) calibration at ITS references (what the prior would still add, and the
objective's own overconfidence: the gain of the best tempered beta over beta 1), the truth readout, and for R / RH the share of
each state's raw deviation energy inside the span and the cross-half FSC of the deviation before and after the projection.
usage: f2r_regref.py --node NODE --dump DUMP --out OUT --cache planes.npy --truth KIND:ARGS [--seed 3 --maxn 15000 --nit 4]"""
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


def project(d, W, hi, mode):
    """the band projection of a deviation onto the orthonormal span W; 'hybrid' keeps the deviation beyond flex's band."""
    fw = G['fw_band']
    p = sum(np.sum(fw * (d.real * w.real + d.imag * w.imag)) * w for w in W)
    return np.where(hi, d, p) if mode == 'hybrid' else p


def band_energy(d, sel):
    return float(np.sum(G['fw_band'] * sel * (d.real ** 2 + d.imag ** 2)))


def references_reg(r, ia, ib, eps, W, kf, mode):
    """tree_estep_v3.references with the per-state deviations projected onto W; returns the references and per-state diagnostics."""
    K = r.shape[1]
    hi = np.r_[G['cell_shell'] > kf, False]
    lo = ~hi
    alls = [TE.nn_norm(*TE.nn_acc(h, np.ones(len(h))), eps, len(h)) for h in (ia, ib)]
    st = [[TE.nn_norm(*TE.nn_acc(h, r[h, c]), eps, r[h, c].sum()) for c in range(K)] for h in (ia, ib)]
    refs = [[], []]
    diag = []
    for c in range(K):
        dev = [st[j][c] - alls[j] for j in range(2)]
        devp = [project(d, W, hi, mode) for d in dev]
        f_raw = np.clip(TE.shell_fsc(dev[0], dev[1]), 0., 1.)
        f = np.clip(TE.shell_fsc(devp[0], devp[1]), 0., 1.)
        fc = np.r_[f[G['cell_shell']], 0.]
        for j in range(2):
            refs[j].append(alls[j] + TE.masked(devp[j] * fc))
        diag.append(dict(share=[band_energy(project(dev[j], W, hi, 'full'), lo) / max(band_energy(dev[j], lo), 1e-30) for j in range(2)],
                         fsc_raw=float(f_raw[1:kf].mean()), fsc=float(f[1:kf].mean()), fsc_hi_raw=float(f_raw[kf:G['kb']].mean())))
    return refs, diag


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
    G['fw_band'] = np.r_[G['cell_fw'] * (G['cell_shell'] <= G['kb']), 0.]
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
    log('input labels, basis %s, flex band to ring %d of %d: %s' % (a.basis, kf, kb, ' '.join('%d' % (lab0 == c).sum() for c in range(K))))
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
    w, wA, wB = np.zeros((n, K)), np.zeros((n, K)), np.zeros((n, K))
    for fs, es in ((ia, ib), (ib, ia)):
        f = LS.xd_fit(b[fs], Gm[fs], LS.onehot(lab0[fs], K), Om1, G['mu1'], fixed=True)
        for dst, bb, gg in ((w, b, Gm), (wA, bA, GA), (wB, bB, GB)):
            dst[es] = LS.xd_estep(bb[es], gg[es], f['pi'], f['xi'], f['Om'])[0]
    log('flex weights (q %d): mean max all / A / B %.3f %.3f %.3f' % (q, w.max(1).mean(), wA.max(1).mean(), wB.max(1).mean()))
    if l_tgt:
        rd = LS.readout(kind, T, lab0)
        log('input labels vs truth: %s' % {kk_: v for kk_, v in rd.items() if kk_ in ('nmi_combo', 'ari_combo', 'best_enrich')})
    lw = np.log(np.maximum(w, 1e-300))
    eps = 0.5 * (TE.density(ia) + TE.density(ib))
    modes = {'O': None, 'C': None, 'R': 'full', 'RH': 'hybrid'}

    def loglik(lab, mode):
        if mode is None:
            refs, _, _ = TE.references(LS.onehot(lab, K), ia, ib, eps)
            diag = None
        else:
            refs, diag = references_reg(LS.onehot(lab, K), ia, ib, eps, W, kf, mode)
        beta, sig2 = TE.ref_scale_noise([refs[1], refs[0]], [ia, ib], lab)
        lA, lB = np.zeros((n, K)), np.zeros((n, K))
        for hh_, other in ((ia, refs[1]), (ib, refs[0])):
            lA[hh_] = FP.score_band(hh_, other, beta, sig2, vA)
            lB[hh_] = FP.score_band(hh_, other, beta, sig2, vB)
        return lA, lB, diag

    metric = {'qtoy': 'nmi_combo', 'closed': 'best_enrich'}.get(kind) if l_tgt else None
    chains = list(modes)
    labs = {ch: lab0 for ch in chains}
    summary = dict(node=a.node, n=n, K=int(K), kind=kind, seed=a.seed, maxn=a.maxn, basis=a.basis, q=q, kf=kf, kb=kb, its=[])
    l_raw1 = loglik(lab0, None)
    for it in range(1, a.nit + 1):
        rec = dict(it=it)
        for ch in chains:
            lA, lB, diag = l_raw1 if (it == 1 and modes[ch] is None) else loglik(labs[ch], modes[ch])
            ltot = lA + lB
            S = {(bt, al): C2.predictive(lA, lB, wA, wB, bt, al) for bt in C2.BETAS for al in C2.ALPHAS}
            bstar, astar = max(S, key=S.get)
            bobj = max(C2.BETAS, key=lambda x: S[(x, 0.)])
            qall = FP.softmax(bstar * ltot + astar * lw) if ch == 'C' else FP.softmax(ltot)
            hl, row = FC.readout_row(qall, lab0, kind, T, l_tgt)
            row.update(beta=bstar, alpha=astar, prior_gain=(S[(bstar, astar)] - S[(bobj, 0.)]) / n,
                       beta_obj=bobj, overconf=(S[(bobj, 0.)] - S[(1., 0.)]) / n)
            if diag is not None:
                row['diag'] = diag
            rec[ch] = row
            labs[ch] = hl
        log('it %d:' % it)
        for ch in chains:
            r_ = rec[ch]
            log('   %-2s best beta %.2f (overconfidence %+.4f), prior would add %+.4f at (%.2f, %.2f) | %s'
                % (ch, r_['beta_obj'], r_['overconf'], r_['prior_gain'], r_['beta'], r_['alpha'], FC.fmt(r_)))
            if 'diag' in r_:
                log('      span share of the deviation (even/odd) %s; deviation FSC in flex band raw -> projected %s; beyond it raw %s'
                    % (' '.join('%.2f/%.2f' % tuple(d['share']) for d in r_['diag']), ' '.join('%.2f->%.2f' % (d['fsc_raw'], d['fsc']) for d in r_['diag']),
                       ' '.join('%.2f' % d['fsc_hi_raw'] for d in r_['diag'])))
        summary['its'].append(rec)
        json.dump(summary, open(os.path.join(a.out, 'f2r_regref_summary.json'), 'w'), indent=1, default=float)
    log('table (%s; overconf = held-out gain of the best tempered beta over beta 1; prior = what alpha > 0 would still add):' % (metric or 'no truth'))
    log('   it |  O: metric agree overconf prior |  C: metric agree |  R: metric agree overconf prior |  RH: metric agree overconf prior')
    for rec in summary['its']:
        def m(ch):
            return ('%.3f' % rec[ch][metric]) if metric else '  -  '
        log('   %d  | %s %.3f %+.4f %+.4f | %s %.3f | %s %.3f %+.4f %+.4f | %s %.3f %+.4f %+.4f'
            % (rec['it'], m('O'), rec['O']['agree_in'], rec['O']['overconf'], rec['O']['prior_gain'], m('C'), rec['C']['agree_in'],
               m('R'), rec['R']['agree_in'], rec['R']['overconf'], rec['R']['prior_gain'],
               m('RH'), rec['RH']['agree_in'], rec['RH']['overconf'], rec['RH']['prior_gain']))
    summary['wall'] = time.time() - t0
    json.dump(summary, open(os.path.join(a.out, 'f2r_regref_summary.json'), 'w'), indent=1, default=float)
    log('done in %.0f s' % summary['wall'])
    FC._POOL[0].shutdown()


if __name__ == '__main__':
    main()
