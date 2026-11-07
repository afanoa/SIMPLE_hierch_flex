#!/usr/bin/env python3
"""land_e1c.py: re-analysis of a land_e1b.py run with a noise calibration from WINDOW MEANS only (no within-window spread, which the unmixing
under-disperses). Per fold (windows on lam_S, maps and lam_test from the other pixel set):
  b1 = slope of the windows' unmixed mean t on their mean lam_S   (shrinkage of the noisy coordinate)
  b2 = slope of the windows' mean lam_test on their unmixed mean t (calibration lam'(t); held-out, so unbiased)
  reliability of lam_S = b1 b2;  noise variance of lam_S = var(lam_S) (1 - b1 b2);  c_reg = that / mean(1 / I_S)
Calibration lam(t) = monotone interpolation of the pooled (mean t, mean lam_test) window points (linear extrapolation to the grid ends).
Deconvolution (land_deconv smooth cross-half fit) of the full lam with noise c / I under: c_reg + window-mean calibration (truth-free),
c_split, and the oracle (true E[lam | t], true residual c). usage: land_e1c.py RUNDIR [--nboot 20]"""
import os, sys, json, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import land_deconv as LD
ap = argparse.ArgumentParser(); ap.add_argument('run'); ap.add_argument('--nboot', type=int, default=20); ap.add_argument('--seed', type=int, default=0); a = ap.parse_args()
def log(*s): print(*s, flush=True)
Z = np.load(a.run + '/e1b_particles.npz'); lam, I, lamA, IA, lamB, IB, t, eo = Z['lam'], Z['I'], Z['lamA'], Z['IA'], Z['lamB'], Z['IB'], Z['t'], Z['eo'].astype(int)
WL, QF, tb, LWF, PWF, zsp = Z['WL'], Z['QF'], Z['tb'], Z['LWF'], Z['PWF'], Z['zsplit']; NG = 31; tg = np.arange(NG, dtype=float); cod = t >= 0
pt = np.bincount(t[cod], minlength=NG) / cod.sum()
LS = [(lamA, IA), (lamB, IB)]; rel, cr, pts = [], [], []
for f_ in range(2):
    ls_, Is_ = LS[f_]; nw = QF[f_].shape[0]; mt = QF[f_] @ tb; mls = np.array([ls_[WL[f_] == w].mean() for w in range(nw)]); mlt = LWF[f_]
    b1 = np.polyfit(mls, mt, 1)[0]; b2 = np.polyfit(mt, mlt, 1)[0]; r_ = b1 * b2; nv = ls_.var() * (1. - r_); c_ = nv / np.mean(1. / np.maximum(Is_, 1e-30))
    rel.append(r_); cr.append(c_); pts += list(zip(mt, mlt))
    log('fold %d: b1 %.2f (t per lam), b2 %.4f (lam per A), reliability %.3f, c_reg (half-pixel lam) %.3f' % (f_, b1, b2, r_, c_))
c_reg = float(np.mean(cr)); pts = np.array(sorted(pts)); xs, ys = pts[:, 0], np.maximum.accumulate(pts[:, 1])   # monotone non-decreasing
slope = np.polyfit(xs, ys, 1)[0]; lam_t = np.interp(tg, xs, ys); lo, hi = tg < xs[0], tg > xs[-1]; lam_t[lo] = ys[0] + slope * (tg[lo] - xs[0]); lam_t[hi] = ys[-1] + slope * (tg[hi] - xs[-1])
Et = np.array([lam[t == k].mean() if (t == k).any() else np.nan for k in range(NG)]); Et = np.where(np.isfinite(Et), Et, np.interp(tg, tg[np.isfinite(Et)], Et[np.isfinite(Et)]))
c_or = float(np.mean(((lam[cod] - Et[t[cod]]) ** 2) * I[cod])); c_split = float(zsp.var())
log('noise scale: c_reg %.3f | split %.3f | oracle %.3f;  oracle reliability of lam %.3f' % (c_reg, c_split, c_or, np.corrcoef(lam[cod], Et[t[cod]])[0, 1] ** 2))
log('calibration (window means) lam(t) at 0/5/10/15/20/25/30: %s | oracle %s' % (' '.join('%+.3f' % lam_t[k] for k in range(0, 31, 5)), ' '.join('%+.3f' % Et[k] for k in range(0, 31, 5))))
def lik(lt, cv): dl = lam[:, None] - lt[None]; L = np.exp(-0.5 * dl ** 2 * (I / cv)[:, None]); return L / np.maximum(L.max(1, keepdims=True), 1e-300)
res = dict(c_reg=c_reg, c_split=c_split, c_oracle=c_or, reliability=rel, lam_t=lam_t.tolist(), truth=LD.metrics(pt, pt))
for vn, lt, cv in (('reg', lam_t, c_reg), ('reg_csplit', lam_t, c_split), ('oracle', Et, c_or), ('oraclecal_creg', Et, c_reg)):
    fr = LD.fit_cv_boot(lik(lt, cv), eo, nboot=a.nboot, seed=a.seed); m_, b_ = LD.metrics(fr['p'], pt), LD.band(fr['boot'], pt)
    res[vn] = dict(m=m_, band=b_, lam=fr['lam'], p=fr['p'].tolist(), lo=np.percentile(fr['boot'], 10, 0).tolist(), hi=np.percentile(fr['boot'], 90, 0).tolist())
    log('RESULT %-14s tv %.3f low %.3f minor %.3f [%.3f %.3f] barrier %+.2f [%+.2f %+.2f] (lam %.0e) | truth low %.3f minor %.3f barrier %+.2f' % (vn, m_['tv'], m_['low'], m_['minor'],
        b_['minor'][0], b_['minor'][2], m_['barrier'], b_['barrier'][0], b_['barrier'][2], fr['lam'], res['truth']['low'], res['truth']['minor'], res['truth']['barrier']))
json.dump(res, open(a.run + '/e1c.json', 'w'), indent=1); log('DONE %s' % a.run)
