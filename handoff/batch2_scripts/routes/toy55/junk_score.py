#!/usr/bin/env python3
"""usage: junk_score.py --sim A,B --poses PA,PB --labels LA,LB --mr MRA,MRB --nmaster N --out DIR [--cache C --ctfabs 1 --em_lim 12 --split 20 --defbins 8 --truth T --tcol junk --guard COL --guardtail COL --jtype COL]
Junk rule = the xd3 cut (junk_flags_xd3.txt); mix (per defocus bin) and xd2 are kept for comparison."""
import os, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--sim', required=True); ap.add_argument('--poses', required=True); ap.add_argument('--labels', required=True); ap.add_argument('--mr', required=True)
ap.add_argument('--nmaster', type=int, required=True); ap.add_argument('--out', required=True); ap.add_argument('--cache', default=''); ap.add_argument('--ctfabs', type=int, default=1)
ap.add_argument('--em_lim', type=float, default=12.); ap.add_argument('--crop', type=int, default=48); ap.add_argument('--mskdiam', type=float, default=230.); ap.add_argument('--smpd', type=float, default=1.464)
ap.add_argument('--kv', type=float, default=300.); ap.add_argument('--cs', type=float, default=2.7); ap.add_argument('--fraca', type=float, default=0.1); ap.add_argument('--eps', type=float, default=1.)
ap.add_argument('--defbins', type=int, default=8); ap.add_argument('--truth', default=''); ap.add_argument('--tcol', default='junk'); ap.add_argument('--guard', default=''); ap.add_argument('--tag', default='')
ap.add_argument('--split', type=float, default=20.); ap.add_argument('--jtype', default=''); ap.add_argument('--guardtail', default='')
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
def log(*s): print(*s, flush=True)
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
lam = 12.2643247 / np.sqrt(a.kv * 1e3 * (1. + a.kv * 1e3 * 0.978466e-6))
def rotm(e1, e2, e3):
    a1, a2, a3 = np.radians(e1), np.radians(e2), np.radians(e3); c1, si1, c2, si2, c3, si3 = np.cos(a1), np.sin(a1), np.cos(a2), np.sin(a2), np.cos(a3), np.sin(a3)
    Z1 = np.zeros((len(e1), 3, 3)); Z1[:, 0, 0] = c1; Z1[:, 0, 1] = -si1; Z1[:, 1, 0] = si1; Z1[:, 1, 1] = c1; Z1[:, 2, 2] = 1
    Y2 = np.zeros((len(e1), 3, 3)); Y2[:, 0, 0] = c2; Y2[:, 0, 2] = si2; Y2[:, 2, 0] = -si2; Y2[:, 2, 2] = c2; Y2[:, 1, 1] = 1
    Z3 = np.zeros((len(e1), 3, 3)); Z3[:, 0, 0] = c3; Z3[:, 0, 1] = -si3; Z3[:, 1, 0] = si3; Z3[:, 1, 1] = c3; Z3[:, 2, 2] = 1
    return Z1 @ Y2 @ Z3
SIMS = a.sim.split(','); h0 = np.frombuffer(open(SIMS[0] + '/all.mrcs', 'rb').read(1024), np.int32); n = int(h0[0]); c = a.crop
x = np.arange(n) - n // 2; R0 = a.mskdiam / 2 / a.smpd; msk2 = np.clip((R0 + 4 - np.hypot(x[:, None], x[None, :])) / 8., 0, 1).astype(np.float32); sl = np.r_[0:c // 2, n - c // 2:n]
fc = np.fft.fftfreq(c) * c / n; fy2, fx2 = np.meshgrid(fc, fc, indexing='ij'); kk0 = np.stack([fx2, fy2, np.zeros_like(fx2)], -1).reshape(-1, 3) * n
f3 = np.fft.fftfreq(c) * c / n; kx, ky, kz = np.meshgrid(f3, f3, f3, indexing='ij'); sh = np.minimum(np.round(np.sqrt(kx ** 2 + ky ** 2 + kz ** 2) * n).astype(int), c // 2).ravel()
ns = c // 2 + 1; m3 = c ** 3 + 1; kb = min(c // 2 - 1, int(round(n * a.smpd / a.em_lim)))
r2 = np.round(np.hypot(kk0[:, 0], kk0[:, 1])).astype(int); cols = np.where((r2 >= 2) & (r2 <= kb))[0]; NPIX = len(cols)
kf = (np.hypot(fx2, fy2) / a.smpd).ravel()[cols]; ang = np.arctan2(fy2, fx2).ravel()[cols]; HI = r2[cols] >= int(round(n * a.smpd / a.split))
log('box %d, crop %d, band shells 2..%d (to %.1f A), %d pixels per image' % (n, c, kb, n * a.smpd / kb, NPIX))
def load_half(simdir, posef, lab, tag):
    rows = np.where(lab > 0)[0]
    Pz = np.array([[float(kv(l).get(k, 0)) for k in ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast')] for l in open(posef) if l.strip()])[rows]
    hh = np.frombuffer(open(simdir + '/all.mrcs', 'rb').read(1024), np.int32); nz, nsym = int(hh[2]), int(hh[23])
    cache = '%s_%s_c%d.npy' % (a.cache, tag, c) if a.cache else ''; Fc = np.load(cache, mmap_mode='r') if cache and os.path.isfile(cache) else None
    I = None if Fc is not None else np.memmap(simdir + '/all.mrcs', np.float32, 'r', offset=1024 + nsym, shape=(nz, n, n))
    N = len(rows); Wv = np.zeros((N, NPIX), np.complex64); Wd = np.zeros((N, NPIX), np.float32); F2 = np.zeros((N, NPIX), np.float32); ix = np.zeros((N, NPIX), np.int32)
    for b in range(0, N, 10000):
        rr = rows[b:b + 10000]; pp = Pz[b:b + 10000]
        if Fc is not None: F = np.asarray(Fc[rr]).reshape(len(rr), -1)[:, cols]
        else: blk = np.array(I[rr], np.float32) * msk2; F = np.fft.fft2(np.fft.ifftshift(blk, axes=(1, 2)))[:, sl][:, :, sl].reshape(len(rr), -1)[:, cols]
        df = 0.5 * ((pp[:, 5] + pp[:, 6])[:, None] + (pp[:, 5] - pp[:, 6])[:, None] * np.cos(2. * (ang[None] - np.radians(pp[:, 7])[:, None]))) * 1e4
        chi = np.pi * lam * df * kf[None] ** 2 - 0.5 * np.pi * a.cs * 1e7 * lam ** 3 * kf[None] ** 4; C = (np.sqrt(1. - a.fraca ** 2) * np.sin(chi) + a.fraca * np.cos(chi)).astype(np.float32)
        if a.ctfabs: C = np.abs(C)
        ramp = np.exp(2j * np.pi * (fx2.ravel()[cols][None] * pp[:, 3, None] + fy2.ravel()[cols][None] * pp[:, 4, None])).astype(np.complex64)
        Wv[b:b + len(rr)] = C * np.conj(ramp) * F; Wd[b:b + len(rr)] = C ** 2; F2[b:b + len(rr)] = np.abs(F) ** 2
        kk = kk0[cols]; p3 = np.einsum('qj,nij->nqi', kk, rotm(pp[:, 0], pp[:, 1], pp[:, 2])); q = np.round(p3).astype(np.int32)
        ok = (np.abs(q) <= c // 2 - 1).all(2); q = q % c; ixb = q[..., 2] * c * c + q[..., 1] * c + q[..., 0]; ixb[~ok] = c ** 3; ix[b:b + len(rr)] = ixb
    log('%s: %d node particles, states %s' % (tag, N, np.bincount(lab[rows]).tolist()[1:]))
    return dict(rows=rows, Wv=Wv, Wd=Wd, F2=F2, ix=ix, shq=np.minimum(np.r_[sh, 0][ix], ns - 1).astype(np.int16), lab=lab[rows], df=0.5 * (Pz[:, 5] + Pz[:, 6]))
LABS = [np.loadtxt(f, dtype=int) for f in a.labels.split(',')]; K = int(max(l.max() for l in LABS)); MR = [np.loadtxt(f).astype(int) - 1 for f in a.mr.split(',')]
D = [load_half(SIMS[h], a.poses.split(',')[h], LABS[h], 'AB'[h]) for h in (0, 1)]
for d in D: W = np.zeros((len(d['rows']), K), np.float32); W[np.arange(len(d['rows'])), d['lab'] - 1] = 1.; d['W'] = W
neq = int(np.median([max((D[h]['W'][:, g] > 0).sum(), 1) for h in (0, 1) for g in range(K)]))
def shell_eps(d):
    den = np.bincount(d['ix'].ravel(), d['Wd'].ravel(), m3)[:-1] * (neq / len(d['W'])); return np.array([np.median(den[(sh == k) & (den > 0)]) if ((sh == k) & (den > 0)).any() else 1. for k in range(ns)]) * a.eps
eps = 0.5 * (shell_eps(D[0]) + shell_eps(D[1])); epsv = np.r_[eps[sh], eps[-1]]
def wmap(d, w):
    ii = d['ix'].ravel(); num = np.bincount(ii, (d['Wv'].real * w[:, None]).ravel(), m3) + 1j * np.bincount(ii, (d['Wv'].imag * w[:, None]).ravel(), m3)
    V = num / (np.bincount(ii, (d['Wd'] * w[:, None]).ravel(), m3) + epsv * max(float(w.sum()), 1.) / neq); V[-1] = 0.; return V
Mi = [(wmap(D[h], np.ones(len(D[h]['W']), np.float32)), [wmap(D[h], D[h]['W'][:, g]) for g in range(K)]) for h in (0, 1)]
def refs(h):
    (cA, MA), (cB, MB) = Mi; R = []
    for g in range(K):
        DA, DB = MA[g][:-1] - cA[:-1], MB[g][:-1] - cB[:-1]
        f = np.clip(np.bincount(sh, np.real(DA * np.conj(DB)), ns) / np.maximum(np.sqrt(np.bincount(sh, np.abs(DA) ** 2, ns) * np.bincount(sh, np.abs(DB) ** 2, ns)), 1e-30), 0., 1.)
        R.append(np.r_[(cB if h == 0 else cA)[:-1] + f[sh] * (DB if h == 0 else DA), 0.].astype(np.complex64))
    return R
OUT = {k: np.full(a.nmaster, np.nan) for k in ('amp', 'amph', 'den', 'q', 'chi', 'gain', 'best', 'df')}
for h in (0, 1):
    d = D[h]; R = refs(h); nn = len(d['W']); mp = d['lab'] - 1; bn = np.zeros(ns); bd = np.zeros(ns); R2 = np.zeros(ns); C2 = np.zeros(ns)
    for g in range(K):
        rr = np.where(mp == g)[0]
        for b in range(0, len(rr), 8000):
            r = rr[b:b + 8000]; Pg = R[g][d['ix'][r]]; qq = d['shq'][r].ravel()
            bn += np.bincount(qq, np.real(Pg * np.conj(d['Wv'][r])).ravel(), ns); bd += np.bincount(qq, (d['Wd'][r] * np.abs(Pg) ** 2).ravel(), ns)
    beta = np.clip(bn / np.maximum(bd, 1e-30), 0., 2.)
    for g in range(K):
        rr = np.where(mp == g)[0]
        for b in range(0, len(rr), 8000):
            r = rr[b:b + 8000]; Pg = R[g][d['ix'][r]] * beta[d['shq'][r]]; res = d['F2'][r] - 2. * np.real(Pg * np.conj(d['Wv'][r])) + d['Wd'][r] * np.abs(Pg) ** 2
            R2 += np.bincount(d['shq'][r].ravel(), res.ravel(), ns); C2 += np.bincount(d['shq'][r].ravel(), minlength=ns)
    wsh = (1. / np.maximum(R2 / np.maximum(C2, 1), 1e-12)).astype(np.float32); NUM, DEN, NUMh, DENh = (np.zeros((nn, K)) for _ in range(4)); E2 = np.zeros(nn)
    for b in range(0, nn, 8000):
        s_ = slice(b, b + 8000); qq = d['shq'][s_]; bw = beta[qq].astype(np.float32); wi = wsh[qq]; E2[s_] = (d['F2'][s_] * wi).sum(1)
        for g in range(K):
            Pg = R[g][d['ix'][s_]] * bw; t1 = np.real(Pg * np.conj(d['Wv'][s_])) * wi; t2 = d['Wd'][s_] * np.abs(Pg) ** 2 * wi
            NUM[s_, g] = t1.sum(1); DEN[s_, g] = t2.sum(1); NUMh[s_, g] = t1[:, HI].sum(1); DENh[s_, g] = t2[:, HI].sum(1)
    S = NUM - 0.5 * DEN; bg = S.argmax(1); i = np.arange(nn); num, den = NUM[i, bg], DEN[i, bg]; mr = MR[h][d['rows']]
    OUT['amp'][mr] = num / np.maximum(den, 1e-12); OUT['q'][mr] = (num / np.sqrt(np.maximum(den * E2, 1e-12))) / np.sqrt(den / (den + NPIX)); OUT['chi'][mr] = (E2 - 2 * num + den) / NPIX
    OUT['gain'][mr] = S[i, bg]; OUT['best'][mr] = bg + 1; OUT['df'][mr] = d['df']; OUT['den'][mr] = den; OUT['amph'][mr] = NUMh[i, bg] / np.maximum(DENh[i, bg], 1e-12)
    log('half %s: best state differs from the label for %.3f; amplitude median %.2f, chi2/pix median %.2f' % ('AB'[h], (bg != mp).mean(), np.median(num / den), np.median((E2 - 2 * num + den) / NPIX)))
def mix2(v, var, iters=200):
    lo, hi = np.percentile(v, [5, 50]); mu = np.array([lo, hi]); sd = np.array([v.std(), v.std()]) * 0.5 + 1e-9; pi = np.array([0.2, 0.8])
    for _ in range(iters):
        L = pi[None] * np.exp(-0.5 * ((v[:, None] - mu[None]) / sd[None]) ** 2) / sd[None]; R_ = L / np.maximum(L.sum(1, keepdims=True), 1e-300)
        pi = R_.mean(0); mu = (R_ * v[:, None]).sum(0) / np.maximum(R_.sum(0), 1e-12); sd = np.sqrt((R_ * (v[:, None] - mu[None]) ** 2).sum(0) / np.maximum(R_.sum(0), 1e-12)) + 1e-9
    return mu, sd ** 2, pi
def xdk(v, var, mu, s2, pi, iters=500):   # heteroscedastic mixture: v_i ~ N(mu_c, s2_c + var_i), var_i = 1/DEN_i known per particle (extreme deconvolution EM)
    for _ in range(iters):
        T = s2[None] + var[:, None]; L = pi[None] * np.exp(-0.5 * (v[:, None] - mu[None]) ** 2 / T) / np.sqrt(T); R_ = L / np.maximum(L.sum(1, keepdims=True), 1e-300)
        bm = mu[None] + s2[None] / T * (v[:, None] - mu[None]); Bv = s2[None] * var[:, None] / T; w = np.maximum(R_.sum(0), 1e-12)
        pi = R_.mean(0); mu = (R_ * bm).sum(0) / w; s2 = np.maximum((R_ * ((bm - mu[None]) ** 2 + Bv)).sum(0) / w, 1e-8)
    return mu, s2, pi
def fit(v, var, kind):
    p2, p5, p50 = np.percentile(v, [2, 5, 50]); vv = v.var()
    if kind == 'mix': return mix2(v, var)
    if kind == 'xd2': return xdk(v, var, np.array([p5, p50]), np.full(2, 0.25 * vv), np.array([0.2, 0.8]))
    return xdk(v, var, np.array([p2, p50, p50]), np.array([0.1, 0.25, 1.]) * vv, np.array([0.05, 0.75, 0.2]))   # junk + genuine core + broad genuine tail
def junk_set(mu, s2, pi, mv):   # genuine = the heaviest component; junk = any other >= 2 sd below it (two Gaussians are bimodal only beyond ~2 sd)
    g = int(np.argmax(pi)); return [c for c in range(len(mu)) if c != g and pi[c] >= 0.005 and (mu[g] - mu[c]) / np.sqrt(0.5 * (s2[g] + s2[c]) + mv) >= 2.]
def post(x, var, mu, s2, pi, J):
    T = s2[None] + var[:, None]; L = pi[None] * np.exp(-0.5 * (x[:, None] - mu[None]) ** 2 / T) / np.sqrt(T); return L[:, J].sum(1) / np.maximum(L.sum(1), 1e-300)
node = np.isfinite(OUT['amp']); AMP = OUT['amp']; VAR = 1. / np.maximum(OUT['den'], 1e-12)
dfb = np.full(a.nmaster, -1); qs = np.quantile(OUT['df'][node], np.linspace(0, 1, a.defbins + 1)); dfb[node] = np.clip(np.searchsorted(qs, OUT['df'][node], side='right') - 1, 0, a.defbins - 1)
def cut(fitm, kind, quiet=True):
    P = np.full(a.nmaster, np.nan)
    for b in (range(a.defbins) if kind == 'mix' else [0]):
        m = node & (dfb == b) if kind == 'mix' else node; mf = m & fitm
        if mf.sum() < 200: P[m] = 0.; continue
        mu, s2, pi = fit(AMP[mf], VAR[mf], kind); mv = 0. if kind == 'mix' else np.median(VAR[mf]); J = junk_set(mu, s2, pi, mv)
        P[m] = post(AMP[m], np.zeros(m.sum()) if kind == 'mix' else VAR[m], mu, s2, pi, J) if J else 0.
        if not quiet: log('%s (%d ptcls): components (mean / sd / weight) %s | junk %s' % ('defocus bin %d (%.2f-%.2f um)' % (b, qs[b], qs[b + 1]) if kind == 'mix' else '%s, per-particle noise 1/DEN (median noise sd %.2f; intrinsic sds)' % (kind, np.sqrt(mv)),
                             mf.sum(), '  '.join('%.2f/%.2f/%.3f' % (m_, np.sqrt(s_), w_) for m_, s_, w_ in sorted(zip(mu, s2, pi))), [round(float(mu[c]), 2) for c in J] or 'none'))
    return P
CUTS = (('mix', 'junk_flags'), ('xd2', 'junk_flags_xd'), ('xd3', 'junk_flags_xd3')); PO = {k: cut(node, k, False) for k, _ in CUTS}; FL = {k: node & (PO[k] >= 0.5) for k, _ in CUTS}
L = np.zeros(a.nmaster, int)
for h in (0, 1): L[MR[h][D[h]['rows']]] = D[h]['lab']
HA = np.zeros(a.nmaster, bool); HA[MR[0][D[0]['rows']]] = True; HB = node & ~HA
def jac(x, y): return (x & y).sum() / max((x | y).sum(), 1)
def qtl(x): return ' '.join('%.2f' % v for v in np.percentile(x, [1, 5, 25, 50, 75, 95, 99])) if len(x) else '-'
log('amplitude quantiles 1/5/25/50/75/95/99 %% per state label: ' + ' | '.join('%d: %s' % (g, qtl(AMP[node & (L == g)])) for g in range(1, K + 1)))
np.savetxt(a.out + '/junk_scores.txt', np.c_[OUT['best'], AMP, OUT['amph'], OUT['q'], OUT['chi'], OUT['gain'], OUT['den'], FL['mix'], FL['xd2'], PO['xd3'], FL['xd3']], fmt='%.4g',
           header='best_state amplitude amp_highband corr_over_expected chi2_per_pixel ll_gain den junk_mix junk_xd2 junk_post_xd3 junk_xd3 (nan outside the node)')
for k, fn in CUTS: np.savetxt(a.out + '/' + fn + '.txt', np.where(node, FL[k].astype(int), -1), fmt='%d')
if a.truth:
    hd = [l.split()[1:] for l in open(a.truth) if l.startswith('#')][0]; T = np.genfromtxt(a.truth, comments='#'); tr = T[:, 0].astype(int) - 1; tj = np.zeros(a.nmaster, bool); tj[tr] = T[:, hd.index(a.tcol)] > 0
    def col(cn): v = np.full(a.nmaster, np.nan); v[tr] = T[:, hd.index(cn)]; return v
    def auc(s, y): rk = np.argsort(np.argsort(s)) + 1.; n1 = y.sum(); n0 = len(y) - n1; return (rk[y].sum() - n1 * (n1 + 1) / 2) / max(n1 * n0, 1)
    m = node; print('  truth: %d junk of %d node particles; AUC (junk LOWER): amplitude %.3f, high-band amplitude %.3f, corr/expected %.3f, ll gain %.3f; chi2 (junk HIGHER) %.3f'
                    % ((tj & m).sum(), m.sum(), 1 - auc(AMP[m], tj[m]), 1 - auc(OUT['amph'][m], tj[m]), 1 - auc(OUT['q'][m], tj[m]), 1 - auc(OUT['gain'][m], tj[m]), auc(OUT['chi'][m], tj[m])))
    print('  amplitude quantiles: genuine %s | junk %s' % (qtl(AMP[m & ~tj]), qtl(AMP[m & tj])) + ''.join(' | junk type %g %s' % (v, qtl(AMP[m & tj & (col(a.jtype) == v)])) for v in (np.unique(col(a.jtype)[m & tj]) if a.jtype else [])))
    print('  high-band amplitude quantiles: genuine %s | junk %s' % (qtl(OUT['amph'][m & ~tj]), qtl(OUT['amph'][m & tj])) + ''.join(' | junk type %g %s' % (v, qtl(OUT['amph'][m & tj & (col(a.jtype) == v)])) for v in (np.unique(col(a.jtype)[m & tj]) if a.jtype else [])))
for nm, _ in CUTS:
    J = FL[nm]
    print('%s [%s cut]: %d node particles, %d flagged junk (%.4f; half A %.4f, half B %.4f); per state (label): %s' % (a.tag or a.out, nm, node.sum(), J.sum(), J.sum() / max(node.sum(), 1), (J & HA).sum() / max(HA.sum(), 1), (J & HB).sum() / max(HB.sum(), 1),
          ' '.join('%d: %d/%d (%.4f)' % (g, (J & (L == g)).sum(), (node & (L == g)).sum(), (J & (L == g)).sum() / max((node & (L == g)).sum(), 1)) for g in range(1, K + 1))))
    JA = node & (cut(HA, nm) >= 0.5); JB = node & (cut(HB, nm) >= 0.5)
    print('  cut fitted on half A alone / half B alone, applied to all node particles: flagged %d / %d; Jaccard A-B %.3f, A-full %.3f, B-full %.3f' % (JA.sum(), JB.sum(), jac(JA, JB), jac(JA, J), jac(JB, J)))
    if a.truth:
        print('  vs truth: precision %.3f recall %.3f; false-junk rate on genuine %.4f' % ((J & tj).sum() / max(J.sum(), 1), (J & tj).sum() / max((tj & m).sum(), 1), (J & ~tj & m).sum() / max((~tj & m).sum(), 1))
              + ''.join('; recall junk type %g %.3f' % (v, (J & tj & (col(a.jtype) == v)).sum() / max((m & tj & (col(a.jtype) == v)).sum(), 1)) for v in (np.unique(col(a.jtype)[m & tj]) if a.jtype else [])))
        for gcol in [x for x in a.guard.split(',') if x]:
            gc = col(gcol); print('  guard %s: ' % gcol + ', '.join('= %g (genuine %d) false-junk %.4f' % (v, (m & ~tj & (gc == v)).sum(), (J & m & ~tj & (gc == v)).sum() / max((m & ~tj & (gc == v)).sum(), 1)) for v in np.unique(gc[m & ~tj])))
        for gcol in [x for x in a.guardtail.split(',') if x]:
            gc = col(gcol); lo, hi = np.nanpercentile(gc[m & ~tj], [2.5, 97.5]); tl = m & ~tj & ((gc <= lo) | (gc >= hi))
            print('  guard %s outer 5 %% (<= %.1f or >= %.1f; genuine %d): false-junk %.4f' % (gcol, lo, hi, tl.sum(), (J & tl).sum() / max(tl.sum(), 1)))
