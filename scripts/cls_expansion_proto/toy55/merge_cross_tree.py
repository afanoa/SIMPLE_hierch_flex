#!/usr/bin/env python3
"""Cross-tree merge: the leaf states of two INDEPENDENT trees (built on disjoint particle halves A and B) are matched
(Hungarian on their deviations); two matched states p, q are one state unless their difference reproduces ACROSS the
trees: FSC(A_p - A_q, B_p - B_q) >= 0.143 in the band 40 A .. the pooled resolution (best 3 shells), inside the union of
their region masks or the whole envelope. Within-tree halves cannot do this: a board selects particles by their match to
its own (noisy) cell map, and that selected noise is in BOTH of the state's halves (the null set keeps 15 states).
A leaf without a partner in the other tree is not reproducible and is dropped.
usage: merge_cross_tree.py RUN_A RUN_B CONSENSUS OUT ; env LEAFGLOB, FSCMAX, MINREC, MASK (region|sphere), NPROC, FAMS (family glob), SPECMATCH, DOCK, NUIS (1, or 2 = + solvent/background terms fitted over the whole sphere), BANDLO (40 A), THRMIN (threshold floor, e.g. from a shuffled-label control), ARCH (1: prune groups whose map is a reproducible blend of others), STAT (window: best 3-shell FSC mean | z: coefficient-weighted FSC sum over the band), MATCHBAND (low: pair the trees' groups on shells coarser than BANDLO only, so the pairing cannot bias the tested band | lp15: the old 15 A low-pass), NULLCAL (merge threshold = that quantile of the statistic on unmatched half-map noise differences)"""
import sys, os, glob, re, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc; from merge_states_lib import leaf_ok
from scipy.optimize import linear_sum_assignment
ap = argparse.ArgumentParser(); ap.add_argument('runA'); ap.add_argument('runB'); ap.add_argument('cons'); ap.add_argument('out'); ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--crop', type=int, default=80)
ap.add_argument('--margin', type=float, default=1.0); ap.add_argument('--minmatch', type=float, default=0.0, help='matched pairs below this cross-tree deviation correlation (15 A low-pass) are unpaired'); a = ap.parse_args()
cons = read_mrc(a.cons).astype(np.float32); n = cons.shape[0]; x = np.arange(n) - n // 2
rr = np.sqrt(x[:, None, None] ** 2 + x[None, :, None] ** 2 + x[None, None, :] ** 2); msk = np.clip((0.45 * n - rr) / (0.05 * n), 0, 1)
f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); G = np.exp(-0.5 * (kx**2 + ky**2 + kz**2) * (a.smpd * 10. / 1.2) ** 2 / a.smpd ** 2)
lpc = np.real(np.fft.ifftn(np.fft.fftn(cons) * G)); supp = lpc > 0.2 * lpc.max(); c = a.crop; sl = np.r_[0:c // 2, n - c // 2:n]
Gm = np.exp(-3. * (kx**2 + ky**2 + kz**2) / (a.smpd / 15.) ** 2)
def crop_fft(v): return np.fft.fftn(v * msk)[np.ix_(sl, sl, sl)]
def crop_real(v): return np.real(np.fft.ifftn(np.fft.fftn(v)[np.ix_(sl, sl, sl)])).astype(np.float32) * (c / n) ** 3
def dev(fn):
    v = read_mrc(fn).astype(np.float32); A = np.stack([lpc[supp], np.ones(supp.sum())], 1); ab = np.linalg.lstsq(A, v[supp], rcond=None)[0]
    return ((v - ab[1]) / max(ab[0], 1e-12) - lpc).astype(np.float32)
fc = np.fft.fftfreq(c) * c / n; kxc, kyc, kzc = np.meshgrid(fc, fc, fc, indexing='ij'); sh = np.minimum(np.round(np.sqrt(kxc**2 + kyc**2 + kzc**2) * n).astype(int), n // 2).ravel(); nsh = n // 2 + 1
nsh_cnt = np.bincount(sh, minlength=nsh).astype(float); Gl = np.exp(-3. * (kxc**2 + kyc**2 + kzc**2) / (a.smpd / 15.) ** 2)
envc = np.clip(crop_real(np.real(np.fft.ifftn(np.fft.fftn((lpc > 0.1 * lpc.max()).astype(np.float32)) * Gm))), 0, 1); _mc = {}
def region_mask(run, name):
    m = re.match(r'f(\d)_(motion_)?_?(\d+)', name)
    if m is None or os.environ.get('MASK', 'region') != 'region': return envc
    fam, mot, reg = m.group(1), m.group(2), m.group(3); fn = '%s/fam%s/%s/region_%s.mrc' % (run, fam, 'regions_motion' if mot else 'regions', reg)
    if fn not in _mc:
        if os.path.isfile(fn): v = read_mrc(fn).astype(np.float32); v = np.real(np.fft.ifftn(np.fft.fftn(v) * Gm)); _mc[fn] = np.clip(crop_real(v), 0, 1) * envc
        else: print('  region mask %s MISSING: whole envelope' % fn); _mc[fn] = envc
    return _mc[fn]
LCs = None
def specmatch(De, Do):   # SPECMATCH=1: leaf scaled per shell to the consensus signal spectrum (half-map cross spectrum), so a blurred copy is not a difference
    global LCs
    if LCs is None: LCs = crop_fft(lpc); LCs = (LCs, np.sqrt(np.maximum(np.bincount(sh, (np.abs(LCs) ** 2).ravel(), nsh) / np.maximum(nsh_cnt, 1), 0)))
    LC, sref = LCs; E_, O_ = De + LC, Do + LC; si = np.sqrt(np.maximum(np.bincount(sh, np.real(E_ * np.conj(O_)).ravel(), nsh) / np.maximum(nsh_cnt, 1), 0))
    w = np.where(si > 0, np.clip(sref / np.maximum(si, 1e-30), 0.2, 5.), 1.)[sh].reshape(E_.shape); return w * E_ - LC, w * O_ - LC
def load(run):
    names, RE, RO, cnt, MK = [], [], [], [], []
    for d in sorted(glob.glob(run + '/fam' + os.environ.get('FAMS', '?') + '/' + os.environ.get('LEAFGLOB', 'region*') + '/level1/maps')):
        m = re.search(r'fam(\d)/(?:zoom)?region([^/]+)/', d); tag = 'f%s_%s' % (m.group(1), m.group(2)); L1 = os.path.dirname(d); nrec = {}
        if os.path.isfile(L1 + '/clusters.txt'):
            for l in open(L1 + '/clusters.txt'):
                if not l.startswith('#'): nrec[int(l.split()[0])] = int(l.split()[3])
        for fe in sorted(glob.glob(d + '/recvol_state??_even.mrc')):
            s = int(re.findall(r'state(\d+)_even', fe)[0]); fo = fe.replace('_even', '_odd')
            if not os.path.isfile(fo) or not leaf_ok(L1, s): continue
            De, Do = crop_fft(dev(fe)), crop_fft(dev(fo))
            if os.environ.get('SPECMATCH', '0') == '1': De, Do = specmatch(De, Do)
            nm = '%s_s%02d' % (tag, s); names.append(nm); RE.append(np.real(np.fft.ifftn(De)).astype(np.float32)); RO.append(np.real(np.fft.ifftn(Do)).astype(np.float32)); cnt.append(nrec.get(s, 1)); MK.append(region_mask(run, nm))
    return names, np.array(RE), np.array(RO), np.array(cnt, float), MK
def _dock(args):   # rigid-body fit of one leaf's 15 A map onto the consensus (6 dof, Powell), applied to both halves
    from scipy.ndimage import affine_transform; from scipy.optimize import minimize; from scipy.spatial.transform import Rotation
    E, O, lc = args; ref = np.real(np.fft.ifftn(np.fft.fftn(lc) * Gl)); ctr = (np.array(ref.shape) - 1) / 2.; env = ref > 0.15 * ref.max(); r_ = ref[env] - ref[env].mean()
    def xf(v, p): Rm = Rotation.from_rotvec(np.radians(p[:3])).as_matrix(); return affine_transform(v, Rm, offset=ctr - Rm @ ctr - p[3:], order=1, mode='constant')
    v = np.real(np.fft.ifftn(np.fft.fftn(0.5 * (E + O) + lc) * Gl))
    def cost(p): w = xf(v, p)[env]; w = w - w.mean(); return -np.dot(w, r_) / max(np.linalg.norm(w) * np.linalg.norm(r_), 1e-30)
    res = minimize(cost, np.zeros(6), method='Powell', options={'xtol': 0.05, 'ftol': 1e-6, 'maxfev': 800}); p = res.x
    return (xf(E + lc, p) - lc).astype(np.float32), (xf(O + lc, p) - lc).astype(np.float32), np.linalg.norm(p[:3]), np.linalg.norm(p[3:]) * a.smpd * n / c, -cost(np.zeros(6)), -res.fun
def dock_all(names, RE, RO):
    from multiprocessing import Pool; lc = np.real(np.fft.ifftn(crop_fft(lpc))).astype(np.float32)
    with Pool(int(os.environ.get('NPROC', '16'))) as pool: out = pool.map(_dock, [(RE[i], RO[i], lc) for i in range(len(names))])
    for i, o in enumerate(out): RE[i], RO[i] = o[0], o[1]; print('  dock %s: rot %.2f deg, shift %.2f A, corr %.4f -> %.4f' % (names[i], o[2], o[3], o[4], o[5]), flush=True)
    return RE, RO
nA, REA, ROA, cA, MA = load(a.runA); nB, REB, ROB, cB, MB = load(a.runB); print('tree A: %d leaf states, tree B: %d leaf states%s' % (len(nA), len(nB), ' (spectral match to the consensus)' if os.environ.get('SPECMATCH', '0') == '1' else ''), flush=True)
if os.environ.get('DOCK', '0') == '1': REA, ROA = dock_all(nA, REA, ROA); REB, ROB = dock_all(nB, REB, ROB)
def fsc_of(E, O):
    num = np.bincount(sh, np.real(E * np.conj(O)).ravel(), nsh); d1 = np.bincount(sh, (np.abs(E) ** 2).ravel(), nsh); d2 = np.bincount(sh, (np.abs(O) ** 2).ravel(), nsh); return num / np.maximum(np.sqrt(d1 * d2), 1e-30)
def res_of(E, O):
    fsc = fsc_of(E, O); below = np.where(fsc[1:n // 2] < 0.143)[0]; k = (below[0] + 1) if len(below) else n // 2 - 1; return n * a.smpd / k
_klo_m = max(int(round(n * a.smpd / float(os.environ.get('BANDLO', '40')))), 2); _mband = ((sh >= 2) & (sh < _klo_m - 1)).reshape(c, c, c)
def corr_lp(u, v, m):   # pairing similarity: Fourier correlation on the shells coarser than the tested band, one shell of gap (default), or the old masked 15 A real-space correlation (MATCHBAND=lp15)
    if os.environ.get('MATCHBAND', 'low') != 'lp15': U, V = np.fft.fftn(u)[_mband], np.fft.fftn(v)[_mband]; return float(np.real(np.vdot(U, V)) / max(np.linalg.norm(U) * np.linalg.norm(V), 1e-30))
    ul = np.real(np.fft.ifftn(np.fft.fftn(u) * Gl)) * m; vl = np.real(np.fft.ifftn(np.fft.fftn(v) * Gl)) * m; ul -= ul.mean(); vl -= vl.mean(); return float(np.dot(ul.ravel(), vl.ravel()) / max(np.linalg.norm(ul) * np.linalg.norm(vl), 1e-30))
# 1. match the trees' leaves (Hungarian on the cross-tree deviation correlation, union of the two region masks)
S = np.array([[corr_lp(0.5 * (REA[i] + ROA[i]), 0.5 * (REB[j] + ROB[j]), np.maximum(MA[i], MB[j])) for j in range(len(nB))] for i in range(len(nA))])
ia, ib = linear_sum_assignment(-S); pairs = [(i, j) for i, j in zip(ia, ib) if S[i, j] >= a.minmatch]
print('matched %d pairs (cross-tree deviation correlation median %.2f, min %.2f); unpaired: A %s, B %s' % (len(pairs), np.median([S[i, j] for i, j in pairs]) if pairs else 0, min([S[i, j] for i, j in pairs]) if pairs else 0,
      ' '.join(nA[i] for i in range(len(nA)) if i not in [p[0] for p in pairs]) or '-', ' '.join(nB[j] for j in range(len(nB)) if j not in [p[1] for p in pairs]) or '-'), flush=True)
for i, j in pairs: print('  %s <-> %s  corr %.2f' % (nA[i], nB[j], S[i, j]))
# 2. pooled groups over matched pairs; the merge statistic is the cross-tree FSC of the difference
P = len(pairs); groups = {p: [p] for p in range(P)}
GA = {p: 0.5 * (REA[i] + ROA[i]) for p, (i, j) in enumerate(pairs)}; GB = {p: 0.5 * (REB[j] + ROB[j]) for p, (i, j) in enumerate(pairs)}
GAE, GAO = {p: REA[i].copy() for p, (i, j) in enumerate(pairs)}, {p: ROA[i].copy() for p, (i, j) in enumerate(pairs)}; GN = {p: cA[i] + cB[j] for p, (i, j) in enumerate(pairs)}; GM = {p: np.maximum(MA[i], MB[j]) for p, (i, j) in enumerate(pairs)}
def band_stat(DA, DB, PA, PB):
    d = max(res_of(PA, PB) * a.margin, 2.5 * a.smpd); fsc = fsc_of(DA, DB); klo, khi = max(int(round(n * a.smpd / float(os.environ.get('BANDLO', '40')))), 2), int(round(n * a.smpd / d))
    if os.environ.get('STAT', 'window') == 'z':   # every shell of the band, weighted by its independent coefficients: ~N(0,1) for unrelated differences
        w = 0.5 * nsh_cnt[klo:khi + 1]; return float(np.sum(w * fsc[klo:khi + 1]) / np.sqrt(max(np.sum(w), 1.))), d
    return max(float(fsc[k:k + 3].mean()) for k in range(klo, max(khi - 2, klo) + 1)), d
_NB = {}
def nuis_resid(D, Pdev):   # NUIS=1: the part of a difference that no pose/shift/blur/scale/spectral-shape nuisance of the pooled map explains (34-map basis, fit over the envelope)
    if not _NB:
        k = np.fft.fftfreq(c); _NB['k'] = np.meshgrid(k, k, k, indexing='ij'); g = np.arange(c) - c // 2; _NB['x'] = np.meshgrid(g, g, g, indexing='ij')
        _NB['lp'] = np.exp(-3. * sum(kk ** 2 for kk in _NB['k']) / (a.smpd * n / c / 8.) ** 2); _NB['env'] = envc.ravel() > 0.5
        _NB['lc'] = np.real(np.fft.ifftn(crop_fft(lpc))).astype(np.float32); edges = (3, 6, 9, 13, 18, 25, 41); _NB['bands'] = [((sh >= lo) & (sh < hi)).reshape(c, c, c) for lo, hi in zip(edges[:-1], edges[1:])]
        mc = crop_real(msk.astype(np.float32)); r2 = sum(xx.astype(np.float32) ** 2 for xx in _NB['x']) / (c / 2.) ** 2; _NB['bg'] = [mc, mc * r2, mc * r2 ** 2, envc]   # solvent level, radial background, particle/solvent step
        if os.environ.get('NUIS', '0') == '2': _NB['env'] = mc.ravel() > 0.5
    K, X = _NB['k'], _NB['x']; F = np.fft.fftn(Pdev + _NB['lc']) * _NB['lp']
    d = lambda G, i: G * (2j * np.pi * K[i]); real = lambda G: np.real(np.fft.ifftn(G))
    def L(G, i):   # rotation generator about the box centre, on a Fourier map
        j, l = (i + 1) % 3, (i + 2) % 3; return np.fft.fftn(X[j] * real(d(G, l)) - X[l] * real(d(G, j)))
    D1 = [d(F, i) for i in range(3)]; L1 = [L(F, i) for i in range(3)]
    basis = [F] + D1 + L1 + [d(D1[i], j) for i in range(3) for j in range(i, 3)] + [0.5 * (L(L1[j], i) + L(L1[i], j)) for i in range(3) for j in range(i, 3)] + [L(D1[j], i) for i in range(3) for j in range(3)] + [F * b for b in _NB['bands']]
    Bm = np.stack([real(b).ravel() for b in basis] + ([b.ravel() for b in _NB['bg']] if os.environ.get('NUIS', '0') == '2' else []), 1); e = _NB['env']; beta = np.linalg.lstsq(Bm[e], D.ravel()[e], rcond=None)[0]
    return (D.ravel() - Bm @ beta).reshape(D.shape).astype(np.float32)
def dres(p, q):
    out = []; wp, wq = GN[p], GN[q]; dA, dB = GA[p] - GA[q], GB[p] - GB[q]
    if os.environ.get('NUIS', '0') in ('1', '2'): dA, dB = nuis_resid(dA, (wp * GA[p] + wq * GA[q]) / (wp + wq)), nuis_resid(dB, (wp * GB[p] + wq * GB[q]) / (wp + wq))
    for m in (np.maximum(GM[p], GM[q]), envc):
        DA, DB = np.fft.fftn(dA * m), np.fft.fftn(dB * m)
        PA, PB = np.fft.fftn(((wp * GA[p] + wq * GA[q]) / (wp + wq)) * m), np.fft.fftn(((wp * GB[p] + wq * GB[q]) / (wp + wq)) * m); out.append(band_stat(DA, DB, PA, PB)[0])
    return max(out)
def rep(p):   # the pooled state's own cross-tree reproducibility in its band (support across trees)
    m = GM[p]; return band_stat(np.fft.fftn(GA[p] * m), np.fft.fftn(GB[p] * m), np.fft.fftn(GA[p] * m), np.fft.fftn(GB[p] * m))[0]
def _dres_pair(pq): return pq, dres(*pq)
from multiprocessing import Pool; nproc = int(os.environ.get('NPROC', '16'))
with Pool(nproc) as pool: R = dict(pool.map(_dres_pair, [(p, q) for p in range(P) for q in range(p + 1, P)]))
print('pairs at the start: cross-tree difference FSC median %.2f, %d of %d pairs unresolved (< 0.143); state support (own cross-tree band FSC) median %.2f min %.2f' % (np.median(list(R.values())) if R else 0, sum(v < 0.143 for v in R.values()), len(R), np.median([rep(p) for p in range(P)]), min(rep(p) for p in range(P))), flush=True)
THR = 0.143
if os.environ.get('NULLCAL', ''):   # null of the band statistic: two independent Gaussian fields under the pair's masks and windows (half maps are not exchangeable: cross-half EM ties the half to the assignment)
    def _null(pr):
        p, r, seed = pr; g = np.random.default_rng(seed); dA, dB = g.standard_normal((c, c, c)).astype(np.float32), g.standard_normal((c, c, c)).astype(np.float32)
        if os.environ.get('NUIS', '0') in ('1', '2'): dA, dB = nuis_resid(dA, GA[p]), nuis_resid(dB, GB[r])
        return max(band_stat(np.fft.fftn(dA * m), np.fft.fftn(dB * m), np.fft.fftn(GA[p] * m), np.fft.fftn(GB[r] * m))[0] for m in (np.maximum(GM[p], GM[r]), envc))
    with Pool(nproc) as pool: NUL = np.array(pool.map(_null, [(p, r, 1000 * p + r) for p in range(P) for r in range(P) if r != p]))
    THR = max(0.143, float(np.quantile(NUL, float(os.environ['NULLCAL'])))); print('null of the statistic (%d independent-noise pairs under the pair masks): median %.3f, 95 %% %.3f, 99 %% %.3f, max %.3f -> merge threshold %.3f' % (len(NUL), np.median(NUL), np.quantile(NUL, 0.95), np.quantile(NUL, 0.99), NUL.max(), THR), flush=True)
THR = max(THR, float(os.environ.get('THRMIN', '0')))
if os.environ.get('DUMPPAIR', ''):   # the most and least extreme pairs' differences (after projection) for inspection
    for tag, (p, q) in (('max', max(R, key=R.get)), ('min', min(R, key=R.get))):
        wp, wq = GN[p], GN[q]; dA, dB = GA[p] - GA[q], GB[p] - GB[q]
        if os.environ.get('NUIS', '0') in ('1', '2'): dA, dB = nuis_resid(dA, (wp * GA[p] + wq * GA[q]) / (wp + wq)), nuis_resid(dB, (wp * GB[p] + wq * GB[q]) / (wp + wq))
        np.savez('%s_%s.npz' % (os.environ['DUMPPAIR'], tag), dA=dA, dB=dB, rawA=GA[p] - GA[q], rawB=GB[p] - GB[q], PA=GA[p], fsc=fsc_of(np.fft.fftn(dA * envc), np.fft.fftn(dB * envc)), fscraw=fsc_of(np.fft.fftn((GA[p] - GA[q]) * envc), np.fft.fftn((GB[p] - GB[q]) * envc)), nsh_cnt=nsh_cnt, stat=R[(p, q)], env=envc)
        print('dumped %s pair %s statistic %.2f' % (tag, (p, q), R[(p, q)]), flush=True)
print('pair statistics at the start (sorted): %s; threshold %.3f' % (' '.join('%.2f' % v for v in sorted(R.values())), THR), flush=True)
merged_at = []
while R and min(R.values()) < THR:
    (p, q), v = min(R.items(), key=lambda kv: kv[1]); merged_at.append(v); wp, wq = GN[p], GN[q]
    for D in (GA, GB, GAE, GAO): D[p] = (wp * D[p] + wq * D[q]) / (wp + wq); del D[q]
    GN[p] = wp + wq; GM[p] = np.maximum(GM[p], GM[q]); groups[p] += groups.pop(q); del GN[q], GM[q]; R = {k: v for k, v in R.items() if q not in k}
    with Pool(nproc) as pool: R.update(dict(pool.map(_dres_pair, [tuple(sorted((p, o))) for o in groups if o != p])))
def _mix_resid(G, p, others):   # residual of group p's map after the best convex blend of the other groups' maps (weights >= 0, sum 1), fit over the envelope
    from scipy.optimize import nnls
    e = envc.ravel() > 0.5; X = np.stack([G[o].ravel()[e] for o in others], 1); y = G[p].ravel()[e]; lam = 10. * max(np.abs(X).max(), 1e-12) * np.sqrt(e.sum())
    w = nnls(np.vstack([X, lam * np.ones((1, len(others)))]), np.r_[y, lam])[0]; return G[p] - sum(wi * G[o] for wi, o in zip(w, others)), w
def _arch_stat(p):
    others = [o for o in groups if o != p]; rA, wA = _mix_resid(GA, p, others); rB, wB = _mix_resid(GB, p, others)
    if os.environ.get('NUIS', '0') in ('1', '2'): rA, rB = nuis_resid(rA, GA[p]), nuis_resid(rB, GB[p])
    v = max(band_stat(np.fft.fftn(rA * m), np.fft.fftn(rB * m), np.fft.fftn(GA[p] * m), np.fft.fftn(GB[p] * m))[0] for m in (GM[p], envc))
    return p, v, others[int(np.argmax(wA + wB))], wA + wB
if os.environ.get('ARCH', '0') == '1':   # a group whose map is reproducibly explained as a blend of the others is a mixture, not a state: its members join its main component
    while len(groups) > 2:
        with Pool(nproc) as pool: AS = pool.map(_arch_stat, list(groups))
        p, v, tgt, w = min(AS, key=lambda t: t[1]); print('  archetype test: least-novel group %d residual FSC %.2f (threshold %.3f), blend weights %s' % (p, v, THR, ' '.join('%.2f' % (x / 2) for x in w)), flush=True)
        if v >= THR: break
        groups[tgt] += groups.pop(p); GN[tgt] += GN.pop(p); GM.pop(p)
        for D in (GA, GB, GAE, GAO): D.pop(p)
        R = {k: val for k, val in R.items() if p not in k}
    print('archetypes: %d' % len(groups))
K = len(groups); lab = np.zeros(P, int)
for g, (p, mem) in enumerate(sorted(groups.items(), key=lambda kv: -GN[kv[0]])): lab[mem] = g + 1
print('%d merges -> %d states; merges at cross-tree difference FSC %s' % (len(merged_at), K, ' '.join('%.2f' % v for v in merged_at)))
print('remaining pairs: cross-tree difference FSC median %.2f, least-resolved pair %.2f (threshold %.3f)' % (np.median(list(R.values())), min(R.values()), THR) if R else 'one state')
labA = np.zeros(len(nA), int); labB = np.zeros(len(nB), int)
for p, (i, j) in enumerate(pairs): labA[i] = lab[p]; labB[j] = lab[p]
np.savez(a.out + '_A.npz', names=np.array(nA), labels=labA, floor=12.); np.savez(a.out + '_B.npz', names=np.array(nB), labels=labB, floor=12.)
with open(a.out + '.txt', 'w') as o:
    o.write('# cross-tree merge: %d + %d leaves, %d matched pairs -> %d states (label 0 = unpaired, dropped)\n' % (len(nA), len(nB), P, K))
    for p, (i, j) in enumerate(pairs): o.write('%d %s %s %.2f\n' % (lab[p], nA[i], nB[j], S[i, j]))
print('-> %s.txt (+_A.npz, _B.npz for consolidate_states.py on each half)' % a.out)
