#!/usr/bin/env python3
"""hetero_tree audit items: (1) mask core by a relative share rule + data-derived falloff, (2) SVD site factors with a share
threshold, (3) a diffuse-variance stopping rule. Python port of half_variance / falloff_mask / discover_sites (validated by the
2D sub-class session, Dice 0.98-1.00) applied to a node's saved half maps: the node-mask maps (flex at a project root, cheap
children at a subtree root, the parent's recB0 state maps below the root), the recB0 state maps and the stage-B state maps.
usage: vmask_proto.py --node DIR --tag T --out DIR --mskdiam A [--parent PDIR] [--labels-cheap F] [--labels-stageb F]
       [--cons C.mrc --closed X.mrc | --truth a.mrc,b.mrc,...] [--no-disc] [--sources mask,recB0,stageB]"""
import argparse, os, re, sys, time, math
import numpy as np
from scipy import ndimage
MLP, MSIG, MFRAC, MDIL, MFALL, MSUP = 15., 1.5, 0.2, 3., 6., 0.1
DMINVOX, DSEEDMIN, DMINSHARE, DFACT_ITS, DFACT_COS, DFACT_SHARE, DDICE, DMAXMOL = 200, 1000, 0.05, 20, 0.8, 0.1, 0.9, 0.3
SHARES = (0.5, 0.6, 0.7, 0.8); FSHARES = (0.1, 0.15, 0.2, 0.25); ZFAC = 3.
ap = argparse.ArgumentParser(); ap.add_argument('--node', required=True); ap.add_argument('--tag', required=True); ap.add_argument('--out', required=True)
ap.add_argument('--mskdiam', type=float, required=True); ap.add_argument('--parent', default=''); ap.add_argument('--labels-cheap', default='')
ap.add_argument('--labels-stageb', default=''); ap.add_argument('--cons', default=''); ap.add_argument('--closed', default=''); ap.add_argument('--truth', default='')
ap.add_argument('--no-disc', action='store_true'); ap.add_argument('--sources', default='mask,recB0,stageB')
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); T0 = time.time(); REP = open(os.path.join(a.out, 'report.txt'), 'w')
def say(*s):
    msg = ' '.join(str(x) for x in s); print('[%6.0f s] %s' % (time.time() - T0, msg), flush=True); REP.write(msg + '\n'); REP.flush()
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); assert h[3] == 2, 'mode %d' % h[3]
    return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n)
def wr(p, x, smpd):
    x = np.ascontiguousarray(x, np.float32); h = np.zeros(256, np.int32); h[0:3] = x.shape; h[3] = 2; h[7:10] = x.shape; hf = h.view(np.float32); hf[10:13] = np.array(x.shape) * smpd; hf[13:16] = 90.
    h[16:19] = [1, 2, 3]; hf[19], hf[20], hf[21] = x.min(), x.max(), x.mean(); h[52] = int.from_bytes(b'MAP ', 'little'); h[53] = 0x4444
    with open(p, 'wb') as f: f.write(h.tobytes()); f.write(x.tobytes())
def hdr(p):
    h = np.fromfile(p, np.int32, 256); return int(h[0]), float(np.fromfile(p, np.float32, 13)[10]) / int(h[0])
# ---------------------------------------------------------------- map sets
def labels_of(path, col=None):
    if col is None: return np.loadtxt(path).astype(int)
    return np.loadtxt(path, comments='#', usecols=col).astype(int)
def map_set(d, sub, lab):
    st = [s for s in range(1, lab.max() + 1) if (lab == s).sum() > 0 and os.path.exists(os.path.join(d, sub, 'recvol_state%02d_even_unfil.mrc' % s))]
    E = [os.path.join(d, sub, 'recvol_state%02d_even_unfil.mrc' % s) for s in st]; O = [e.replace('_even_', '_odd_') for e in E]
    return E, O, np.array([(lab == s).sum() for s in st], float), np.array(st)
SETS = {}
src = a.sources.split(',')
if 'mask' in src:
    if a.parent:
        plab = labels_of(os.path.join(a.parent, 'sites', 'disc_labels.txt')); prows = np.loadtxt(os.path.join(a.parent, 'rows.txt')).astype(int)
        crows = np.loadtxt(os.path.join(a.node, 'rows.txt')).astype(int); pos = {r: i for i, r in enumerate(prows)}
        lab = np.array([plab[pos[r]] for r in crows]); SETS['mask'] = map_set(a.parent, 'recB0', lab) + ('parent recB0',)
    elif os.path.exists(os.path.join(a.node, 'cheap')) and a.labels_cheap:
        SETS['mask'] = map_set(a.node, 'cheap', labels_of(a.labels_cheap)) + ('cheap children',)
    else:
        SETS['mask'] = map_set(a.node, 'flex', labels_of(os.path.join(a.node, 'flex', 'flex_pca_coordinates.txt'), 2)) + ('flex states',)
if 'recB0' in src and os.path.exists(os.path.join(a.node, 'sites', 'disc_labels.txt')):
    SETS['recB0'] = map_set(a.node, 'recB0', labels_of(os.path.join(a.node, 'sites', 'disc_labels.txt'))) + ('recB0 states',)
if 'stageB' in src and os.path.exists(os.path.join(a.node, 'stageB')):
    if a.labels_stageb and os.path.exists(a.labels_stageb): lab = labels_of(a.labels_stageb)
    else:
        ks = sorted(int(f[12:14]) for f in os.listdir(os.path.join(a.node, 'stageB')) if re.match(r'recvol_state\d\d_even_unfil\.mrc', f)); lab = np.array(ks)
    SETS['stageB'] = map_set(a.node, 'stageB', lab) + ('stage-B states',)
GEOM = {}
def geometry(nn, sp):
    """per box: the 15 A low-pass and 1.5 px smoothing filters, the radius grid and the sphere (set as globals for the active map set)"""
    global n, smpd, LPF, SMF, R, sph
    if (nn, sp) not in GEOM:
        f = np.fft.fftfreq(nn); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(nn), indexing='ij'); f2 = kx ** 2 + ky ** 2 + kz ** 2; del kx, ky, kz
        g = np.arange(nn) - nn // 2; Rr = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2)
        GEOM[(nn, sp)] = (np.exp(-f2 / sp ** 2 * MLP ** 2 * np.log(2.)).astype(np.float32), np.exp(-2. * np.pi ** 2 * MSIG ** 2 * f2).astype(np.float32), Rr, Rr <= a.mskdiam / 2. / sp)
    n, smpd = nn, sp; LPF, SMF, R, sph = GEOM[(nn, sp)]
lp = lambda x: np.fft.irfftn(np.fft.rfftn(x) * LPF, s=(n, n, n)).astype(np.float32); smooth = lambda x: np.fft.irfftn(np.fft.rfftn(x) * SMF, s=(n, n, n))
say('%s: sets: %s' % (a.tag, ', '.join('%s=%d %s (box %d, %.3f A/px)' % (k, len(v[0]), v[4], *hdr(v[0][0])) for k, v in SETS.items())))
def load_lp(files):
    X = np.empty((len(files), n, n, n), np.float32)
    for j, p in enumerate(files): X[j] = lp(rd(p))
    return X
def half_variance(E, O, w):
    wn = (w / w.sum()).astype(np.float64); em = np.tensordot(wn, E, 1); om = np.tensordot(wn, O, 1); mean = 0.5 * (em + om)
    sup = (mean > MSUP * mean[sph].max()) & sph
    def sc(x, ref): return float(x[sup] @ ref[sup]) / max(float(x[sup] @ x[sup]), 1e-30)
    for j in range(len(E)): E[j] *= sc(E[j], em); O[j] *= sc(O[j], om)
    em = np.tensordot(wn, E, 1); om = np.tensordot(wn, O, 1); v = np.zeros((n, n, n))
    for j in range(len(E)): E[j] -= em; O[j] -= om; v += wn[j] * E[j].astype(np.float64) * O[j]
    vs = smooth(v).astype(np.float32); vs[~sph] = 0.; return vs, sup, wn
# ---------------------------------------------------------------- cores, falloff, correlation length
def share_core(vs, sup, S):
    """the level set of vs inside the support holding share S of the positive variance there: core and its cut"""
    idx = np.flatnonzero(sup); vals = np.clip(vs.ravel()[idx].astype(np.float64), 0., None); o = np.argsort(-vals, kind='stable'); cs = np.cumsum(vals[o])
    k = min(int(np.searchsorted(cs, S * cs[-1])) + 1, len(o)); core = np.zeros(vs.shape, bool); core.ravel()[idx[o[:k]]] = True
    return core, float(vals[o[k - 1]])
def top_share(vs, sup, fr):
    vals = np.clip(vs[sup].astype(np.float64), 0., None); o = np.sort(vals)[::-1]; k = max(1, int(fr * len(o))); return float(o[:k].sum() / max(o.sum(), 1e-30))
def falloff(core, mfall):
    if not core.any(): return np.where(sph, 1., 0.).astype(np.float32)
    core2 = ndimage.distance_transform_edt(~core) <= MDIL; m = np.exp(-ndimage.distance_transform_edt(~core2) ** 2 / (2. * mfall ** 2)).astype(np.float32); m[~sph] = 0.; return m
def acf_length(vs, sup):
    """mask-corrected radial autocorrelation of vs inside the support: the lags where it falls to 1/2 and 1/e (px)"""
    x = np.where(sup, vs - vs[sup].mean(), 0.).astype(np.float64); m = sup.astype(np.float64)
    ax = np.fft.fftshift(np.fft.irfftn(np.abs(np.fft.rfftn(x)) ** 2, s=x.shape)); am = np.fft.fftshift(np.fft.irfftn(np.abs(np.fft.rfftn(m)) ** 2, s=x.shape))
    ok = am > 0.25 * am.max(); rb = np.round(R).astype(int); prof = []
    for i in range(0, 40):
        q = ok & (rb == i)
        if q.sum() == 0: break
        prof.append(float((ax[q] / am[q]).mean()))
    prof = np.array(prof) / prof[0]
    def lag(level):
        for i in range(1, len(prof)):
            if prof[i] < level: return i - 1 + (prof[i - 1] - level) / max(prof[i - 1] - prof[i], 1e-12)
        return float(len(prof))
    return lag(0.5), lag(math.exp(-1.)), prof
# ---------------------------------------------------------------- regions
REG, RNAME, NREG = [], [], -1
if a.cons and a.closed:
    NREG, spR = hdr(a.cons); fz = np.fft.fftfreq(NREG, d=spR); KX, KY, KZ = np.meshgrid(fz, fz, fz, indexing='ij'); G = np.exp(-0.5 * (KX ** 2 + KY ** 2 + KZ ** 2) * (10. / 1.2) ** 2); del KX, KY, KZ
    l10 = lambda v: np.real(np.fft.ifftn(np.fft.fftn(v.astype(np.float32)) * G)); cons = l10(rd(a.cons)); supp = cons > 0.2 * cons.max(); outside = cons < 0.02 * cons.max()
    std = lambda v: (v - v[outside].mean()) / max(np.sqrt(np.mean((v - v[outside].mean())[supp] ** 2)), 1e-12); cons = std(cons); closed = std(l10(rd(a.closed))); d = closed - cons
    thr = 0.5 * np.percentile(np.abs(d[supp | (closed > 0.2 * closed.max())]), 99); REG = [d > thr, d < -thr]; RNAME = ['closed-arm', 'open-arm']; del G
if a.truth:
    REG = [rd(p) > 0.5 for p in a.truth.split(',')]; RNAME = [os.path.basename(p).replace('.mrc', '').replace('_travel', '') for p in a.truth.split(',')]; NREG = REG[0].shape[0]
ANY = np.any(REG, 0) if REG else None
def regs_ok(): return bool(REG) and n == NREG
def cover(m, reg): return float(((m > 0.5) & reg).sum()) / max(1, reg.sum())
def reg_line(m):
    if not regs_ok(): return ''
    return 'cover [%s] prec %.2f' % (' '.join('%.2f' % cover(m, r) for r in REG), float(((m > 0.5) & ANY).sum()) / max(1, (m > 0.5).sum()))
def pct(m, sup): return 100. * float(((m > 0.5) & sup).sum()) / max(1, sup.sum())
# ---------------------------------------------------------------- (1) + (3) per map set
def top_dirs(p, m):
    lam, u = np.linalg.eigh(p.T.astype(np.float64) @ p.astype(np.float64)); return u[:, np.argsort(-lam)[:m]].T
KEEP = {}
for key, (Ef, Of, w, names, what) in SETS.items():
    geometry(*hdr(Ef[0])); E = load_lp(Ef); O = load_lp(Of); vs, sup, wn = half_variance(E, O, w); nsup = int(sup.sum())
    if regs_ok(): say('%s %s: regions [%s] hold %s voxels = %s %% of the support' % (a.tag, key, ', '.join(RNAME), ', '.join(str(int(r.sum())) for r in REG), ', '.join('%.1f' % (100. * (r & sup).sum() / nsup) for r in REG)))
    core_old = vs > MFRAC * vs.max(); m_old = falloff(core_old, MFALL)
    if key == 'mask' and os.path.exists(os.path.join(a.node, 'mask.mrc')):
        mf = rd(os.path.join(a.node, 'mask.mrc')); say('%s %s: port vs tree mask.mrc Dice %.3f (tree %.1f %% of the molecule)' % (a.tag, key, 2. * ((m_old > 0.5) & (mf > 0.5)).sum() / max(1, (m_old > 0.5).sum() + (mf > 0.5).sum()), pct(mf, sup)))
    h2, he, prof = acf_length(vs, sup)
    say('%s %s: %d maps (%s), peak %.3e, variance share of the top 1/2/5/10/20 %% of support voxels %s; acf half-lag %.1f px, 1/e lag %.1f px (%.1f / %.1f A)'
        % (a.tag, key, len(E), what, vs.max(), ' '.join('%.2f' % top_share(vs, sup, fr) for fr in (0.01, 0.02, 0.05, 0.1, 0.2)), h2, he, h2 * smpd, he * smpd))
    say('%s %s: OLD core 0.2 x peak: core %.1f %% of molecule, mask > 0.5 %.1f %% %s' % (a.tag, key, pct(core_old, sup), pct(m_old, sup), reg_line(m_old)))
    for S in SHARES:
        core, cut = share_core(vs, sup, S); lab, nl = ndimage.label(core); sz = np.sort(np.bincount(lab.ravel())[1:])[::-1]
        m6 = falloff(core, MFALL); ma = falloff(core, max(h2, 1.))
        say('%s %s: share %.1f core: cut %.3f x peak, core %.1f %% of molecule, %d pieces >= %d (largest %s); mask(fall 6) %.1f %% %s; mask(fall acf %.1f) %.1f %% %s'
            % (a.tag, key, S, cut / max(vs.max(), 1e-30), pct(core, sup), int((sz >= DMINVOX).sum()), DMINVOX, sz[:4].tolist(), pct(m6, sup), reg_line(m6), h2, pct(ma, sup), reg_line(ma)))
    KEEP[key] = dict(E=E, O=O, vs=vs, sup=sup, wn=wn, w=w, names=names, n=n)
    if key != 'recB0': del E, O; KEEP[key].pop('E'); KEEP[key].pop('O')
# ---------------------------------------------------------------- (2) discovery variants on the recB0 maps
def discover(vs, sup, dE, dO, wn, counts, names, mode, fshare, core_mode, S, tag):
    """port of discover_sites (c1); mode 'kmeans' = the 2-means rule, 'svd' = weighted SVD factors kept at Fisher z >= ZFAC and share >= fshare"""
    nk = len(dE); nsup = int(sup.sum()); okk = counts >= DSEEDMIN
    if core_mode == 'old': core = vs > MFRAC * vs.max(); cut0 = MFRAC * vs.max()
    else: core, cut0 = share_core(vs, sup, S)
    lab, nl = ndimage.label(core); tot = float(vs[core].sum())
    regs = [dict(v=np.flatnonzero(lab.ravel() == q), cut=cut0, top=True, comp=q) for q in range(1, nl + 1)]; cands = []; np2 = nl; log = []
    neff_per_vox = 2.0944 * (smpd / MLP) ** 3
    def add(c, comp, fac, prof=None):
        idx = np.argwhere(c); lo = np.maximum(idx.min(0) - (int(MDIL) + 8 * int(MFALL)), 0); hi = np.minimum(idx.max(0) + (int(MDIL) + 8 * int(MFALL)) + 1, n); sl = tuple(slice(lo[i], hi[i]) for i in range(3))
        m = np.zeros((n, n, n), np.float32); cb = c[sl]; m[sl] = np.exp(-(ndimage.distance_transform_edt(~ndimage.binary_dilation(cb, iterations=int(MDIL))) * 1.) ** 2 / (2. * MFALL ** 2)); m[~sph] = 0.
        mb = m[sl]; pe = dE[(slice(None),) + sl] * mb; po = dO[(slice(None),) + sl] * mb
        ek = np.einsum('kxyz,kxyz->k', pe, po, dtype=np.float64); e2 = np.einsum('kxyz,kxyz->k', pe, pe, dtype=np.float64); o2 = np.einsum('kxyz,kxyz->k', po, po, dtype=np.float64); del pe, po
        num, dsq, osq = float(wn @ ek), float(wn @ e2), float(wn @ o2); on = m > 0.5
        d = dict(m=m, on=on, r=num / max(np.sqrt(dsq * osq), 1e-300), comp=comp, fac=fac, ncore=int(c.sum()), vox=int(on.sum()), mol=float((on & sup).sum()) / max(1, nsup), share=float(vs[c].sum()) / tot)
        key_ = np.abs(prof) if prof is not None else ek; d['seed'] = int(names[int(np.argmax(np.where(okk, key_, -np.inf)))]); d['tag'] = str(comp) + ('f%d' % fac if fac > 0 else ''); cands.append(d)
    ir = 0
    while ir < len(regs):
        reg = regs[ir]; ir += 1
        if len(reg['v']) < DMINVOX: continue
        cm = np.zeros(n ** 3, bool); cm[reg['v']] = True; cm = cm.reshape(n, n, n)
        if reg['top'] and vs[cm].sum() / tot < DMINSHARE: continue
        ic = reg['comp']; nc0 = len(cands); l_fac = False
        if cm.sum() >= 2 * DMINVOX:
            idx = np.nonzero(cm); nv = len(idx[0])
            if mode == 'kmeans':
                pn = (np.sqrt(wn)[:, None] * 0.5 * (dE[:, idx[0], idx[1], idx[2]].astype(np.float64) + dO[:, idx[0], idx[1], idx[2]])).T
                pn /= np.maximum(np.linalg.norm(pn, axis=1, keepdims=True), 1e-30); vv = vs[idx]; cc = top_dirs(pn, 2)
                for it in range(DFACT_ITS):
                    asg = np.where(np.abs(pn @ cc[0]) >= np.abs(pn @ cc[1]), 1, 2)
                    for j in range(2):
                        if (asg == j + 1).sum() > 2: cc[j] = top_dirs(pn[asg == j + 1], 1)[0]
                cosv = abs(float(cc[0] @ cc[1])); sh = [float(vv[asg == j].sum() / vv.sum()) for j in (1, 2)]
                if cosv < DFACT_COS and min(sh) >= DFACT_SHARE and min((asg == 1).sum(), (asg == 2).sum()) >= DMINVOX:
                    log.append('comp %d: 2 factors (|cos| %.2f, shares %.2f/%.2f)' % (ic, cosv, sh[0], sh[1])); profs = [cc[0], cc[1]]; kept = [1, 2]
                    for j in (1, 2):
                        b = np.zeros((n, n, n), bool); b[tuple(ix[asg == j] for ix in idx)] = True; lb, nb = ndimage.label(b); ps = np.bincount(lb.ravel(), minlength=nb + 1)
                        b = np.isin(lb, 1 + np.flatnonzero(ps[1:] >= max(20, DMINVOX // 4))); add(b, ic, j, prof=cc[j - 1] / np.sqrt(wn))
                    l_fac = True
                else: log.append('comp %d: 1 factor (|cos| %.2f, shares %.2f/%.2f)' % (ic, cosv, sh[0], sh[1]))
            else:
                De = (np.sqrt(wn)[:, None] * dE[:, idx[0], idx[1], idx[2]].astype(np.float64)).T; Do = (np.sqrt(wn)[:, None] * dO[:, idx[0], idx[1], idx[2]].astype(np.float64)).T
                lam, V = np.linalg.eigh(De.T @ De); o = np.argsort(-lam); lam = np.clip(lam[o], 0, None); V = V[:, o]; shv = lam / max(lam.sum(), 1e-30)
                ue = De @ V; uo = Do @ V; neff = nv * neff_per_vox
                rr = np.array([float(np.corrcoef(ue[:, k], uo[:, k])[0, 1]) if lam[k] > 0 else 0. for k in range(nk)])
                zz = np.arctanh(np.clip(rr, -0.999999, 0.999999)) * math.sqrt(max(neff - 3., 1.))
                keep = np.flatnonzero((zz >= ZFAC) & (shv >= fshare))
                log.append('comp %d (%d vox, n_eff %.0f): shares %s, r %s, z %s -> %d factor(s)' % (ic, nv, neff, ' '.join('%.2f' % x for x in shv[:4]), ' '.join('%.2f' % x for x in rr[:4]), ' '.join('%.1f' % x for x in zz[:4]), len(keep)))
                if len(keep) >= 2:
                    mag = np.abs(0.5 * (ue[:, keep] + uo[:, keep])); asg = keep[np.argmax(mag, axis=1)]
                    for j, kf in enumerate(keep, 1):
                        b = np.zeros((n, n, n), bool); b[tuple(ix[asg == kf] for ix in idx)] = True; lb, nb = ndimage.label(b); ps = np.bincount(lb.ravel(), minlength=nb + 1)
                        b = np.isin(lb, 1 + np.flatnonzero(ps[1:] >= max(20, DMINVOX // 4)))
                        if b.sum() == 0: continue
                        add(b, ic, j, prof=V[:, kf] / np.sqrt(wn))
                    l_fac = len(cands) > nc0
        if not l_fac: add(cm, ic, 0)
        if all(c['mol'] <= DMAXMOL for c in cands[nc0:]): continue
        cut2 = reg['cut'] + MFRAC * (vs[cm].max() - reg['cut']); lb, nb = ndimage.label(cm & (vs > cut2)); ps = np.bincount(lb.ravel(), minlength=nb + 1)[1:]
        log.append('comp %d not local: re-cut at %.3f x max -> %d pieces' % (ic, cut2 / vs.max(), int((ps >= DMINVOX).sum())))
        if (ps >= DMINVOX).sum() < 2: continue
        for k in np.flatnonzero(ps >= DMINVOX) + 1: np2 += 1; regs.append(dict(v=np.flatnonzero(lb.ravel() == k), cut=cut2, top=False, comp=np2))
    ord_ = sorted([i for i in range(len(cands)) if cands[i]['mol'] <= DMAXMOL], key=lambda i: (-cands[i]['share'], -cands[i]['r']))
    def dice(i, j): return 2. * float((cands[i]['on'] & cands[j]['on']).sum()) / max(1, cands[i]['vox'] + cands[j]['vox'])
    reps, memb = [], {}
    for i in ord_:
        k = 0
        for jr, rp in enumerate(reps, 1):
            if dice(i, rp) < DDICE: continue
            if any(cands[i]['fac'] > 0 and cands[it]['fac'] > 0 and cands[it]['comp'] == cands[i]['comp'] for it in memb if memb[it] == jr): continue
            k = jr; break
        if k == 0: reps.append(i); k = len(reps)
        memb[i] = k
    sites = [cands[i] for i in reps]
    desc = '%d sites from %d candidates (%d not local)' % (len(sites), len(cands), len(cands) - len(ord_))
    if regs_ok():
        best = []
        for r, nm in zip(REG, RNAME):
            cv = [cover(s['m'], r) for s in sites]; j = int(np.argmax(cv)) if cv else -1
            best.append('%s: site %d cov %.2f' % (nm, j + 1, cv[j] if cv else 0.))
        prec = [float(((s['m'] > 0.5) & ANY).sum()) / max(1, (s['m'] > 0.5).sum()) for s in sites]
        desc += '; recall >= 0.8: %d/%d, >= 0.5: %d/%d; sites >= 50 %% inside a region: %d/%d; %s' % (sum(1 for r in REG if sites and max(cover(s['m'], r) for s in sites) >= 0.8), len(REG),
            sum(1 for r in REG if sites and max(cover(s['m'], r) for s in sites) >= 0.5), len(REG), sum(p >= 0.5 for p in prec), len(sites), '; '.join(best))
    say('%s DISC %s: %s' % (a.tag, tag, desc))
    for l in log: say('      ' + l)
    say('      sites: ' + ' | '.join('%s share %.2f mol %.2f r %.3f seed %d' % (s['tag'], s['share'], s['mol'], s['r'], s['seed']) for s in sites))
    return sites
if 'recB0' in KEEP and not a.no_disc:
    K = KEEP['recB0']; geometry(*hdr(SETS['recB0'][0][0])); ts = os.path.join(a.node, 'sites', 'sites.txt')
    if os.path.exists(ts): say('%s DISC tree sites.txt: %d sites' % (a.tag, sum(1 for l in open(ts) if not l.startswith('#'))))
    discover(K['vs'], K['sup'], K['E'], K['O'], K['wn'], K['w'], K['names'], 'kmeans', 0., 'old', 0., 'current 2-means, old core')
    for fs in FSHARES: discover(K['vs'], K['sup'], K['E'], K['O'], K['wn'], K['w'], K['names'], 'svd', fs, 'old', 0., 'svd fshare %.2f, old core' % fs)
    for S in (0.6, 0.7): discover(K['vs'], K['sup'], K['E'], K['O'], K['wn'], K['w'], K['names'], 'svd', 0.15, 'share', S, 'svd fshare 0.15, share-%.1f core' % S)
say('done in %.0f s' % (time.time() - T0)); REP.close()
