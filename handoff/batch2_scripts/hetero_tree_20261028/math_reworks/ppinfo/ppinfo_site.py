#!/usr/bin/env python3
"""ppinfo_site.py: numpy port of site_rank's site EM (simple_state_mergeup::site_em and the hold-out acceptance of exec_site_rank)
with the per-particle information levers switchable: the box / band (--box, --emlim), the mask (--mask soft|core|oracle:...), the
poses (--poses input|stageB|true) and the labels (--labels em|truth|file:<post>). One call = one site under one configuration;
writes OUT/summary.json (per endpoint: N*, class share, hold-out peak r, n_eff, SE, z against HALVES_KSIG, rare AUC, purity and
recall of the class at posterior >= 0.5 and of the top-N* set, the deviation FSC per shell) and site_post_<tag>.txt in site_rank's
format (project-row order), so lever/site_score.py can cross-check the scores.
The port follows the Fortran step by step: planes = the shifted observation times the signed CTF on the native lattice of the
cropped box, nearest-neighbour class maps per half with the median-density regulariser, cross-half E-steps on FSC-weighted
site-masked deviations with per-shell scale and noise, the top-N hold-out curve of the masked deviation over shells kb+2..khi,
Fisher z at mask_neff. The conventions (origin at the box centre, R = Rz(e3) Ry(e2) Rz(e1) rows as the plane basis, CTF
sin(chi + asin(fraca))) are those of simple_ori_utils / simple_math_ctf; --geomcheck VOL verifies them against a saved map.
usage: ppinfo_site.py --node NODE --dump DUMPDIR --mask site_05.mrc --dirs 1,3 --out OUT [--truth truth.txt --rare 18]
       [--box 48 --emlim 15 --poses input --labels em --maskmode soft --cache planes.npy --nproc 4 --geomcheck VOL]"""
import argparse, json, os, sys, time
import numpy as np
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mrcio import read_mrc, stack, fcrop3

SITE_GRID = [500, 1000, 1500, 2500, 4000, 6000, 10000, 15000, 25000, 40000, 60000, 90000]
SITE_ITS, SITE_ITS_MIN, SITE_CONV, SITE_PFLOOR, RB_MARGIN, HALVES_KSIG, COSMSKHALFWIDTH = 50, 10, 0.005, 1e-3, 2, 3., 6.
KEYS = ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast', 'phshift', 'state', 'eo', 'indstk', 'stkind', 'contrast')
G = {}   # the shared read-only arrays of the forked start workers


def log(*a):
    print('[%s] %s' % (time.strftime('%H:%M:%S'), ' '.join(str(x) for x in a)), flush=True)


def read_dump(fn):
    """print_project_field text -> dict of arrays (one row per record; missing keys 0)."""
    out = {k: [] for k in KEYS}
    for ln in open(fn):
        if 'e1=' not in ln:
            continue
        d = dict(t.split('=', 1) for t in ln.split() if '=' in t)
        for k in KEYS:
            out[k].append(float(d.get(k, 0.)))
    return {k: np.array(v) for k, v in out.items()}


def read_stk(fn):
    """the stk segment: one record per stack (path, fromp, top); box, smpd and the microscope from the first."""
    recs = []
    for ln in open(fn):
        if 'stk=' not in ln:
            continue
        d = dict(t.split('=', 1) for t in ln.split() if '=' in t)
        recs.append((d['stk'], int(float(d.get('fromp', 1))), int(float(d.get('top', 0)))))
    d = dict(t.split('=', 1) for t in open(fn).readline().split() if '=' in t)
    return recs, int(float(d['box'])), float(d['smpd']), float(d['kv']), float(d['cs']), float(d['fraca'])


def read_oritab(fn):
    """an oritab with e1 e2 e3 (and x y) per line, in stack order: the true poses of the pool."""
    e = []
    for ln in open(fn):
        d = dict(t.split('=', 1) for t in ln.split() if '=' in t)
        if 'e1' in d:
            e.append([float(d['e1']), float(d['e2']), float(d['e3']), float(d.get('x', 0.)), float(d.get('y', 0.))])
    return np.array(e)


def rotmat(ang, axis):
    c, s = np.cos(np.radians(ang)), np.sin(np.radians(ang))
    r = np.zeros((len(ang), 3, 3))
    if axis == 2:
        r[:, 0, 0] = c; r[:, 0, 2] = -s; r[:, 1, 1] = 1.; r[:, 2, 0] = s; r[:, 2, 2] = c
    else:
        r[:, 0, 0] = c; r[:, 0, 1] = s; r[:, 1, 0] = -s; r[:, 1, 1] = c; r[:, 2, 2] = 1.
    return r


def euler2m(e):   # simple_ori_utils::euler2m: r3 r2 r1
    return rotmat(e[:, 2], 3) @ rotmat(e[:, 1], 2) @ rotmat(e[:, 0], 3)


def soft_sphere(B, rad):
    """simple_image::mask3D_soft(rad, backgr=0): cosine edge of half-width COSMSKHALFWIDTH about rad, origin at B/2."""
    x = np.arange(B) - B // 2
    r = np.sqrt(x[:, None, None] ** 2 + x[None, :, None] ** 2 + x[None, None, :] ** 2)
    hw = COSMSKHALFWIDTH
    w = np.clip(0.5 * (1. + np.cos(np.pi * (r - (rad - hw)) / (2. * hw))), 0., 1.)
    w[r <= rad - hw] = 1.
    w[r >= rad + hw] = 0.
    return w


def crop_mask(fn, smpd_full, B):
    v, smpd = read_mrc(fn)
    if v.shape[0] != B:
        v = fcrop3(v, B)
    return np.clip(v, 0., 1.)


def prep_crop(recs, box, stkind, indstk, msk_px, Bc, cache):
    """the Fourier planes of the particles at the native lattice of box Bc (origin at the box centre, background mean removed,
    no shift, no CTF): crop[n, Bc (k, fft order), Bc/2 (h >= 0)] complex64; cached as .npy. recs = the stk records, the particle
    i is image indstk[i] of stack stkind[i]."""
    n = len(indstk)
    if cache and os.path.exists(cache):
        c = np.load(cache, mmap_mode='r')
        if c.shape == (n, Bc, Bc // 2):
            log('planes from', cache)
            return np.ascontiguousarray(c)
    out = np.empty((n, Bc, Bc // 2), np.complex64)
    yy, xx = np.mgrid[:box, :box]
    bg = np.sqrt((yy - box // 2) ** 2 + (xx - box // 2) ** 2) > msk_px
    rows = np.r_[0:Bc // 2, box - Bc // 2:box]
    t0 = time.time()
    order = np.lexsort((indstk, stkind))
    for s0 in range(0, n, 2000):
        sel = order[s0:s0 + 2000]
        img = np.empty((len(sel), box, box), np.float64)
        for sk in np.unique(stkind[sel]):
            w = np.where(stkind[sel] == sk)[0]
            st = stack(recs[sk - 1][0])
            img[w] = np.asarray(st[indstk[sel[w]] - 1], np.float32)
        img -= img[:, bg].mean(1)[:, None, None]
        F = np.fft.rfft2(np.fft.ifftshift(img, axes=(1, 2)))
        out[sel] = F[:, rows, :Bc // 2].astype(np.complex64)
    log('planes prepared for %d particles at box %d in %.0f s' % (n, Bc, time.time() - t0))
    if cache:
        np.save(cache, out)
    return out


def ctf_eval(h, k, dfx, dfy, angast, phshift, kv, cs, fraca, box_A):
    """simple_math_ctf::ft_map_ctf_kernel in Angstrom units: sin(pi wl s2 (df - 0.5 wl^2 s2 Cs) + phshift + asin(fraca));
    h, k [m], the particle parameters [n] -> [n, m]."""
    V = kv * 1000.
    wl = 12.26 / np.sqrt(V + 0.9784 * V * V / 1e6)
    s2 = (h * h + k * k)[None, :] / box_A ** 2
    ang = np.arctan2(k, h)[None, :]
    df = 0.5 * ((dfx + dfy)[:, None] * 1e4 + np.cos(2. * (ang - np.radians(angast)[:, None])) * (dfx - dfy)[:, None] * 1e4)
    chi = np.pi * wl * s2 * (df - 0.5 * wl * wl * s2 * cs * 1e7)
    return np.sin(chi + phshift[:, None] + np.arcsin(fraca)).astype(np.float32)


def auc(score, pos):
    o = np.argsort(score, kind='stable'); r = np.empty(len(score)); r[o] = np.arange(1, len(score) + 1)
    npos = pos.sum(); nneg = len(pos) - npos
    return float((r[pos].sum() - npos * (npos + 1) / 2.) / max(npos * nneg, 1))


def setup(a):
    """build the shared arrays: planes P (CTF x shifted observation), C (CTF^2), f2 (data power per shell), the cell indices of
    every plane sample under every particle's rotation, the shells and weights of the cells, the masks."""
    nd = read_dump(os.path.join(a.dump, 'n_ptcl3D.txt'))
    recs, box, smpd, kv, cs, fraca = read_stk(os.path.join(a.dump, 'n_stk.txt'))
    n0 = len(nd['e1'])
    rows = np.where((nd['state'] > 0) & ((nd['eo'] == 0) | (nd['eo'] == 1)))[0]
    n = len(rows)
    G['n0'], G['rows'], G['n'] = n0, rows, n
    indstk = nd['indstk'][rows].astype(int)
    stkind = np.maximum(nd['stkind'][rows].astype(int), 1)
    eo = nd['eo'][rows].astype(int)
    log('%d stacks; %d of %d rows' % (len(recs), n, n0))
    # exec_site_rank: box_crop = min(box, max(26, 2 ceil(box smpd / lp))), lp 10 by default
    B = a.box if a.box > 0 else min(box, max(26, 2 * int(np.ceil(box * smpd / a.lp))))
    smpd_crop = box * smpd / B
    box_A = box * smpd
    nyq, nh, khi = B // 2, B // 2 + 1, B // 2 - 1
    kb = min(khi, int(round(box_A / a.emlim)))
    klo = max(int(round(box_A / 40.)), 3)
    msk_px = a.mskdiam / 2. / smpd
    msk_crop = a.mskdiam / 2. / smpd_crop
    log('node %s: %d of %d rows, box %d smpd %.4f -> box_crop %d smpd_crop %.4f, EM shells 2..%d, test band %d..%d, hold-out %d..%d'
        % (a.node, n, n0, box, smpd, B, smpd_crop, kb, klo, khi, kb + RB_MARGIN, khi))
    # poses
    if a.poses == 'input':
        e = np.c_[nd['e1'], nd['e2'], nd['e3'], nd['x'], nd['y']][rows]
    elif a.poses == 'stageB':
        bd = read_dump(os.path.join(a.dump, 'b_ptcl3D.txt'))
        assert len(bd['e1']) == n0, 'stage-B dump has %d rows for %d' % (len(bd['e1']), n0)
        e = np.c_[bd['e1'], bd['e2'], bd['e3'], bd['x'], bd['y']][rows]
    elif a.poses.startswith('true:'):
        t = read_oritab(a.poses[5:])
        e = t[indstk - 1]
    else:
        raise SystemExit('poses: input | stageB | true:<oritab>')
    R = euler2m(e[:, :3])
    if a.posenoise > 0:
        # a random rotation of angle ~ |N(0, posenoise deg)| about a random axis composed on every pose, shifts jittered by 0.4 px/deg
        rng = np.random.default_rng(7)
        ax = rng.normal(size=(n, 3)); ax /= np.linalg.norm(ax, axis=1)[:, None]
        th = np.radians(np.abs(rng.normal(0., a.posenoise, n)))
        K = np.zeros((n, 3, 3)); K[:, 0, 1] = -ax[:, 2]; K[:, 0, 2] = ax[:, 1]; K[:, 1, 0] = ax[:, 2]; K[:, 1, 2] = -ax[:, 0]; K[:, 2, 0] = -ax[:, 1]; K[:, 2, 1] = ax[:, 0]
        Rn = np.eye(3)[None] + np.sin(th)[:, None, None] * K + (1 - np.cos(th))[:, None, None] * (K @ K)
        R = Rn @ R
        e = e.copy(); e[:, 3:5] += rng.normal(0., 0.4 * a.posenoise, (n, 2))
        log('poses perturbed: rotation |N(0, %.1f deg)| about random axes (median %.2f deg), shifts N(0, %.2f px)' % (a.posenoise, np.degrees(np.median(th)), 0.4 * a.posenoise))
    R = R.astype(np.float32)
    # planes
    Bc = max(B, a.cachebox)
    assert Bc >= B, 'cache box below the working box'
    crop = prep_crop(recs, box, stkind, indstk, msk_px, Bc, a.cache)
    if Bc != B:
        crop = np.ascontiguousarray(crop[:, np.r_[0:B // 2, Bc - B // 2:Bc], :B // 2])
    # plane samples: h >= 0 half, 3 <= h^2 + k^2 <= khi (khi + 1); weight 1 at h = 0 (both k signs present), 2 otherwise
    kk, hh = np.meshgrid(np.r_[0:B // 2, -B // 2:0], np.arange(B // 2), indexing='ij')
    r2 = hh * hh + kk * kk
    sel = (r2 >= 3) & (r2 <= khi * (khi + 1))
    h, k = hh[sel].astype(np.int32), kk[sel].astype(np.int32)
    krow = np.where(k >= 0, k, k + B)
    sw = np.where(h == 0, 1., 2.).astype(np.float32)
    shell = (np.rint(np.sqrt(h * h + k * k)).astype(int) - 1)          # emsh / fbsh: ring nint(r) - 1 (1..khi-1)
    em = (h * h + k * k) <= kb * (kb + 1)
    nsh = shell.max() + 1
    # X: shifted observation; P = CTF X; f2 = per-shell |X|^2
    X = crop[:, krow, h]
    del crop
    ph = np.exp(-2j * np.pi * a.shiftsign * (h[None, :] * e[:, 3:4] + k[None, :] * e[:, 4:5]) / box).astype(np.complex64)
    X *= ph
    del ph
    ctf = ctf_eval(h, k, nd['dfx'][rows], nd['dfy'][rows], nd['angast'][rows], nd['phshift'][rows], kv, cs, fraca, box_A)
    f2 = np.zeros((n, nsh), np.float64)
    for s in range(nsh):
        m = shell == s
        if m.any():
            f2[:, s] = (sw[m][None, :] * (X[:, m].real ** 2 + X[:, m].imag ** 2)).sum(1)
    X *= ctf
    cn = np.ones(n, np.float32)
    cf = os.path.join(a.site, 'site_contrast.txt') if a.site else ''
    if a.contrast and os.path.exists(cf):
        t = np.loadtxt(cf, ndmin=2)
        cmap = np.zeros(n0 + 1); cmap[t[:, 0].astype(int)] = t[:, 1]
        cn = cmap[rows + 1].astype(np.float32)
        log('contrast weights from site_contrast.txt: median %.3f' % np.median(cn))
    # geometry: p = h R[0,:] + k R[1,:] rounded; Friedel to h >= 0; cell (l, k, h) of the half volume V[B, B, nh]
    NC = B * B * nh
    idx = np.empty((n, len(h)), np.int32)
    cnj = np.empty((n, len(h)), bool)
    mate = np.full((n, len(h)), -1, np.int32)
    for b0 in range(0, n, 4000):
        b1 = min(n, b0 + 4000)
        p = np.rint(h[None, :, None] * R[b0:b1, None, 0, :] + k[None, :, None] * R[b0:b1, None, 1, :]).astype(np.int32)
        c = p[..., 0] < 0
        p[c] = -p[c]
        ok = (np.abs(p) < nyq).all(-1)
        ii = (p[..., 2] % B) * B * nh + (p[..., 1] % B) * nh + p[..., 0]
        idx[b0:b1] = np.where(ok, ii, NC)
        cnj[b0:b1] = c
        m0 = ok & (p[..., 0] == 0)
        mate[b0:b1] = np.where(m0, ((-p[..., 2]) % B) * B * nh + ((-p[..., 1]) % B) * nh, -1)
    # cell shells and Friedel weights
    l = np.arange(B); ls = np.where(l < nyq, l, l - B)
    hc = np.arange(nh)
    cs_ = np.rint(np.sqrt(ls[:, None, None] ** 2 + ls[None, :, None] ** 2 + hc[None, None, :] ** 2))
    cell_shell = np.minimum(cs_, nyq).astype(np.int32).ravel()
    cell_fw = np.where((hc > 0) & (hc < nh - 1), 2., 1.)[None, None, :].repeat(B, 0).repeat(B, 1).ravel()
    # masks: the EM / hold-out mask sm = site mask x soft sphere; the sphere alone for mask_neff's reference
    sph = soft_sphere(B, msk_crop)
    if a.maskmode == 'soft':
        sm = crop_mask(a.mask, smpd, B) * sph
    elif a.maskmode == 'core':
        sm = (crop_mask(a.mask, smpd, B) >= 0.5).astype(float) * sph
    elif a.maskmode.startswith('oracle:'):
        # the true rare-minus-majority difference, low-passed to the EM limit, thresholded at a fraction of its peak
        fa, fb, frac = a.maskmode[7:].split(',')
        va, _ = read_mrc(fa); vb, _ = read_mrc(fb)
        d = fcrop3(va - vb, B)
        d = np.abs(d)
        core = d >= float(frac) * d.max()
        from scipy.ndimage import binary_dilation, gaussian_filter
        sm = gaussian_filter(binary_dilation(core, iterations=1).astype(float), 1.0) * sph
        sm = np.clip(sm / max(sm.max(), 1e-9), 0., 1.)
    elif a.maskmode == 'sphere':
        sm = sph.copy()
    else:
        raise SystemExit('maskmode: soft | core | sphere | oracle:<rare.mrc>,<maj.mrc>,<frac>')
    frac_sphere = float((sm > 0.5).sum()) / max(1, (sph > 0.5).sum())
    # mask_neff over the hold-out band
    s2, s4 = (sm ** 2).sum(), (sm ** 4).sum()
    sh3 = np.minimum(np.rint(np.sqrt(ls[:, None, None] ** 2 + ls[None, :, None] ** 2 + ls[None, None, :] ** 2)), nyq)
    nband = int(((sh3 >= kb + RB_MARGIN) & (sh3 <= khi)).sum())
    neff = 0.5 * nband * s2 * s2 / max(s4, 1e-30) / B ** 3
    log('mask %s (%s): %.1f %% of the sphere, n_eff %.1f over hold-out shells %d..%d (bar r %.3f)'
        % (os.path.basename(a.mask), a.maskmode, 100. * frac_sphere, neff, kb + RB_MARGIN, khi, np.tanh(HALVES_KSIG / np.sqrt(max(neff - 3., 1.)))))
    G.update(dict(P=X, C=ctf * ctf, f2=f2, cn=cn, idx=idx, cnj=cnj, mate=mate, sw=sw, shell=shell, em=em, nsh=nsh, NC=NC, B=B, nh=nh, nyq=nyq,
                  kb=kb, khi=khi, cell_shell=cell_shell, cell_fw=cell_fw, sm_c=np.fft.ifftshift(sm), eo=eo, indstk=indstk, e=e,
                  smpd_crop=smpd_crop, neff=neff, frac_sphere=frac_sphere, nband=nband))
    return nd


# ---------------------------------------------------------------- nearest-neighbour maps

def nn_acc(pidx, pw, sub):
    """num / den of the nearest-neighbour insertion of the samples sub of the particles pidx with weights pw (Fortran nn_acc)."""
    P, C, cn, idx, cnj, mate, sw, NC = G['P'], G['C'], G['cn'], G['idx'], G['cnj'], G['mate'], G['sw'], G['NC']
    num = np.zeros(NC + 1, np.complex128)
    den = np.zeros(NC + 1, np.float64)
    for b0 in range(0, len(pidx), 8000):
        pi = pidx[b0:b0 + 8000]
        w = pw[b0:b0 + 8000]
        if not np.any(w > 0):
            continue
        ii = idx[pi][:, sub]
        v = P[pi][:, sub] * (w[:, None] * cn[pi][:, None] * sw[None, sub])
        v = np.where(cnj[pi][:, sub], np.conj(v), v)
        d = C[pi][:, sub] * (w[:, None] * cn[pi][:, None] ** 2 * sw[None, sub])
        mm = mate[pi][:, sub]
        half = mm >= 0
        wf = np.where(half, 0.5, 1.)
        f = ii.ravel()
        num.real += np.bincount(f, (v.real * wf).ravel(), NC + 1)
        num.imag += np.bincount(f, (v.imag * wf).ravel(), NC + 1)
        den += np.bincount(f, (d * wf).ravel(), NC + 1)
        if half.any():
            fm = mm[half]
            vm = v[half]
            num.real += np.bincount(fm, 0.5 * vm.real, NC + 1)
            num.imag += np.bincount(fm, -0.5 * vm.imag, NC + 1)
            den += np.bincount(fm, 0.5 * d[half], NC + 1)
    return num[:NC], den[:NC]


def density(pidx, sub):
    """per-shell median of the unit-weight den per particle (Fortran density)."""
    _, den = nn_acc(pidx, np.ones(len(pidx), np.float32), sub)
    eps = np.ones(G['nyq'] + 1)
    for s in range(G['nyq'] + 1):
        v = den[(G['cell_shell'] == s) & (den > 0)]
        if len(v):
            eps[s] = np.sort(v)[(len(v) + 1) // 2 - 1] / len(pidx)
    return eps


def nn_norm(num, den, eps, wsum):
    return (num / (den + eps[G['cell_shell']] * max(wsum, 1.))).astype(np.complex64)


def to_real(V):
    return np.fft.irfftn(V.reshape(G['B'], G['B'], G['nh']), s=(G['B'],) * 3)


def masked(V):
    return np.fft.rfftn(to_real(V) * G['sm_c']).astype(np.complex64).ravel()


def shell_weighted(V, f):
    return (V * f[G['cell_shell']]).astype(np.complex64)


def shell_fsc(x, y):
    fw, cs = G['cell_fw'], G['cell_shell']
    sxy = np.bincount(cs, fw * (x.real * y.real + x.imag * y.imag), G['nyq'] + 1)
    sxx = np.bincount(cs, fw * (np.abs(x) ** 2), G['nyq'] + 1)
    syy = np.bincount(cs, fw * (np.abs(y) ** 2), G['nyq'] + 1)
    return sxy / np.maximum(np.sqrt(sxx * syy), 1e-30)


def shell_corr_dot(x, y, shlo, shhi):
    m = (G['cell_shell'] >= shlo) & (G['cell_shell'] <= shhi)
    fw = G['cell_fw'][m]
    sxy = float((fw * (x[m].real * y[m].real + x[m].imag * y[m].imag)).sum())
    sxx = float((fw * np.abs(x[m]) ** 2).sum()); syy = float((fw * np.abs(y[m]) ** 2).sum())
    return sxy / max(np.sqrt(sxx * syy), 1e-30), sxy


# ---------------------------------------------------------------- the E-step

def nn_terms(pidx, refs, sub):
    """aa[nsh, 2, np], bb[nsh, 2, np] of the particles pidx against the two reference maps (Fortran nn_terms)."""
    P, C, cn, idx, cnj, sw, shell, nsh = G['P'], G['C'], G['cn'], G['idx'], G['cnj'], G['sw'], G['shell'], G['nsh']
    npart = len(pidx)
    aa = np.zeros((nsh, 2, npart)); bb = np.zeros((nsh, 2, npart))
    shs = shell[sub]
    for b0 in range(0, npart, 8000):
        pi = pidx[b0:b0 + 8000]
        ii = idx[pi][:, sub]
        ok = ii < G['NC']
        ii = np.where(ok, ii, 0)
        c = cnj[pi][:, sub]
        dv = P[pi][:, sub]
        c2v = C[pi][:, sub]
        key = (np.arange(len(pi))[:, None] * nsh + shs[None, :]).ravel()
        wc = (cn[pi][:, None] * sw[None, sub]) * ok
        for g in range(2):
            pv = refs[g][ii]
            pv = np.where(c, np.conj(pv), pv)
            a = (pv.real * dv.real + pv.imag * dv.imag) * wc
            b = (pv.real ** 2 + pv.imag ** 2) * c2v * wc * cn[pi][:, None]
            aa[:, g, b0:b0 + len(pi)] = np.bincount(key, a.ravel(), len(pi) * nsh).reshape(len(pi), nsh).T
            bb[:, g, b0:b0 + len(pi)] = np.bincount(key, b.ravel(), len(pi) * nsh).reshape(len(pi), nsh).T
    return aa, bb


def nn_post(aa, bb, beta, sig2, pr, emsh):
    s = 0.5 * ((2. * beta[:, None, None] * aa[emsh] - beta[:, None, None] ** 2 * bb[emsh]) / sig2[:, None, None]).sum(0) + pr[:, None]
    s -= s.max(0)[None, :]
    w = np.exp(s)
    return (w / w.sum(0)[None, :]).T


def nn_estep(pidx, w, refs, sub, emsh, it, tag):
    aa, bb = nn_terms(pidx, refs, sub)
    mpc = np.argmax(w, 1)
    pr = np.log(np.maximum(w.sum(0) / len(pidx), SITE_PFLOOR))
    ar = np.arange(len(pidx))
    sa = aa[emsh][:, mpc, ar].sum(1); sb = bb[emsh][:, mpc, ar].sum(1)
    sf = G['f2'][pidx][:, emsh].sum(0)
    cw = np.bincount(G['shell'][sub], G['sw'][sub], G['nsh'])[emsh]
    sc = len(pidx) * cw
    beta = np.clip(sa / np.maximum(sb, 1e-30), 0., 2.)
    sig2 = np.maximum((sf - 2. * beta * sa + beta * beta * sb) / np.maximum(sc, 1.), 1e-12)
    wn = nn_post(aa, bb, beta, sig2, pr, emsh)
    chg = float(((wn[:, 1] > 0.5) != (w[:, 1] > 0.5)).mean())
    llr = 0.5 * ((2. * beta[:, None] * (aa[emsh][:, 1, :] - aa[emsh][:, 0, :]) - beta[:, None] ** 2 * (bb[emsh][:, 1, :] - bb[emsh][:, 0, :])) / sig2[:, None]).sum(0)
    log('  %s iteration %d: argmax changed %.3f, site class weight %.3f, model scale at shells 2 / %d %.2f / %.2f'
        % (tag, it, chg, wn[:, 1].mean(), G['kb'], beta[0], beta[-1]))
    return wn, chg, llr


# ---------------------------------------------------------------- one start

def site_em(job):
    tag, seed, given = job
    n, eo, kb, khi, nyq = G['n'], G['eo'], G['kb'], G['khi'], G['nyq']
    sub_em = np.where(G['em'])[0]
    sub_f = np.arange(len(G['sw']))
    emsh = np.arange(1, kb)          # EM shells 1..kb-1 (rings nint(r) 2..kb)
    ia = np.where(eo == 0)[0]; ib = np.where(eo == 1)[0]
    t0 = time.time()
    post = np.zeros(n)
    llr = np.zeros(n)
    nit = 0
    if given is None:
        epsb = 0.5 * (density(ia, sub_em) + density(ib, sub_em))
        wa = np.c_[1. - seed[ia], seed[ia]].astype(np.float64); wb = np.c_[1. - seed[ib], seed[ib]].astype(np.float64)
        nma, dna = nn_acc(ia, np.ones(len(ia), np.float32), sub_em)
        nmb, dnb = nn_acc(ib, np.ones(len(ib), np.float32), sub_em)
        all_a = nn_norm(nma, dna, epsb, len(ia)); all_b = nn_norm(nmb, dnb, epsb, len(ib))
        for it in range(1, SITE_ITS + 1):
            n2, d2 = nn_acc(ia, wa[:, 1].astype(np.float32), sub_em)
            mp_a = [nn_norm(nma - n2, dna - d2, epsb, wa[:, 0].sum()), nn_norm(n2, d2, epsb, wa[:, 1].sum())]
            n2, d2 = nn_acc(ib, wb[:, 1].astype(np.float32), sub_em)
            mp_b = [nn_norm(nmb - n2, dnb - d2, epsb, wb[:, 0].sum()), nn_norm(n2, d2, epsb, wb[:, 1].sum())]
            rf_a, rf_b = [], []
            for g in range(2):
                f = np.clip(shell_fsc(mp_a[g] - all_a, mp_b[g] - all_b), 0., 1.)
                rf_a.append(all_a + masked(shell_weighted(mp_a[g] - all_a, f)))
                rf_b.append(all_b + masked(shell_weighted(mp_b[g] - all_b, f)))
            wa, chga, la = nn_estep(ia, wa, rf_b, sub_em, emsh, it, tag + ' even')
            wb, chgb, lb = nn_estep(ib, wb, rf_a, sub_em, emsh, it, tag + ' odd')
            nit = it
            if it >= G['maxits'] or (it >= SITE_ITS_MIN and max(chga, chgb) < SITE_CONV):
                break
        post[ia] = wa[:, 1]; post[ib] = wb[:, 1]
        llr[ia] = la; llr[ib] = lb
        log('  %s: site EM %s after %d iterations (%.0f s)' % (tag, 'converged' if nit < G['maxits'] else 'stopped at the cap', nit, time.time() - t0))
        fsc_cls = np.clip(shell_fsc(mp_a[1] - all_a, mp_b[1] - all_b), -1., 1.)
    else:
        post[:] = given
        fsc_cls = None
    # hold-out curve: the top-N/2 of each half, full-band nearest-neighbour maps, masked deviation from the half's all-map
    epsf = 0.5 * (density(ia, sub_f) + density(ib, sub_f))
    nma, dna = nn_acc(ia, np.ones(len(ia), np.float32), sub_f)
    nmb, dnb = nn_acc(ib, np.ones(len(ib), np.float32), sub_f)
    tm1 = nn_norm(nma, dna, epsf, len(ia)); tm2 = nn_norm(nmb, dnb, epsf, len(ib))
    oa = ia[np.argsort(-post[ia], kind='stable')]; ob = ib[np.argsort(-post[ib], kind='stable')]
    grid = [g for g in SITE_GRID if g <= int(0.7 * n)]
    curve = np.zeros((4, len(grid)))
    npa = npb = 0
    nma = np.zeros(G['NC'], np.complex128); dna = np.zeros(G['NC']); nmb = nma.copy(); dnb = dna.copy()
    best = (-np.inf, None, None, grid[0])
    fsc_dev = None
    for ig, N in enumerate(grid):
        nh2 = N // 2
        ka, kb_ = min(nh2, len(ia)), min(nh2, len(ib))
        if ka > npa:
            n2, d2 = nn_acc(oa[npa:ka], np.ones(ka - npa, np.float32), sub_f); nma += n2; dna += d2; npa = ka
        if kb_ > npb:
            n2, d2 = nn_acc(ob[npb:kb_], np.ones(kb_ - npb, np.float32), sub_f); nmb += n2; dnb += d2; npb = kb_
        da = masked(nn_norm(nma, dna, epsf, npa) - tm1)
        db = masked(nn_norm(nmb, dnb, epsf, npb) - tm2)
        r, dot = shell_corr_dot(da, db, kb + RB_MARGIN, khi)
        curve[:, ig] = [N, r, np.sqrt(max(dot, 0.)), 0.]
        log('  %s N %d: hold-out cross-half correlation of the site deviation %.3f' % (tag, N, r))
        if r > best[0]:
            best = (r, da, db, N)
            fsc_dev = shell_fsc(da, db)
    curve[3] = curve[2] * curve[0] / curve[2, 0] if curve[2, 0] > 0 else best[3]
    r, da, db, nstar = best
    # the pole rule's seed test: the particles at posterior >= 0.5, their masked deviation per half, cross-half r over the hold-out band
    hi = post >= 0.5
    sa, sb = np.where(hi & (eo == 0))[0], np.where(hi & (eo == 1))[0]
    seed_r = seed_z = 0.
    if len(sa) >= 3 and len(sb) >= 3:
        ea = masked(nn_norm(*nn_acc(sa, np.ones(len(sa), np.float32), sub_f), epsf, len(sa)) - tm1)
        eb = masked(nn_norm(*nn_acc(sb, np.ones(len(sb), np.float32), sub_f), epsf, len(sb)) - tm2)
        seed_r = shell_corr_dot(ea, eb, kb + RB_MARGIN, khi)[0]
        seed_z = float(np.arctanh(np.clip(seed_r, -0.999999, 0.999999)) * np.sqrt(max(G['neff'] - 3., 1.)))
    neff = G['neff']
    z = float(np.arctanh(np.clip(r, -0.999999, 0.999999)) * np.sqrt(max(neff - 3., 1.)))
    se = (1. - r * r) / np.sqrt(max(neff - 3., 1.))
    share = float(post.mean())
    log('  %s: N* %d (hold-out shells %d..%d, correlation %.3f, z %.2f at n_eff %.1f), class share %.3f, equivalent pure members %d'
        % (tag, nstar, kb + RB_MARGIN, khi, r, z, neff, share, int(round(curve[3].max()))))
    log('  %s: seeds at posterior >= 0.5: %d (%d / %d per half), their deviation r %.3f z %.2f' % (tag, int(hi.sum()), len(sa), len(sb), seed_r, seed_z))
    return dict(tag=tag, post=post, llr=llr, nstar=int(nstar), r=float(r), z=z, se=float(se), share=share, curve=curve.tolist(),
                seed_n=int(hi.sum()), seed_r=float(seed_r), seed_z=seed_z,
                epm=float(curve[3].max()), nit=nit, fsc_dev=None if fsc_dev is None else fsc_dev.tolist(),
                fsc_cls=None if fsc_cls is None else fsc_cls.tolist(), dev=(da, db), wall=time.time() - t0)


def geomcheck(vol_fn):
    """FSC of the all-particle nearest-neighbour maps against a saved map of the node, and between the two halves."""
    v, smpd = read_mrc(vol_fn)
    B = G['B']
    if v.shape[0] != B:
        v = fcrop3(v, B)
    V = np.fft.rfftn(np.fft.ifftshift(v)).astype(np.complex64).ravel()
    ia = np.where(G['eo'] == 0)[0]; ib = np.where(G['eo'] == 1)[0]
    sub = np.arange(len(G['sw']))
    eps = 0.5 * (density(ia, sub) + density(ib, sub))
    ma = nn_norm(*nn_acc(ia, np.ones(len(ia), np.float32), sub), eps, len(ia))
    mb = nn_norm(*nn_acc(ib, np.ones(len(ib), np.float32), sub), eps, len(ib))
    f1, f2 = shell_fsc(ma, V), shell_fsc(ma, mb)
    res = [B * G['smpd_crop'] / max(s, 1) for s in range(G['nyq'] + 1)]
    log('geomcheck vs %s: shell  res_A  fsc(even map, vol)  fsc(even, odd)' % vol_fn)
    for s in range(1, G['nyq'] + 1):
        log('  %2d %6.1f %7.3f %7.3f' % (s, res[s], f1[s], f2[s]))
    return f1.tolist(), f2.tolist()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--node', required=True); ap.add_argument('--dump', required=True); ap.add_argument('--mask', required=True)
    ap.add_argument('--out', required=True); ap.add_argument('--site', default='', help='the saved site_rank dir (site_contrast.txt, saved posteriors)')
    ap.add_argument('--dirs', default='1,3'); ap.add_argument('--latent', default='')
    ap.add_argument('--truth', default=''); ap.add_argument('--rare', type=int, default=0)
    ap.add_argument('--box', type=int, default=0, help='box_crop; 0 = the Fortran formula from --lp'); ap.add_argument('--lp', type=float, default=10.)
    ap.add_argument('--emlim', type=float, default=15.)
    ap.add_argument('--mskdiam', type=float, default=225.456)
    ap.add_argument('--poses', default='input'); ap.add_argument('--maskmode', default='soft')
    ap.add_argument('--labels', default='em', help='em | truth | seedtruth | file:<site_post.txt> (truth: no EM, the hold-out curve of the rare indicator; seedtruth: the EM started from the truth, with --maxits 1 one M + E step = the posterior of perfect references)')
    ap.add_argument('--maxits', type=int, default=SITE_ITS)
    ap.add_argument('--cache', default=''); ap.add_argument('--cachebox', type=int, default=0)
    ap.add_argument('--contrast', type=int, default=1); ap.add_argument('--nproc', type=int, default=4)
    ap.add_argument('--geomcheck', default=''); ap.add_argument('--shiftsign', type=float, default=1.)
    ap.add_argument('--posenoise', type=float, default=0., help='perturb the poses: rotation |N(0, deg)| about a random axis, shifts 0.4 px/deg')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    nd = setup(a)
    G['maxits'] = max(1, a.maxits)
    n, n0, rows = G['n'], G['n0'], G['rows']
    summary = dict(args=vars(a), n=n, n0=n0, box=G['B'], kb=G['kb'], khi=G['khi'], neff=G['neff'], frac_sphere=G['frac_sphere'],
                   nband=G['nband'], smpd_crop=G['smpd_crop'], bar_r=float(np.tanh(HALVES_KSIG / np.sqrt(max(G['neff'] - 3., 1.)))))
    if a.geomcheck:
        summary['geomcheck'] = geomcheck(a.geomcheck)
    # truth
    israre = None
    if a.truth and a.rare > 0:
        t = np.loadtxt(a.truth, ndmin=2)
        conf = np.zeros(int(t[:, 0].max()) + 1, int); conf[t[:, 0].astype(int)] = t[:, 7].astype(int)
        israre = conf[np.minimum(G['indstk'], len(conf) - 1)] == a.rare
        summary['nrare'] = int(israre.sum())
        log('truth: %d rare of %d (%.2f %%)' % (israre.sum(), n, 100. * israre.mean()))
    # starts
    jobs = []
    if a.labels == 'em':
        lat = a.latent or os.path.join(a.node, 'fflex', 'flex_pca_coordinates.txt')
        dirs = [int(x) for x in a.dirs.split(',')]
        z = np.zeros((n0 + 1, max(dirs)))
        for ln in open(lat):
            if ln.startswith('#') or not ln.split():
                continue
            w = ln.split()
            pind = int(w[0])
            if 1 <= pind <= n0:
                z[pind] = [float(x) for x in w[5:5 + max(dirs)]]
        ntail = min(15000, n // 6)
        for d in dirs:
            key = z[rows + 1, d - 1]
            o = np.argsort(key, kind='stable')
            lo = np.zeros(n, np.float32); lo[o[:ntail]] = 1.
            hi = np.zeros(n, np.float32); hi[o[n - ntail:]] = 1.
            jobs += [('z%dlo' % d, lo, None), ('z%dhi' % d, hi, None)]
    elif a.labels == 'truth':
        assert israre is not None, '--labels truth needs --truth/--rare'
        jobs = [('truth', None, israre.astype(float))]
    elif a.labels == 'seedtruth':
        assert israre is not None, '--labels seedtruth needs --truth/--rare'
        jobs = [('truth%d' % a.maxits, israre.astype(np.float32), None)]
    elif a.labels.startswith('file:'):
        for fn in a.labels[5:].split(','):
            pv = np.loadtxt(fn, ndmin=1)
            assert len(pv) == n0
            jobs.append((os.path.basename(fn).replace('site_post_', '').replace('.txt', ''), None, pv[rows]))
    else:
        raise SystemExit('labels: em | truth | file:<post>[,<post>]')
    log('%d starts; %d processes' % (len(jobs), a.nproc))
    if a.nproc > 1 and len(jobs) > 1:
        with Pool(min(a.nproc, len(jobs))) as pool:
            res = pool.map(site_em, jobs, chunksize=1)
    else:
        res = [site_em(j) for j in jobs]
    # endpoint table, cross-endpoint Gram over the hold-out band
    ends = []
    devs = []
    for r in res:
        post = r['post']
        e = {k: r[k] for k in ('tag', 'nstar', 'r', 'z', 'se', 'share', 'epm', 'nit', 'curve', 'fsc_dev', 'fsc_cls', 'wall', 'seed_n', 'seed_r', 'seed_z')}
        e['accepted'] = bool(r['z'] >= HALVES_KSIG and r['share'] <= 0.5)
        e['kept'] = bool(e['accepted'] and r['seed_n'] >= 6 and r['seed_z'] >= HALVES_KSIG)   # accepted endpoint whose own seeds reproduce
        if israre is not None:
            o = np.argsort(-post, kind='stable'); top = o[:r['nstar']]; hi = post >= 0.5
            e.update(auc=auc(post, israre), top_P=float(israre[top].mean()), top_R=float(israre[top].sum() / max(israre.sum(), 1)),
                     hi_n=int(hi.sum()), hi_P=float(israre[hi].sum() / max(hi.sum(), 1)), hi_R=float(israre[hi].sum() / max(israre.sum(), 1)),
                     auc_llr=auc(r['llr'], israre) if r['llr'].any() else None)
        ends.append(e)
        devs.append(r['dev'])
        full = np.full(n0, -1.); full[rows] = post
        np.savetxt(os.path.join(a.out, 'site_post_%s.txt' % r['tag']), full, fmt='%8.5f')
        np.save(os.path.join(a.out, 'llr_%s.npy' % r['tag']), r['llr'].astype(np.float32))
    m = len(devs)
    geo = np.zeros((m, m)); gee = np.zeros(m); goo = np.zeros(m)
    for i in range(m):
        gee[i] = shell_corr_dot(devs[i][0], devs[i][0], G['kb'] + RB_MARGIN, G['khi'])[1]
        goo[i] = shell_corr_dot(devs[i][1], devs[i][1], G['kb'] + RB_MARGIN, G['khi'])[1]
        for j in range(m):
            geo[i, j] = shell_corr_dot(devs[i][0], devs[j][1], G['kb'] + RB_MARGIN, G['khi'])[1]
    summary['cross'] = (geo / np.sqrt(np.outer(gee, goo)).clip(1e-300)).tolist()
    summary['endpoints'] = ends
    summary['wall'] = time.time() - t0
    json.dump(summary, open(os.path.join(a.out, 'summary.json'), 'w'), indent=1)
    log('summary: bar r %.3f (n_eff %.1f)' % (summary['bar_r'], G['neff']))
    for e in ends:
        s = '%-6s N* %6d share %.3f r %.3f z %5.2f %s | seeds %d z %.2f %s' % (e['tag'], e['nstar'], e['share'], e['r'], e['z'], 'accepted' if e['accepted'] else 'no',
                                                                          e['seed_n'], e['seed_z'], 'KEPT POLE' if e['kept'] else '')
        if 'auc' in e:
            s += ' | AUC %.3f top-N* P %.3f R %.3f | >=0.5: n %d P %.3f R %.3f' % (e['auc'], e['top_P'], e['top_R'], e['hi_n'], e['hi_P'], e['hi_R'])
        log(s)
    log('done in %.0f s' % summary['wall'])


if __name__ == '__main__':
    main()
