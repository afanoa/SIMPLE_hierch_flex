#!/usr/bin/env python3
"""flexsub.py: remove each particle's within-state continuous deformation (the node's own flex model, at the node's frozen
stage-B poses) from its image and write corrected stacks.
Model (latent_states_v4): y_i = a_i beta CTF_i P_i(V0 + W z_i) + noise over the band to --lp, site_em's per-shell scale beta and
noise, V0 = flex's mean, W = flex's modes cut to the flex band (--lpf), orthonormalised. Slices are 2x-padded trilinear with the
gridding correction (not the nearest-neighbour cells of the prototypes): the subtraction needs accurate projections.
Class Gaussians (tied Omega) of the stage-B labels are fitted on one particle half and evaluated on the other:
  dz_i = E[z_i | y_i, c_i] - xi_{c_i}
Arms: full = W dz_i; perp = W dz_i with the between-state discriminant directions projected out (Omega^-1 metric), so the
subtraction leaves the particle's latent state evidence unchanged.
The subtracted image is alpha_i beta CTF_i P_i(W dz_i), alpha_i refitted per particle against the image with the full latent
model, padded to the box and shifted back to the raw frame. Checks: the band residual against the particle's state model before /
after (in sample) and on checkerboard-held-out cells (dz from the other parity), and an exact re-read of written particles.
usage: flexsub.py --node NODE --dump DUMP --out OUT --stkroot NEWROOT [--arms full,perp --nproc 8 --truth qtoy:TRUTH]"""
import argparse, glob, json, os, sys, time
import numpy as np

LAT = '/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent'
PP = os.environ.get('PPINFO', '/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo')
sys.path.insert(0, LAT)
sys.path.insert(0, PP)
import latent_states_v4 as LS
from ppinfo_site import read_dump, read_stk, euler2m, prep_crop, ctf_eval
from mrcio import read_mrc, fcrop3, stack

G = LS.G
log = LS.log
CH = 400


# ---------------------------------------------------------------- setup (LS.setup, keeping the unweighted planes and the CTF)

def setup(a):
    nd = read_dump(os.path.join(a.dump, 'n_ptcl3D.txt'))
    bd = read_dump(os.path.join(a.dump, 'b_ptcl3D.txt'))
    recs, box, smpd, kv, cs, fraca = read_stk(os.path.join(a.dump, 'n_stk.txt'))
    n0 = len(nd['e1'])
    assert len(bd['e1']) == n0
    didx, prow = LS.node_select(n0, a.node)
    ok = (nd['eo'][didx] == 0) | (nd['eo'][didx] == 1)
    didx, prow = didx[ok], prow[ok]
    if a.maxn and a.maxn < len(didx):
        s = np.sort(np.random.default_rng(3).choice(len(didx), a.maxn, replace=False))
        didx, prow = didx[s], prow[s]
    n = len(didx)
    vol1, mskdiam = LS.flex_keys(a.node)
    box_A = box * smpd
    kb = int(round(box_A / a.lp))
    B = min(box, max(26, 2 * kb + 4)); B += B % 2
    nyq, nh, khi = B // 2, B // 2 + 1, B // 2 - 1
    kb = min(kb, khi)
    kf = int(round(box_A / a.lpf))
    log('node %s: %d particles (dump %d rows), box %d smpd %.3f, mskdiam %.1f -> box_crop %d, band rings 2..%d (%.1f A), flex band '
        '<= ring %d (%.1f A), mean %s' % (a.node, n, n0, box, smpd, mskdiam, B, kb, box_A / kb, kf, box_A / kf, vol1))
    e = np.c_[bd['e1'], bd['e2'], bd['e3'], bd['x'], bd['y']][didx]
    R = euler2m(e[:, :3])
    indstk = nd['indstk'][didx].astype(int)
    stkind = np.maximum(nd['stkind'][didx].astype(int), 1)
    crop = prep_crop(recs, box, stkind, indstk, mskdiam / 2. / smpd, B, a.cache)
    kk, hh = np.meshgrid(np.r_[0:B // 2, -B // 2:0], np.arange(B // 2), indexing='ij')
    r2 = hh * hh + kk * kk
    sel = (r2 >= 3) & (r2 <= kb * (kb + 1))
    h, k = hh[sel].astype(np.int32), kk[sel].astype(np.int32)
    krow = np.where(k >= 0, k, k + B)
    sw = np.where(h == 0, 1., 2.)
    shell = np.rint(np.sqrt(h * h + k * k)).astype(int) - 1
    nsh = shell.max() + 1
    fb = np.rint(np.sqrt(h * h + k * k)) <= kf
    X = crop[:, krow, h]
    del crop
    ph = np.exp(-2j * np.pi * (h[None, :] * e[:, 3:4] + k[None, :] * e[:, 4:5]) / box).astype(np.complex64)
    X *= ph
    del ph
    ctf = ctf_eval(h, k, nd['dfx'][didx], nd['dfy'][didx], nd['angast'][didx], nd['phshift'][didx], kv, cs, fraca, box_A)
    f2 = np.zeros((n, nsh))
    for s in range(nsh):
        m = shell == s
        f2[:, s] = (sw[m][None, :] * (X[:, m].real ** 2 + X[:, m].imag ** 2)).sum(1)
    l = np.arange(B); ls = np.where(l < nyq, l, l - B); hc = np.arange(nh)
    cell_shell = np.minimum(np.rint(np.sqrt(ls[:, None, None] ** 2 + ls[None, :, None] ** 2 + hc[None, None, :] ** 2)), nyq).astype(np.int32).ravel()
    cell_fw = np.where((hc > 0) & (hc < nh - 1), 2., 1.)[None, None, :].repeat(B, 0).repeat(B, 1).ravel()
    G.update(dict(X=X, P=X * ctf, C=ctf * ctf, ctf=ctf, f2=f2, R=R, h=h, k=k, fb=fb, e=e, sw=sw, shell=shell, nsh=nsh, NC=B * B * nh,
                  B=B, nh=nh, nyq=nyq, kb=kb, kf=kf, box=box, smpd=smpd, box_A=box_A, mskdiam=mskdiam, cell_shell=cell_shell,
                  cell_fw=cell_fw, eo=nd['eo'][didx].astype(int), n=n, n0=n0, didx=didx, prow=prow, vol1=vol1, recs=recs,
                  indstk=indstk, stkind=stkind, cn=np.ones(n), lab=bd['state'][didx].astype(int)))


# ---------------------------------------------------------------- 2x-padded trilinear slices

class PadVol:
    """a real volume at box_crop B (origin B//2), gridding-corrected for trilinear interpolation, zero-padded to 2B, full FFT."""

    def __init__(self, v, cut=False):
        B = G['B']; Bp = 2 * B
        vp = np.zeros((Bp, Bp, Bp))
        o = Bp // 2 - B // 2
        vp[o:o + B, o:o + B, o:o + B] = v
        s = np.sinc((np.arange(Bp) - Bp // 2) / Bp) ** 2
        vp /= s[:, None, None] * s[None, :, None] * s[None, None, :]
        self.F = np.fft.fftn(np.fft.ifftshift(vp)).ravel()
        self.cut = cut      # zero the plane cells beyond the flex band


TRI = {}


def tri(sl):
    key = (sl.start, sl.stop)
    if TRI.get('key') == key:
        return TRI['idx'], TRI['w']
    Bp = 2 * G['B']
    R = G['R'][sl]
    h, k = G['h'].astype(float), G['k'].astype(float)
    p = 2. * (h[None, :, None] * R[:, None, 0, :] + k[None, :, None] * R[:, None, 1, :])
    f = np.floor(p)
    t = p - f
    f = f.astype(np.int64)
    idx = np.empty((8,) + p.shape[:2], np.int64)
    w = np.empty((8,) + p.shape[:2])
    c = 0
    for dz in (0, 1):
        for dy in (0, 1):
            for dx in (0, 1):
                idx[c] = (((f[..., 2] + dz) % Bp) * Bp + (f[..., 1] + dy) % Bp) * Bp + (f[..., 0] + dx) % Bp
                w[c] = (t[..., 0] if dx else 1. - t[..., 0]) * (t[..., 1] if dy else 1. - t[..., 1]) * (t[..., 2] if dz else 1. - t[..., 2])
                c += 1
    TRI.update(key=key, idx=idx, w=w)
    return idx, w


def slices_tri(V, sl):
    idx, w = tri(sl)
    s = (V.F[idx] * w).sum(0)
    return s * G['fb'][None] if V.cut else s


LS.slices = slices_tri


def real_from_ft(V):
    """an LS.vol_ft half volume (raveled, trailing zero) back to the real volume at box_crop."""
    B, nh = G['B'], G['nh']
    return np.fft.fftshift(np.fft.irfftn(V[:G['NC']].reshape(B, B, nh), s=(B, B, B)))


# ---------------------------------------------------------------- statistics

def stats(V0, Ws, beta, sig2, amp, par):
    """latent sufficient statistics over all band cells and over the checkerboard cells of parity par only."""
    n, q = G['n'], len(Ws)
    w1 = G['sw'] * (beta / sig2)[G['shell']]
    w2 = G['sw'] * (beta ** 2 / sig2)[G['shell']]
    out = []
    for msk in (np.ones(len(w1), bool), par):
        out.append((np.zeros((n, q)), np.zeros((n, q, q)), msk))
    t0 = time.time()
    for c0 in range(0, n, CH):
        sl = slice(c0, min(n, c0 + CH))
        a = amp[sl][:, None]
        S = np.stack([slices_tri(V, sl) for V in Ws], 1)
        s0 = slices_tri(V0, sl)
        for b, Gm, msk in out:
            T = G['C'][sl] * (w2 * msk)[None] * a ** 2
            Rr = G['P'][sl] * ((w1 * msk)[None] * a) - s0 * T
            b[sl] = np.einsum('iqm,im->iq', S.real, Rr.real) + np.einsum('iqm,im->iq', S.imag, Rr.imag)
            Gm[sl] = np.matmul(S.real * T[:, None], S.real.transpose(0, 2, 1)) + np.matmul(S.imag * T[:, None], S.imag.transpose(0, 2, 1))
    log('sufficient statistics: %d particles x %d modes (all cells + parity half) in %.0f s' % (n, q, time.time() - t0))
    return out[0][:2], out[1][:2]


def posterior_dz(b, Gm, f, lab):
    """E[z | y, c] - xi_c for each particle at its label c, under the class Gaussians f."""
    Oi = np.linalg.inv(f['Om'])
    H = Gm + Oi[None]
    xi = f['xi'][lab]
    m = np.linalg.solve(H, (b + xi @ Oi.T)[:, :, None])[:, :, 0]
    return m - xi, m


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--node'); ap.add_argument('--dump'); ap.add_argument('--out'); ap.add_argument('--stkroot')
    ap.add_argument('--cache', default=''); ap.add_argument('--lp', type=float, default=10.); ap.add_argument('--lpf', type=float, default=15.)
    ap.add_argument('--arms', default='full,perp'); ap.add_argument('--nproc', type=int, default=8); ap.add_argument('--truth', default='none')
    ap.add_argument('--maxn', type=int, default=0)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    setup(a)
    n, B = G['n'], G['B']
    lab0 = G['lab']
    states = np.unique(lab0)
    assert states.min() >= 1, 'unassigned particles in the node'
    lab = np.searchsorted(states, lab0)
    K = len(states)
    log('stage-B labels: %s' % ' '.join('%d:%d' % (s, (lab0 == s).sum()) for s in states))
    # volumes: the mean in the full band, the modes cut to the flex band and orthonormalised there
    V0r = read_mrc(G['vol1'])[0]
    V0 = PadVol(fcrop3(V0r, B) if V0r.shape[0] != B else V0r)
    pcs = sorted(glob.glob(os.path.join(a.node, 'fflex', 'flex_pca_pc[0-9][0-9][0-9].mrc')))
    Vf = []
    for f in pcs:
        V = LS.vol_ft(read_mrc(f)[0])
        Vf.append(np.r_[np.where(G['cell_shell'] <= G['kf'], V[:G['NC']], 0.), 0.])
    W, _ = LS.orthonormal(Vf, ['pc%d' % (j + 1) for j in range(len(Vf))])
    Ws = [PadVol(real_from_ft(Wl), cut=True) for Wl in W]
    q = len(Ws)
    beta, sig2, amp = LS.scale_noise(V0)
    par = ((G['h'] + G['k']) % 2) == 0
    (b, Gm), (bA, GA) = stats(V0, Ws, beta, sig2, amp, par)
    G.update(b=b, Gm=Gm, mu1=np.zeros(q))
    f1 = LS.fit_one(np.arange(n), np.eye(q) / max(np.median(np.trace(Gm, axis1=1, axis2=2)) / q, 1e-30) * 10.)
    Om1, G['mu1'] = f1['Om'], f1['xi'][0]
    log('K=1: Omega eigenvalues %s' % np.array2string(np.sort(np.linalg.eigvalsh(Om1))[::-1], precision=3))
    eo = G['eo']
    ia, ib = np.where(eo == 0)[0], np.where(eo == 1)[0]
    dz = np.zeros((n, q)); mz = np.zeros((n, q)); dzA = np.zeros((n, q)); dzp = np.zeros((n, q))
    fits = {}
    for fs, es, nm in ((ia, ib, 'even->odd'), (ib, ia, 'odd->even')):
        f = LS.xd_fit(b[fs], Gm[fs], LS.onehot(lab[fs], K), Om1, G['mu1'], fixed=True)
        fits[nm] = f
        dz[es], mz[es] = posterior_dz(b[es], Gm[es], f, lab[es])
        dzA[es] = posterior_dz(bA[es], GA[es], f, lab[es])[0]
        # perp: project out the between-state discriminant directions in the Omega^-1 metric
        Oi = np.linalg.inv(f['Om'])
        D = (f['xi'][1:] - f['xi'][0]).T
        Pm = D @ np.linalg.solve(D.T @ Oi @ D, D.T @ Oi)
        dzp[es] = dz[es] - dz[es] @ Pm.T
        ev = np.linalg.eigvalsh(f['Om'])
        sep = [float(np.sqrt((f['xi'][i] - f['xi'][j]) @ Oi @ (f['xi'][i] - f['xi'][j]))) for i in range(K) for j in range(i + 1, K)]
        log('%s class Gaussians: pi %s, Omega trace %.3g (top eig %.3g), between-state Mahalanobis %s'
            % (nm, np.round(f['pi'], 3), ev.sum(), ev.max(), ' '.join('%.2f' % s for s in sep)))
    # within-state spread actually removed, in the latent metric
    Oi1 = np.linalg.inv(Om1)
    md = lambda x: np.einsum('iq,qp,ip->i', x, Oi1, x)
    log('latent: |dz|^2_Om1 median %.3f, |dz_perp|^2 %.3f, |m - mu1|^2 %.3f' % (np.median(md(dz)), np.median(md(dzp)), np.median(md(mz - G['mu1']))))
    # pass 2: amplitude refit, subtracted components, residual checks
    arms = a.arms.split(',')
    dzs = {'full': dz, 'perp': dzp}
    Dc = {arm: np.zeros((n, len(G['h'])), np.complex64) for arm in arms}
    w = G['sw'] / sig2[G['shell']]
    fbm = G['fb'] & (G['shell'] >= 1)
    alpha = np.zeros(n)
    res = {k_: np.zeros(n) for k_ in ('pow_state', 'drop_in', 'drop_held', 'drop_held_in', 'pow_D', 'pow_X', 'drop_perp')}
    for c0 in range(0, n, CH):
        sl = slice(c0, min(n, c0 + CH))
        nc = sl.stop - sl.start
        S = np.stack([slices_tri(V, sl) for V in Ws], 1)
        s0 = slices_tri(V0, sl)
        bc = beta[G['shell']][None] * G['ctf'][sl]
        X = G['X'][sl].astype(np.complex128)
        xi = np.stack([fits['even->odd' if eo[i] == 1 else 'odd->even']['xi'][lab[i]] for i in range(sl.start, sl.stop)])
        Mf = bc * (s0 + np.einsum('iq,iqm->im', mz[sl], S))
        al = (w[None] * (Mf.conj() * X).real)[:, G['shell'] >= 1].sum(1) / np.maximum((w[None] * np.abs(Mf) ** 2)[:, G['shell'] >= 1].sum(1), 1e-30)
        alpha[sl] = al
        Ms = al[:, None] * bc * (s0 + np.einsum('iq,iqm->im', xi, S))
        Rs = X - Ms
        for arm in arms:
            Dc[arm][sl] = (al[:, None] * bc * np.einsum('iq,iqm->im', dzs[arm][sl], S)).astype(np.complex64)
        D = Dc['full'][sl] if 'full' in Dc else al[:, None] * bc * np.einsum('iq,iqm->im', dz[sl], S)
        DA = al[:, None] * bc * np.einsum('iq,iqm->im', dzA[sl], S)
        wr = w[None] * fbm[None]
        res['pow_state'][sl] = (wr * np.abs(Rs) ** 2).sum(1)
        res['drop_in'][sl] = (wr * (np.abs(Rs) ** 2 - np.abs(Rs - D) ** 2)).sum(1)
        wB = wr * (~par)[None]
        res['drop_held'][sl] = (wB * (np.abs(Rs) ** 2 - np.abs(Rs - DA) ** 2)).sum(1)
        res['drop_held_in'][sl] = (wB * (np.abs(Rs) ** 2 - np.abs(Rs - D) ** 2)).sum(1)
        res['pow_D'][sl] = (wr * np.abs(D) ** 2).sum(1)
        res['pow_X'][sl] = (wr * np.abs(X) ** 2).sum(1)
        if 'perp' in Dc:
            res['drop_perp'][sl] = (wr * (np.abs(Rs) ** 2 - np.abs(Rs - Dc['perp'][sl]) ** 2)).sum(1)
    log('amplitude: refit alpha / scale-noise a: median %.3f (5-95 %% %.3f %.3f); corr(alpha, a) %.3f'
        % (np.median(alpha / amp), *np.percentile(alpha / amp, [5, 95]), np.corrcoef(alpha, amp)[0, 1]))
    nrm = res['pow_state']
    log('flex band (rings 2..%d), noise-whitened, per particle (sum, then share of the residual against the state model):' % G['kf'])
    log('   |X|^2 %.1f, |X - state model|^2 %.1f, |D|^2 %.2f (%.2f %% of the residual)' % (np.median(res['pow_X']), np.median(nrm), np.median(res['pow_D']), 100 * np.median(res['pow_D'] / nrm)))
    log('   residual drop in sample: mean %.3f nats (%.3f %% of the residual); particles with a drop %.3f'
        % (0.5 * res['drop_in'].mean(), 100 * res['drop_in'].sum() / nrm.sum(), (res['drop_in'] > 0).mean()))
    log('   held-out cells (dz from the other parity): drop %.3f nats/particle vs %.3f in sample on the same cells (ratio %.2f); particles with a drop %.3f'
        % (0.5 * res['drop_held'].mean(), 0.5 * res['drop_held_in'].mean(), res['drop_held'].sum() / max(res['drop_held_in'].sum(), 1e-30), (res['drop_held'] > 0).mean()))
    if 'perp' in Dc:
        log('   perp arm in-sample drop %.3f nats/particle' % (0.5 * res['drop_perp'].mean()))
    # truth diagnostics: does dz track the arm continuum within states?
    tr = {}
    if a.truth.startswith('qtoy:'):
        tf = a.truth.split(':', 1)[1]
        hdr = open(tf).readline().lstrip('#').split()
        t = np.loadtxt(tf, ndmin=2)[G['prow'] - 1]
        if 's1' in hdr:
            s1 = t[:, hdr.index('s1')]
            arm = np.isfinite(s1) & (t[:, hdr.index('junk')] == 0) & (t[:, hdr.index('C_closed')] == 0) & (t[:, hdr.index('holo')] == 0)
            def r2(x, y):
                A = np.c_[x, np.ones(len(x))]
                c = np.linalg.lstsq(A, y, rcond=None)[0]
                return 1. - ((y - A @ c) ** 2).sum() / max(((y - y.mean()) ** 2).sum(), 1e-30)
            tr['r2_s1_m_all'] = r2(mz[arm], s1[arm])
            for s in range(K):
                mm = arm & (lab == s)
                if mm.sum() > 200:
                    tr['state%d' % (s + 1)] = dict(n_arm=int(mm.sum()), s1_sd=float(s1[mm].std()), r2_s1_dz=r2(dz[mm], s1[mm]),
                                                   r2_s1_dzperp=r2(dzp[mm], s1[mm]))
            log('truth: R^2 of the arm s1 on the latent posterior mean (all arm particles) %.3f; within states (s1 sd, R^2 on dz, on dz_perp): %s'
                % (tr['r2_s1_m_all'], '; '.join('%s n %d sd %.2f %.3f %.3f' % (k_, v['n_arm'], v['s1_sd'], v['r2_s1_dz'], v['r2_s1_dzperp'])
                                                 for k_, v in tr.items() if k_.startswith('state'))))
    np.savez(os.path.join(a.out, 'flexsub.npz'), prow=G['prow'], lab=lab0, dz=dz, dzp=dzp, dzA=dzA, mz=mz, alpha=alpha, amp=amp, beta=beta,
             sig2=sig2, **res)
    # write the corrected stacks: every image copied, the node's particles corrected
    recs, box = G['recs'], G['box']
    paths = [r[0] for r in recs]
    old_root = os.path.commonpath([os.path.dirname(p) for p in paths])
    krowf = np.where(G['k'] >= 0, G['k'], G['k'] + box)
    summary = dict(node=a.node, n=n, K=K, q=q, old_root=old_root, stacks={}, tr=tr,
                   drop_in_nats=float(0.5 * res['drop_in'].mean()), drop_held_nats=float(0.5 * res['drop_held'].mean()),
                   drop_held_in_nats=float(0.5 * res['drop_held_in'].mean()), powD_share=float(np.median(res['pow_D'] / nrm)))
    for arm in arms:
        root = os.path.join(a.stkroot, arm)
        for si, p in enumerate(paths):
            dst = os.path.join(root, os.path.relpath(p, old_root))
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            src = stack(p)
            nx_, ny_, nz_ = src.shape[2], src.shape[1], src.shape[0]
            off = src.offset
            with open(p, 'rb') as fi, open(dst, 'wb') as fo:
                fo.write(fi.read(off))
            out = np.memmap(dst, np.float32, 'r+', offset=off, shape=(nz_, ny_, nx_))
            for i0 in range(0, nz_, 5000):
                out[i0:i0 + 5000] = src[i0:i0 + 5000]
            ii = np.where(G['stkind'] == si + 1)[0]
            ph = None
            for c0 in range(0, len(ii), 2000):
                sel = ii[c0:c0 + 2000]
                e = G['e'][sel]
                Dr = Dc[arm][sel].astype(np.complex128) * np.exp(2j * np.pi * (G['h'][None, :] * e[:, 3:4] + G['k'][None, :] * e[:, 4:5]) / box)
                Ff = np.zeros((len(sel), box, box // 2 + 1), np.complex128)
                Ff[:, krowf, G['h']] = Dr
                img = np.fft.fftshift(np.fft.irfft2(Ff, s=(box, box)), axes=(1, 2))
                j = G['indstk'][sel] - 1
                out[j] = (np.asarray(out[j], np.float64) - img).astype(np.float32)
            out.flush()
            del out
            # exact check on a few particles: re-read, re-plane, compare with X - D
            rng = np.random.default_rng(5)
            chk = rng.choice(ii, min(300, len(ii)), replace=False)
            st2 = stack(dst)
            img = np.asarray(st2[G['indstk'][chk] - 1], np.float64)
            F = np.fft.rfft2(np.fft.ifftshift(img, axes=(1, 2)))[:, krowf, G['h']]
            e = G['e'][chk]
            Xc = F * np.exp(-2j * np.pi * (G['h'][None, :] * e[:, 3:4] + G['k'][None, :] * e[:, 4:5]) / box)
            want = G['X'][chk].astype(np.complex128) - Dc[arm][chk]
            err = np.abs(Xc - want).sum() / np.abs(want).sum()
            unchanged = np.abs(Xc[:, ~G['fb']] - G['X'][chk][:, ~G['fb']]).sum() / np.abs(G['X'][chk][:, ~G['fb']]).sum()
            log('%s: wrote %s (%d images, %d corrected); re-read check: |planes - (X - D)| / |X - D| %.2e, outside the flex band %.2e'
                % (arm, dst, nz_, len(ii), err, unchanged))
            summary['stacks'][arm] = dict(root=root, path=dst, check=float(err), check_out=float(unchanged))
    summary['wall'] = time.time() - t0
    json.dump(summary, open(os.path.join(a.out, 'summary.json'), 'w'), indent=1, default=float)
    log('done in %.0f s' % summary['wall'])


if __name__ == '__main__':
    main()
