#!/usr/bin/env python3
"""Kernel-latent EM (generative-topographic-mapping style) for cryo-EM particles: a 2-D latent grid of n1 x n2 cells,
each cell a FREE volume, the manifold made smooth by a radial-basis kernel on the grid (the M-step reconstructs every
cell from kernel-smoothed responsibilities), plus nfree off-grid components for junk. Nonlinear in volume space (large
motions are curves through the grid) and agnostic to composition (a cell may lack a domain). flex_pca is the special
case with a linear map from the latent to the volume.
  E-step: responsibilities from the band-limited image likelihood against the cross-half cell volumes (particles of one
          half are scored against volumes built from the other half)
  M-step: soft Wiener backprojection of every particle into every component on a coarse Fourier grid (nearest voxel,
          Friedel mates), weights = responsibilities convolved with the grid kernel; shrinkage toward the consensus
  schedule: bands 20 -> 14 -> 10 A, kernel width annealed; init from two contrast-free components of a flex latent
  readouts (truth): posterior-mean (z1, z2) vs the conformer axis, 15-NN accuracy, pairwise d', closed AUC, cell purity,
  junk capture by the free components, per-conformer mean position
usage: klem.py STK TABLE GT OUT --cons SIMPLE_consensus.mrc --flexdir DIR --states .. --axis .. [--n1 7 --n2 7 --nfree 2]
       [--lps "20 20 20 14 14 14 10 10 10"] [--sigmas "1.5 1.5 1.2 1.0 1.0 0.8 0.7 0.7 0.7"] [--nsub 0] [--nproc 8] [--init pca|random]"""
import sys, os, argparse, time, json
import numpy as np
os.environ.setdefault('TOY_SMPD', '1.464')
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/routes/pvbag/subclass_pool'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from view_tangent_toy import euler2m, rotmat, ctf_eval, Band, fft2c, SMPD
from latent_windows import read_gt, read_table, slice_batch, Data, calibrate, shift_phase, parse_axis
from model_class import read_mrc, write_mrc, Shells, fftc
from scipy.stats import spearmanr
from scipy.spatial import cKDTree
import multiprocessing as mp

G = {}   # shared state for forked workers


def _load_chunk(args):
    s0, s1 = args; D, tab, ss, fs, box, Bd = G['D'], G['tab'], G['ss'], G['fs'], G['box'], G['Bd']; idx = G['sel'][s0:s1]
    F, C, R = D.chunk(idx); F = F * np.conj(shift_phase(Bd, tab[idx, 3], tab[idx, 4], box, ss, fs))   # shift folded into the data
    return F.astype(np.complex64), C.astype(np.float32), R.astype(np.float32)


def _estep_chunk(args):
    s0, s1 = args; F, C, R, eo = G['F'][s0:s1], G['C'][s0:s1], G['R'][s0:s1], G['eo'][s0:s1]; REF = G['REF']; Bd = G['Bd']; conv = G['conv']; SC = G['SC']; sig2 = G['sig2']
    K = len(REF[0]); ll = np.zeros((s1 - s0, K)); NB = G.get('NB'); ZP = None
    if NB is not None: q = len(NB[0]); ZP = np.zeros((s1 - s0, K, q), np.float32); ws = (Bd.wq / sig2)[None, :]
    for h in (0, 1):
        m = eo == h
        if not m.any(): continue
        Fm, Cm, Rm = F[m].astype(complex), C[m].astype(float), R[m].astype(float)
        if NB is not None:   # node basis of the other half projected at the images' poses: Q (m, q, nband); M = I + 2 Re(Q^H W Q)
            Q = np.stack([Cm * slice_batch(re, im, Rm, Bd.kx, Bd.ky, conv) * SC for (re, im) in NB[1 - h]], 1)
            Mq = np.eye(q)[None] + 2. * np.real(np.einsum('mjb,mlb->mjl', np.conj(Q) * ws[:, None, :], Q)); Mi = np.linalg.inv(Mq); ldet = np.linalg.slogdet(Mq)[1]
        for k in range(K):
            re, im = REF[1 - h][k]; P = Cm * slice_batch(re, im, Rm, Bd.kx, Bd.ky, conv) * SC
            ll[m, k] = -np.sum(Bd.wq * np.abs(Fm - P) ** 2 / sig2, 1)
            if NB is not None:   # evidence marginalised over z ~ N(0, I): + 1/2 beta' M^-1 beta - 1/2 log det M; posterior mean z = M^-1 beta
                beta = 2. * np.real(np.einsum('mjb,mb->mj', np.conj(Q) * ws[:, None, :], Fm - P)); z = np.einsum('mjl,ml->mj', Mi, beta); ZP[m, k] = z
                if G.get('NB_MARG', 1) and not G.get('NB_INIT', False): ll[m, k] += 0.5 * np.sum(beta * z, 1) - 0.5 * ldet
    return ll if ZP is None else (ll, ZP)


def _class_frame(R, x, y, e3, x2, y2, combo, Bd, box, ss, fs):
    """slice rotations and shift phases that render a volume as the members' CLASS-FRAME images: image rotated by the
    2D in-plane angle (sign sr), its 3D-pose shift rotated along (rot3: 0 none, +-1), plus the 2D shift (sign s2)"""
    sr, s2, rot3 = combo; Rz = np.array([rotmat(sr * t, 3) for t in e3]); Rn = np.einsum('mij,mjk->mik', Rz, R)
    if rot3 == 0: xs, ys = x, y
    else:
        c, s_ = np.cos(np.radians(rot3 * e3)), np.sin(np.radians(rot3 * e3)); xs, ys = c * x - s_ * y, s_ * x + c * y
    return Rn, shift_phase(Bd, xs + s2 * x2, ys + s2 * y2, box, ss, fs)


def _estep_bag_chunk(bags):
    """sub-class E-step for a list of bags: the bag's half-h average is scored against pose-mean renders of the other
    half's component volumes over a sample of the bag's half-h MEMBERS' poses (all members, not only the EM subset),
    as residual (average minus consensus render, per-shell scaled) vs deviation (component minus consensus render),
    noise per shell from the even/odd averages. Returns {bag: [log-lik half 0, log-lik half 1]}."""
    REF, GC, Bd, conv, box = G['REF'], G['GC'], G['Bd'], G['conv'], G['box']; K = len(REF[0]); out = {}
    AE, AO, ASIG, arow, BP = G['AE'], G['AO'], G['ASIG'], G['AROW'], G['BAGPOSE']; shell = Bd.shell; nsh = shell.max() + 1
    for b in bags:
        row = arow.get(int(b), -1)
        if row < 0 or int(b) not in BP: continue
        res = []
        for h in (0, 1):
            Rm, xs, ys = BP[int(b)][h]
            if len(Rm) == 0: res.append(None); continue
            phm = shift_phase(Bd, xs, ys, box, G['ss'], G['fs'])
            Pc = np.mean(slice_batch(GC[1 - h][0], GC[1 - h][1], Rm, Bd.kx, Bd.ky, conv) * phm, 0)
            A = (AE if h == 0 else AO)[row].astype(complex); sig2 = np.maximum(ASIG[row], 1e-12)[shell]
            num = np.bincount(shell, Bd.wq * np.real(np.conj(Pc) * A), minlength=nsh); den = np.bincount(shell, Bd.wq * np.abs(Pc) ** 2, minlength=nsh)
            sc = np.where(den > 0, num / np.maximum(den, 1e-30), 0.); sc = np.convolve(np.pad(sc, 1, mode='edge'), np.ones(3) / 3., mode='valid'); sc = np.maximum(sc, 0.)[shell]
            resid = A - sc * Pc; llb = np.zeros(K); D = np.zeros((K, Bd.n), complex)
            for k in range(K):
                Pk = np.mean(slice_batch(REF[1 - h][k][0], REF[1 - h][k][1], Rm, Bd.kx, Bd.ky, conv) * phm, 0); D[k] = sc * (Pk - Pc)
                llb[k] = -np.sum(Bd.wq * np.abs(resid - D[k]) ** 2 / sig2)
            mode_ = G.get('AVGMODE', 'pool')
            if mode_ == 'unmix':   # average residual = sum_k pi_k x deviation_k, pi >= 0 (whitened, real + imag stacked)
                from scipy.optimize import nnls
                w_ = np.sqrt(Bd.wq / sig2); X = np.concatenate([(D * w_).real, (D * w_).imag], 1).T; y = np.concatenate([(resid * w_).real, (resid * w_).imag])
                pi, _ = nnls(X, y, maxiter=200 * K); pi = pi / pi.sum() if pi.sum() > 0 else np.ones(K) / K
                res.append(np.log(pi + 1e-3))   # log mixture prior (floored: a member may still go where its average has no weight)
            elif mode_ == 'unmixc':   # parent-centred: the within-parent difference of averages (same view, the view part cancels) regressed on the deviations
                pp = G['PAROF'].get(int(b)); PA = G.get('PARAVG', {})
                if pp is None or pp not in PA: res.append(np.zeros(K)); continue
                Ap, per = PA[pp]; Rp, xp, yp = per[h]
                if len(Rp) == 0: res.append(np.zeros(K)); continue
                php = shift_phase(Bd, xp, yp, box, G['ss'], G['fs']); Pcp = np.mean(slice_batch(GC[1 - h][0], GC[1 - h][1], Rp, Bd.kx, Bd.ky, conv) * php, 0); A_p = Ap[h].astype(complex)
                nump = np.bincount(shell, Bd.wq * np.real(np.conj(Pcp) * A_p), minlength=nsh); denp = np.bincount(shell, Bd.wq * np.abs(Pcp) ** 2, minlength=nsh)
                scp_ = np.where(denp > 0, nump / np.maximum(denp, 1e-30), 0.); scp_ = np.convolve(np.pad(scp_, 1, mode='edge'), np.ones(3) / 3., mode='valid'); scp_ = np.maximum(scp_, 0.)[shell]
                dres = resid - (A_p - scp_ * Pcp)   # sub-class residual minus parent residual
                w_ = np.sqrt(Bd.wq / sig2); X = np.concatenate([(D * w_).real, (D * w_).imag], 1).T; y = np.concatenate([(dres * w_).real, (dres * w_).imag])
                tr_ = np.trace(X.T @ X)
                if not np.isfinite(tr_) or tr_ <= 0: res.append(np.zeros(K)); continue
                lam = 1e-2 * tr_ / K; c = np.linalg.lstsq(np.vstack([X, np.sqrt(lam) * np.eye(K)]), np.concatenate([y, np.zeros(K)]), rcond=None)[0]; c = (c - c.mean()) / max(c.std(), 1e-12)
                res.append(c)   # standardised within-parent preference among the components (zero-mean: no absolute, view-carrying part)
            else: res.append(llb)
        out[int(b)] = res
    return out


_FORT = {}


def fort_lib(nproc):
    """the OpenMP kernels (libklem.so beside klem.py or $KLEM_FORT); None -> numpy paths"""
    if 'lib' in _FORT: return _FORT['lib']
    import ctypes, ctypes.util
    path = os.environ.get('KLEM_FORT', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'libklem.so'))
    lib = None
    if os.path.exists(path) and os.environ.get('KLEM_NOFORT', '') == '':
        try:
            lib = ctypes.CDLL(path); dp = ctypes.POINTER(ctypes.c_double); fp = ctypes.POINTER(ctypes.c_float); ip = ctypes.POINTER(ctypes.c_int); ci = ctypes.c_int
            lib.klem_set_threads.argtypes = [ci]; lib.klem_set_threads(int(nproc))
            lib.klem_estep.argtypes = [ci, ci, ci, ci, dp, dp, dp, dp, ip, fp, fp, dp, dp, dp, dp, dp, ci, dp]; lib.klem_estep.restype = None
            lib.klem_mstep.argtypes = [ci, ci, ci, ci, ci, dp, dp, dp, dp, ip, dp, dp, dp, ci, dp, dp, dp]; lib.klem_mstep.restype = None
            print('klem: Fortran kernels from %s (%d threads)' % (path, nproc), flush=True)
        except Exception as e: print('klem: Fortran kernels unavailable (%s) -> numpy' % e, flush=True); lib = None
    _FORT['lib'] = lib; return lib


def _fort_refs():
    REF = G['REF']; K = len(REF[0]); n3 = REF[0][0][0].shape[0]
    re = np.ascontiguousarray(np.array([[REF[h][k][0] for k in range(K)] for h in (0, 1)]), np.float32); im = np.ascontiguousarray(np.array([[REF[h][k][1] for k in range(K)] for h in (0, 1)]), np.float32)
    return re, im, K, n3


def fort_estep(lib, chunks):
    """ll (n, K) over all chunks with the Fortran kernel (OpenMP over particles)"""
    import ctypes
    re, im, K, n3 = _fort_refs(); Bd = G['Bd']; n = G['F'].shape[0]; ll = np.zeros((n, K)); dp = ctypes.POINTER(ctypes.c_double); fp = ctypes.POINTER(ctypes.c_float); ip = ctypes.POINTER(ctypes.c_int)
    kx, ky, wq, sig2, SC = (np.ascontiguousarray(v, np.float64) for v in (Bd.kx, Bd.ky, Bd.wq, G['sig2'], G['SC']))
    for s0, s1 in chunks:
        F = G['F'][s0:s1]; fre = np.ascontiguousarray(F.real, np.float64); fim = np.ascontiguousarray(F.imag, np.float64); C = np.ascontiguousarray(G['C'][s0:s1], np.float64); R = np.ascontiguousarray(G['R'][s0:s1], np.float64); eo = np.ascontiguousarray(G['eo'][s0:s1], np.int32)
        out = np.zeros((s1 - s0, K))
        lib.klem_estep(s1 - s0, K, Bd.n, n3, fre.ctypes.data_as(dp), fim.ctypes.data_as(dp), C.ctypes.data_as(dp), R.ctypes.data_as(dp), eo.ctypes.data_as(ip), re.ctypes.data_as(fp), im.ctypes.data_as(fp),
                       kx.ctypes.data_as(dp), ky.ctypes.data_as(dp), wq.ctypes.data_as(dp), sig2.ctypes.data_as(dp), SC.ctypes.data_as(dp), int(G['conv']), out.ctypes.data_as(dp))
        ll[s0:s1] = out
    return ll


def fort_mstep(lib, chunks, Wfun, K):
    """num (2, K, nv) complex and den (2, K, nv) accumulated over all chunks with the Fortran kernel (OpenMP over components)"""
    import ctypes
    kmax, Bd = G['kmax'], G['Bd']; n3 = 2 * kmax + 1; dp = ctypes.POINTER(ctypes.c_double); ip = ctypes.POINTER(ctypes.c_int)
    numre = np.zeros((2, K, n3, n3, n3)); numim = np.zeros((2, K, n3, n3, n3)); den = np.zeros((2, K, n3, n3, n3)); kx, ky = (np.ascontiguousarray(v, np.float64) for v in (Bd.kx, Bd.ky))
    for s0, s1 in chunks:
        F = G['F'][s0:s1]; fre = np.ascontiguousarray(F.real, np.float64); fim = np.ascontiguousarray(F.imag, np.float64); C = np.ascontiguousarray(G['C'][s0:s1], np.float64); R = np.ascontiguousarray(G['R'][s0:s1], np.float64); eo = np.ascontiguousarray(G['eo'][s0:s1], np.int32)
        W = np.ascontiguousarray(Wfun(s0, s1), np.float64); assert W.shape == (s1 - s0, K)
        lib.klem_mstep(s1 - s0, K, Bd.n, n3, kmax, fre.ctypes.data_as(dp), fim.ctypes.data_as(dp), C.ctypes.data_as(dp), R.ctypes.data_as(dp), eo.ctypes.data_as(ip), W.ctypes.data_as(dp), kx.ctypes.data_as(dp), ky.ctypes.data_as(dp), int(G['conv']),
                       numre.ctypes.data_as(dp), numim.ctypes.data_as(dp), den.ctypes.data_as(dp))
    return (numre + 1j * numim).reshape(2, K, -1), den.reshape(2, K, -1)


def _consensus_stats(chunks):
    """per particle: <P, F>, <P, P>, |F - P|^2 / sig2 summed over the band against the OTHER half's consensus (SC applied)"""
    Bd, conv, SC, sig2 = G['Bd'], G['conv'], G['SC'], G['sig2']; n = G['F'].shape[0]; num = np.zeros(n); den = np.zeros(n); res = np.zeros(n)
    for s0, s1 in chunks:
        F, C, R, eo = G['F'][s0:s1].astype(complex), G['C'][s0:s1].astype(float), G['R'][s0:s1].astype(float), G['eo'][s0:s1]
        for h in (0, 1):
            m = np.where(eo == h)[0]
            if len(m) == 0: continue
            re, im = G['GC'][1 - h]; P = C[m] * slice_batch(re, im, R[m], Bd.kx, Bd.ky, conv) * SC
            num[s0 + m] = np.sum(Bd.wq * np.real(np.conj(P) * F[m]), 1); den[s0 + m] = np.sum(Bd.wq * np.abs(P) ** 2, 1); res[s0 + m] = np.sum(Bd.wq * np.abs(F[m] - P) ** 2 / sig2, 1)
    return num, den, res


def _mstep_chunk(args):
    s0, s1, W = args; F, C, R, eo = G['F'][s0:s1].astype(complex), G['C'][s0:s1].astype(float), G['R'][s0:s1].astype(float), G['eo'][s0:s1]
    kmax, Bd, conv = G['kmax'], G['Bd'], G['conv']; n3 = 2 * kmax + 1; nv = n3 ** 3; K = W.shape[1]
    hk = np.stack([Bd.kx, Bd.ky, np.zeros_like(Bd.kx)], 1); M = R if conv == 0 else np.transpose(R, (0, 2, 1)); loc = np.einsum('kj,mji->mki', hk, M)
    ir = np.rint(loc).astype(int); idx = ir + kmax; idxm = -ir + kmax
    ok = np.all((idx >= 0) & (idx < n3), axis=2); okm = np.all((idxm >= 0) & (idxm < n3), axis=2)
    flat = (idx[..., 2] * n3 + idx[..., 1]) * n3 + idx[..., 0]; flatm = (idxm[..., 2] * n3 + idxm[..., 1]) * n3 + idxm[..., 0]
    CF = C * F; C2 = C * C; num = np.zeros((2, K, nv), complex); den = np.zeros((2, K, nv))
    for h in (0, 1):
        m = np.where(eo == h)[0]
        if len(m) == 0: continue
        fl, flm, o, om = flat[m], flatm[m], ok[m], okm[m]; cf, c2 = CF[m], C2[m]
        for k in range(K):
            w = W[m, k]
            if w.max() < 1e-4: continue
            wk = w[:, None] * cf; wd = w[:, None] * c2
            a = np.bincount(fl[o], weights=wk[o].real, minlength=nv) + 1j * np.bincount(fl[o], weights=wk[o].imag, minlength=nv)
            b = np.bincount(flm[om], weights=wk[om].real, minlength=nv) - 1j * np.bincount(flm[om], weights=wk[om].imag, minlength=nv)
            num[h, k] += a + b; den[h, k] += np.bincount(fl[o], weights=wd[o], minlength=nv) + np.bincount(flm[om], weights=wd[om], minlength=nv)
    return num, den


def m2euler_simple(R):
    """inverse of euler2m for R = r3(e3) r2(e2) r1(e1) (SIMPLE), returns (e1, e2, e3) in degrees"""
    e2 = np.degrees(np.arccos(np.clip(R[2, 2], -1, 1)))
    if abs(R[2, 2]) < 1 - 1e-8:
        e1 = np.degrees(np.arctan2(-R[2, 1], R[2, 0])); e3 = np.degrees(np.arctan2(-R[1, 2], -R[0, 2]))   # R = Rz(e3) Ry'(e2) Rz(e1): R[2,:] = (s2 c1, -s2 s1, c2), R[:,2] = (-s2 c3, -s2 s3, c2)
    else:
        e1 = np.degrees(np.arctan2(R[1, 0], R[0, 0])) * (1 if R[2, 2] > 0 else -1); e3 = 0.
    return np.array([e1 % 360., e2, e3 % 360.])


def rotvec_mats(step_deg):
    """27 small rotations: rotation vectors on the {-1,0,1}^3 grid times step (degrees)"""
    from scipy.spatial.transform import Rotation as Rot
    g = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)], float) * np.radians(step_deg)
    return Rot.from_rotvec(g).as_matrix()


def _pose_chunk(args):
    s0, s1, am, step, prior = args; F, C, R, eo = G['F'][s0:s1].astype(complex), G['C'][s0:s1].astype(float), G['R'][s0:s1].astype(float), G['eo'][s0:s1]
    REF = G['REF']; Bd = G['Bd']; conv = G['conv']; SC = G['SC']; sig2 = G['sig2']; dR = rotvec_mats(step); nc = len(dR)
    Rn = R.copy(); gain = np.zeros(s1 - s0)
    for k in np.unique(am):
        for h in (0, 1):
            m = np.where((am == k) & (eo == h))[0]
            if len(m) == 0: continue
            re, im = REF[1 - h][k]; Rcur = R[m].copy(); ll0 = None
            from scipy.spatial.transform import Rotation as Rot
            for stp in (step, step / 2.):     # posterior-MEAN rotation over the candidate grid (an argmax over candidates follows the noise)
                g = np.array([[i, j, kk] for i in (-1, 0, 1) for j in (-1, 0, 1) for kk in (-1, 0, 1)], float) * np.radians(stp); dR = Rot.from_rotvec(g).as_matrix()
                L = np.zeros((len(m), len(dR)))
                for q, d in enumerate(dR):
                    Rc = Rcur @ d; P = C[m] * slice_batch(re, im, Rc, Bd.kx, Bd.ky, conv) * SC; L[:, q] = -np.sum(Bd.wq * np.abs(F[m] - P) ** 2 / sig2, 1)
                if ll0 is None: ll0 = L[:, 13].copy()
                Lp = L - 0.5 * np.sum(g ** 2, 1)[None, :] / np.radians(prior) ** 2      # Gaussian prior on the update: flat evidence -> no move
                Wd = np.exp(Lp - Lp.max(1, keepdims=True)); Wd /= Wd.sum(1, keepdims=True); dv = Wd @ g
                Rcur = Rcur @ Rot.from_rotvec(dv).as_matrix()
            P = C[m] * slice_batch(re, im, Rcur, Bd.kx, Bd.ky, conv) * SC; ll1 = -np.sum(Bd.wq * np.abs(F[m] - P) ** 2 / sig2, 1)
            Rn[m] = Rcur; gain[m] = ll1 - ll0
    return Rn.astype(np.float32), gain


def grid_kernel(n1, n2, sigma, n3=1):
    """row-normalised Gaussian kernel on an n1 x n2 (x n3) lattice and the lattice coordinates of every cell"""
    i1, rem = np.divmod(np.arange(n1 * n2 * n3), n2 * n3); i2, i3 = np.divmod(rem, n3)
    d2 = (i1[:, None] - i1[None, :]) ** 2 + (i2[:, None] - i2[None, :]) ** 2 + (i3[:, None] - i3[None, :]) ** 2
    K = np.exp(-0.5 * d2 / max(sigma, 1e-6) ** 2); return K / K.sum(1, keepdims=True), i1, i2, i3


def _ginit(zz, n1, n2, n3):
    """initial grid positions from the first 2-3 columns of a latent: 5th-95th percentiles mapped onto the grid ranges"""
    def sc(v, nn): return np.clip((v - np.percentile(v, 5)) / max(np.percentile(v, 95) - np.percentile(v, 5), 1e-9), 0, 1) * (nn - 1)
    g3 = sc(zz[:, 2], n3) if (n3 > 1 and zz.shape[1] > 2) else np.zeros(len(zz))
    return sc(zz[:, 0], n1), sc(zz[:, 1], n2), g3


def knn_cv(z, y, k=15, nmax=60000, rng=None):
    idx = np.arange(len(y)) if len(y) <= nmax else rng.choice(len(y), nmax, replace=False); z, y = z[idx], y[idx]; h = np.arange(len(y)) % 2; acc = []
    for a in (0, 1):
        tr, te = h == a, h != a; kk = int(min(k, tr.sum())); nn = cKDTree(z[tr]).query(z[te], k=kk)[1]; nn = nn if kk > 1 else nn[:, None]; lab = y[tr][nn]
        acc.append(np.mean(np.array([np.bincount(r).argmax() for r in lab]) == y[te]))
    return float(np.mean(acc))


def auc(score, pos):
    r = np.argsort(np.argsort(score)) + 1; npos = pos.sum(); nneg = len(pos) - npos; return float((r[pos].sum() - npos * (npos + 1) / 2) / max(npos * nneg, 1))


def readout(tag, r, i1, i2, ngrid, conf, closed, axis, rng, cells_dom=None, layer_of=None, i3=None):
    rg = r[:, :ngrid]; mass = rg.sum(1); rn = rg / np.maximum(mass[:, None], 1e-12); z1 = rn @ i1; z2 = rn @ i2; z3 = (rn @ i3) if (i3 is not None and i3.max() > 0) else None
    good = conf > 0; ongrid = mass > 0.5; g = good & ongrid; ax = np.array([axis.get(c, 0.) for c in conf]); am = r.argmax(1)
    rho1, rho2 = spearmanr(z1[g], ax[g])[0], spearmanr(z2[g], ax[g])[0]
    Z = np.c_[z1, z2] if z3 is None else np.c_[z1, z2, z3]; nd = Z.shape[1]; h = np.arange(len(conf)) % 2; pred = np.zeros(len(conf))
    for s in (0, 1):
        tr, te = g & (h == s), g & (h != s); X = np.c_[Z[tr], np.ones(tr.sum())]; w = np.linalg.lstsq(X, ax[tr], rcond=None)[0]; pred[te] = np.c_[Z[te], np.ones(te.sum())] @ w
    rho_fit = spearmanr(pred[g], ax[g])[0]
    acc = knn_cv(Z[g], conf[g], nmax=60000, rng=rng); cls = np.unique(conf[g]); mu = np.array([Z[g][conf[g] == c].mean(0) for c in cls])
    Wc = sum(np.cov(Z[g][conf[g] == c].T) * np.sum(conf[g] == c) for c in cls) / g.sum(); Wi = np.linalg.inv(Wc + 1e-9 * np.eye(nd))
    dp = [np.sqrt((mu[a] - mu[b]) @ Wi @ (mu[a] - mu[b])) for a in range(len(cls)) for b in range(a + 1, len(cls))]
    cl = closed & g
    if cl.any():
        mu1, mu0 = Z[cl].mean(0), Z[g & ~closed].mean(0); wl = np.linalg.solve(Wc + 1e-9 * np.eye(nd), mu1 - mu0); sc = Z @ wl; a_cl = auc(sc[g], closed[g])
    else: a_cl = 0.
    # argmax-cell purity over conformer particles on the grid
    pur = sum(np.bincount(conf[g][am[g] == k]).max() for k in range(ngrid) if np.any(am[g] == k)) / max(g.sum(), 1)
    junk = ~good; jcap = np.mean(am[junk] >= ngrid) if junk.any() else 0.; ccap = np.mean(am[good] >= ngrid)
    ncells = np.sum(np.bincount(am[g], minlength=ngrid)[:ngrid] > 0.002 * g.sum())
    print('%s: Spearman(z1, axis) %+.3f (z2 %+.3f, fitted %+.3f) | 15-NN %.3f (chance %.3f) | pairwise d\' median %.2f | closed AUC %.3f | cell purity %.3f over %d used cells | junk in free comps %.2f (conformer particles %.2f) | mean max r %.2f' %
          (tag, rho1, rho2, rho_fit, acc, np.bincount(conf[g]).max() / g.sum(), np.median(dp), a_cl, pur, ncells, jcap, ccap, r.max(1).mean()), flush=True)
    print('%s   per-conformer mean (z1, z2%s): %s' % (tag, ', z3' if nd == 3 else '', ' '.join('%d:(%s)' % (c, ','.join('%.2f' % v for v in m)) for c, m in zip(cls, mu))), flush=True)
    if layer_of is not None and layer_of.max() > 0:
        lmass = np.stack([rg[:, layer_of == l].sum(1) for l in range(layer_of.max() + 1)], 1); lmass /= np.maximum(lmass.sum(1, keepdims=True), 1e-12)
        print('%s   layer mass: overall %s | per conformer (layer 2+ share) %s | junk %s' % (tag, np.array2string(lmass[good].mean(0), precision=2), ' '.join('%d:%.2f' % (c, lmass[good & (conf == c), 1:].sum(1).mean()) for c in cls), ('%.2f' % lmass[~good, 1:].sum(1).mean()) if (~good).any() else '-'), flush=True)
    if r.shape[1] > ngrid:
        fr = []
        for k in range(ngrid, r.shape[1]):
            m = am == k; nk = m.sum(); mj = m & good
            fr.append('free%d: n %d, junk %.2f, dominant conformer %s' % (k - ngrid + 1, nk, (m & ~good).sum() / max(nk, 1), ('%d (%.2f)' % (np.bincount(conf[mj]).argmax(), np.bincount(conf[mj]).max() / mj.sum())) if mj.any() else '-'))
        print('%s   %s' % (tag, ' | '.join(fr)), flush=True)
    return dict(rho1=float(rho1), rho2=float(rho2), rho_fit=float(rho_fit), knn=acc, dprime=float(np.median(dp)), closed_auc=a_cl, purity=float(pur), junk_capture=float(jcap), z1=z1, z2=z2, z3=z3, am=am)


def flex_basis(spec, n3, G):
    """orthonormal basis (n3^3 x npc, complex) of the flex PC maps on this band's coarse Fourier cube, cached per n3"""
    key = ('FLEXB', n3)
    if key in G: return G[key]
    d, npc = (spec.split(':') + ['16'])[:2]; npc = int(npc); cols = []
    for i in range(1, npc + 1):
        f = '%s/flex_pca_fitB_pc%03d.mrc' % (d, i)
        if not os.path.exists(f): break
        v = read_mrc(f); nb = v.shape[0]; Fm = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v))); c_ = nb // 2; h_ = n3 // 2
        if nb < n3: Fp = np.zeros((n3, n3, n3), complex); Fp[h_ - c_:h_ - c_ + nb, h_ - c_:h_ - c_ + nb, h_ - c_:h_ - c_ + nb] = Fm; Fm = Fp
        else: Fm = Fm[c_ - h_:c_ - h_ + n3, c_ - h_:c_ - h_ + n3, c_ - h_:c_ - h_ + n3]
        cols.append(Fm.ravel())
    Q, _ = np.linalg.qr(np.array(cols).T); G[key] = Q
    print('   flex prior: %d PC maps -> orthonormal basis on the %d^3 cube' % (Q.shape[1], n3), flush=True); return Q

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('stk'); ap.add_argument('table'); ap.add_argument('gt'); ap.add_argument('out')
    ap.add_argument('--cons', required=True); ap.add_argument('--flexdir', default=''); ap.add_argument('--states', required=True); ap.add_argument('--axis', required=True)
    ap.add_argument('--n1', type=int, default=7); ap.add_argument('--n2', type=int, default=7); ap.add_argument('--n3', type=int, default=1, help='third latent dimension (1 = the 2-D grid)'); ap.add_argument('--nfree', type=int, default=2)
    ap.add_argument('--lps', default='20 20 20 14 14 14 10 10 10'); ap.add_argument('--sigmas', default='1.5 1.5 1.2 1.0 1.0 0.8 0.7 0.7 0.7')
    ap.add_argument('--nsub', type=int, default=0); ap.add_argument('--nproc', type=int, default=8); ap.add_argument('--init', default='pca'); ap.add_argument('--mskdiam', type=float, default=200.)
    ap.add_argument('--node_basis', type=int, default=0, help='q > 0: mixture of PPCAs: one q-dim basis per board (= node) fitted by PPCA-EM to the per-image residuals after the cells (own half), cross-half validated; the E-step uses the evidence marginalised over the basis coordinate (--node_basis_marg 1) and the posterior-mean coordinate per image is written as zb columns')
    ap.add_argument('--node_basis_after', type=int, default=3); ap.add_argument('--node_basis_marg', type=int, default=1); ap.add_argument('--node_basis_rep', type=float, default=0.5); ap.add_argument('--node_basis_deflate', type=int, default=1)
    ap.add_argument('--flex_prior', default='', help='DIR[:NPC]: flex_pca PC maps (flex_pca_fitB_pc*.mrc) as the shrinkage target: every cell is shrunk toward the consensus + its own projection onto the flex subspace instead of toward the consensus alone (PPCA prior on the free states)');
    ap.add_argument('--shrink', type=float, default=0.02); ap.add_argument('--tag', default='klem'); ap.add_argument('--chunk', type=int, default=1000)
    ap.add_argument('--nkseed', type=int, default=0, help='extra free components seeded from k-means states on the flex latent (candidates for discrete / compositional states)')
    ap.add_argument('--nlayer', type=int, default=1, help='parallel latent layers (grids) with no kernel smoothing between them: a compositional state can occupy its own layer')
    ap.add_argument('--bags', default='', help='2D sub-class labels per project particle (print_project_field ptcl2D with class=, or one int per row): members of a bag share a soft prior = their mean log-likelihood x bag_beta')
    ap.add_argument('--bag_beta', type=float, default=2., help='strength of the bag prior (0 = per-particle E-step)')
    ap.add_argument('--split_rounds', type=int, default=0, help='after the schedule: duplicate the grid into 2 layers with a RANDOM per-particle split and run this many more rounds at the last band; a layer difference that reproduces across halves is a discrete factor (composition) the smooth grid cannot hold')
    ap.add_argument('--avgs', default='', help='cls_expansion_cavgs.mrc of the --bags sub-classes (even/odd stacks beside it): sub-class-level E-step, members inherit the responsibilities')
    ap.add_argument('--avg_mode', default='unmix', help="'pool': members take the sub-class's hard responsibilities; 'unmix': the average is unmixed onto the cell deviations (NNLS) and the proportions act as a prior on the members' own likelihoods")
    ap.add_argument('--avg_prior', type=float, default=1., help='weight of the log mixture prior in unmix mode')
    ap.add_argument('--avg_npose', type=int, default=120, help='members per bag (from ALL members) whose poses define the pose-mean renders of the sub-class E-step')
    ap.add_argument('--avg_mix', type=float, default=0., help='0 = pure sub-class responsibilities for bagged particles; m in (0,1] = geometric mix with the particle responsibilities')
    ap.add_argument('--kernel', default='lattice', help="'lattice' (the n1 x n2 board) or 'graph': from iteration --graph_after on, the smoothing kernel is computed from the states themselves (cross-half correlation of the state deviations), the board survives only as a spectral embedding for display")
    ap.add_argument('--graph_after', type=int, default=3); ap.add_argument('--graph_nn', type=float, default=4., help='effective neighbourhood size (states) of the graph kernel at the end of the schedule; scaled up early in proportion to the lattice sigma')
    ap.add_argument('--devw', type=float, default=0., help='> 0: weight every shell of the E-step likelihood by the cross-half FSC of the cell deviations, 2F/(1+F) floored at this value (0 = off)')
    ap.add_argument('--focus', type=float, default=0., help='> 0: E-step references = consensus + deviation masked by the reproducible variability of the cells (soft mask from the cross-half product, scaled by this percentile of its 99th; 0 = off)')
    ap.add_argument('--focus_mask', default='', help='focal classification: a real-space mask (MRC, full box, 0-1) resampled to the coarse grid at every band and applied to the E-step deviations instead of the data-derived focus mask (consensus outside it): the cells then differ only inside the mask')
    ap.add_argument('--pose_marg', type=float, default=0., help='> 0: marginalise the E-step likelihood over the pose and six tilts of this many degrees about three axes (logsumexp)')
    ap.add_argument('--pose_marg_n', type=int, default=7, help='7 (pose + six tilts) or 27 (the full {-1,0,1}^3 rotation-vector grid)')
    ap.add_argument('--shift_marg', type=float, default=0., help='> 0: also marginalise over five shifts (0, +-x, +-y) of this many pixels')
    ap.add_argument('--pose_marg_mstep', action='store_true', help='backproject every particle at its neighbourhood poses with the joint pose-state posterior (pose-marginalised M-step)')
    ap.add_argument('--temp_start', type=float, default=1., help='> 1: deterministic annealing, E-step temperature from this value to 1 over the first two thirds of the schedule')
    ap.add_argument('--pnoise_m', action='store_true', help='use the per-particle noise factor as the Wiener weight in the M-step too')
    ap.add_argument('--temp_final', type=float, default=1., help='< 1: temper (sharpen) the E-step in the last three iterations, log-likelihoods divided by this')
    ap.add_argument('--pnoise', action='store_true', help='per-particle noise scale from the residual against the consensus')
    ap.add_argument('--pscale', action='store_true', help='per-particle contrast scale against the consensus, folded into the CTF weights')
    ap.add_argument('--mixprop', action='store_true', help='estimate the mixing proportions (state fractions) instead of uniform mixing')
    ap.add_argument('--adapt', type=float, default=0., help='adaptive kernel width: cells whose deviation reproduces across halves get sigma x (1 - adapt), the least reproducible sigma x (1 + adapt), linear in the cross-half correlation between 0.3 and 0.8 (0 = fixed width)')
    ap.add_argument('--subset', default='', help='file with one 1-based particle index per line, or a klem_latent.txt: restricts the EM to these particles (instead of a random --nsub draw)')
    ap.add_argument('--bag_parents', default='', help='parent-class labels per project particle (same format as --bags): the bag prior is centred on the parent mean, so only the within-parent (within-view) part of the bag evidence is shared')
    ap.add_argument('--deconv', type=float, default=0., help='one-step deconvolution of the grid kernel on the cell volumes after every M-step: V + deconv x (V - K V)')
    ap.add_argument('--weights', default='', help='soft families: file with "pind weight" lines, the parent family responsibility of every particle; weights multiply the responsibilities in the M-step (consensus and cells), the E-step is unweighted')
    ap.add_argument('--balance', type=float, default=0., help='> 0: optimal-transport balanced E-step (Sinkhorn in log space): the grid cells are constrained to equal shares of the grid mass (free components unconstrained), the column potentials scaled by this value (1 = fully balanced)')
    ap.add_argument('--balance_until', type=int, default=0, help='apply the balance only in the first this many iterations (0 = every iteration)')
    ap.add_argument('--final_all', action='store_true', help='after an EM on a --nsub subset, one E-step over ALL particles at the last band -> latent for every particle')
    ap.add_argument('--tie_layers', action='store_true', help='with nlayer 2: layer 2 = layer 1 + ONE shared additive component (occupancy model): the base map pools all particles, the component is learnt once')
    ap.add_argument('--tie_support', type=float, default=0., help='tied layers: the component is a POSITIVE added density (negative real-space values clipped); from round tie_support_after its support is restricted to the largest connected components of the cross-half reproducible positive part above this fraction of the maximum (0 = unconstrained)')
    ap.add_argument('--tie_support_after', type=int, default=2); ap.add_argument('--tie_ncomp', type=int, default=1, help='tied layers: number of connected components kept in the support')
    ap.add_argument('--tie_support_mask', default='', help='tied layers: fixed full-box support mask of the component (a candidate site); positivity applies inside it')
    ap.add_argument('--pose_every', type=int, default=0, help='refine poses every k iterations (0 = never) against the argmax component volume')
    ap.add_argument('--pose_step', type=float, default=2., help='rotation grid step (deg); candidates = 3x3x3 rotation-vector grid x 2 refinements'); ap.add_argument('--pose_truth', default='', help='toy: true-pose table for the pose-error readout')
    ap.add_argument('--pose_prior', type=float, default=3., help='Gaussian prior width (deg) on the pose update per step: a flat likelihood then leaves the pose where it is')
    ap.add_argument('--pose_mask', default='', help="node alignment policy: mask applied in real space (coarse grid) to the references of the pose step: '' = none (everything), 'core' = density x (1 - reproducible variability of the cells), 'region' = that variability, or an MRC file (full box) Fourier-cropped to the grid")
    ap.add_argument('--pose_ref', default='cell', choices=['cell', 'consensus'], help='reference of the pose step: the particle argmax cell (lever arm, region-aware) or the consensus (core alignment)')
    ap.add_argument('--meta', default='', help='real data: metadata npz (e1 e2 e3 x y dfx dfy angast eo flags stkind indstk [z key]); STK then = comma list of stack files, TABLE/GT ignored')
    ap.add_argument('--zkey', default='z692', help='latent key in --meta for --init npz (rows = flags==1 particles)'); ap.add_argument('--ctf', default='auto', help='signed|abs|auto (auto: tested at 12 A)')
    ap.add_argument('--labels', default='', help='external per-particle labels (one int per flags==1 particle) for sanity readouts on real data')
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); t0 = time.time(); T = a.tag; rng = np.random.default_rng(0)
    axis = parse_axis(a); sid = iid = None; MZ = None
    if a.meta:
        M = np.load(a.meta); flags = M['flags'].astype(int); N = len(flags)
        tab = np.c_[M['e1'], M['e2'], M['e3'], M['x'], M['y'], M['dfx'], M['dfy'], M['angast'], M['eo'] if 'eo' in M.files else -np.ones(N), (flags == 1).astype(float)]
        files = a.stk.split(','); stk = []
        for f in files:
            h = np.fromfile(f, np.int32, 256); box = int(h[0]); stk.append(np.memmap(f, np.float32, 'r', offset=1024 + int(h[23]), shape=(int(h[2]), box, box)))
        sid = M['stkind'].astype(int) - 1; iid = M['indstk'].astype(int) - 1; assert sid.max() < len(stk) and all(iid[sid == k].max() < stk[k].shape[0] for k in range(len(stk)))
        gt = {}; lab_ext = np.loadtxt(a.labels).astype(int) if a.labels else None
        if lab_ext is not None and len(lab_ext) == N: lab_ext = lab_ext[flags == 1]   # labels given for all particles -> the refined subset
        if lab_ext is not None and len(lab_ext) != int((flags == 1).sum()): print('%s: label file has %d rows, %d refined particles -> labels ignored' % (T, len(lab_ext), int((flags == 1).sum())), flush=True); lab_ext = None
        if a.zkey in M.files: MZ = M[a.zkey]
    else:
        gt = read_gt(a.gt); tab = read_table(a.table, ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast', 'eo', 'state')); N = len(tab)
        h = np.fromfile(a.stk, np.int32, 256); box = int(h[0]); stk = np.memmap(a.stk, np.float32, 'r', offset=1024 + int(h[23]), shape=(int(h[2]), box, box)); lab_ext = None
    eo_all = tab[:, 8].astype(int); noeo = eo_all < 0; eo_all[noeo] = np.arange(N)[noeo] % 2
    if len(np.unique(eo_all)) < 2: eo_all = np.arange(N) % 2; print('%s: even/odd flag constant in the table -> particle parity used' % T, flush=True)
    sel = np.where(tab[:, 9] > 0)[0]; sel_full = sel.copy()
    if a.subset:
        sub = np.loadtxt(a.subset, comments='#'); sub = (sub[:, 0] if sub.ndim == 2 else sub).astype(int) - 1; sel = np.intersect1d(sel, sub); print('%s: subset %s: %d particles' % (T, a.subset, len(sel)), flush=True)
    if a.nsub and len(sel) > a.nsub: sel = np.sort(rng.choice(sel, a.nsub, replace=False))
    n = len(sel); bag = None; par = None; pw = np.ones(n)
    if a.weights:
        wmap = {}
        for l in open(a.weights):
            if l.strip() and not l.startswith('#'): t = l.split(); wmap[int(float(t[0]))] = float(t[1])
        pw = np.array([wmap.get(int(i) + 1, 1.) for i in sel]); print('%s: soft-family weights from %s: mean %.2f, 5-95 pct %.2f-%.2f, %d of %d particles below 0.5' % (T, a.weights, pw.mean(), np.percentile(pw, 5), np.percentile(pw, 95), int((pw < 0.5).sum()), n), flush=True)
    def read_label_file(path):   # print_project_field ptcl2D (class=) or one int per row
        bl = []
        for l in open(path):
            if not l.strip() or l.startswith('#'): continue
            kv = dict(t.split('=', 1) for t in l.split() if '=' in t); bl.append(int(float(kv['class'])) if 'class' in kv else int(float(l.split()[0])))
        bl = np.array(bl); return bl[sel] if len(bl) >= sel.max() + 1 else None
    if a.bags:
        bag = read_label_file(a.bags)
        if bag is not None: print('%s: bags from %s: %d bags over %d particles (bag_beta %.1f)' % (T, a.bags, len(np.unique(bag[bag > 0])), int((bag > 0).sum()), a.bag_beta), flush=True)
    if a.avgs and bag is not None:   # sub-class E-step: particles ordered by bag (unbagged first) so that every chunk holds whole bags
        order = np.argsort(np.where(bag > 0, bag, 0), kind='stable'); sel = sel[order]; bag = bag[order]; pw = pw[order]
    if a.bag_parents and bag is not None:
        par = read_label_file(a.bag_parents)
        if par is not None: print('%s: bag parents from %s: %d parents; the bag prior is centred on the parent mean (within-view sharing only)' % (T, a.bag_parents, len(np.unique(par[par > 0]))), flush=True)
    if a.meta:
        pos_in_full = {int(i): j for j, i in enumerate(sel_full)}; rows = np.array([pos_in_full[int(i)] for i in sel])
        conf = (lab_ext[rows] if lab_ext is not None else np.ones(n, int)); closed = np.zeros(n, bool)   # label 0 = unlabelled (treated like junk in the readouts)
        if lab_ext is not None: print('%s: external labels used as pseudo-conformers for the sanity readouts (sizes %s)' % (T, np.bincount(conf)[1:]), flush=True)
    else:
        conf = np.array([gt.get(i + 1, (0, 0))[0] for i in sel]); closed = np.array([gt.get(i + 1, (0, 0))[1] for i in sel]) == 1
    eo = eo_all[sel]
    ncell = a.n1 * a.n2 * a.n3; ngrid = ncell * a.nlayer; nfree_tot = a.nfree + a.nkseed; K = ngrid + nfree_tot; lps = [float(v) for v in a.lps.split()]; sigmas = [float(v) for v in a.sigmas.split()]; assert len(lps) == len(sigmas)
    print('%s: %d particles (%d conformer, %d junk), grid %dx%d + %d free, schedule %s' % (T, n, (conf > 0).sum(), (conf == 0).sum(), a.n1, a.n2, a.nfree, ' '.join('%g@%g' % (l, s) for l, s in zip(lps, sigmas))), flush=True)
    # ---- calibration against the SIMPLE consensus (projection convention, shift sign/scale, contrast sign) at lp 20
    S = Shells(box, SMPD); msk = S.soft_mask(a.mskdiam / 2 / SMPD); VC = fftc(read_mrc(a.cons) * msk); VC = (np.ascontiguousarray(VC.real, np.float32), np.ascontiguousarray(VC.imag, np.float32))
    kmax0 = int(np.floor(box * SMPD / 20.)); Bd0 = Band(box, kmax0); c = box // 2
    cal = rng.choice(sel[conf > 0], min(300, n), replace=False); conv, ss, fs, sign = calibrate(stk, tab, box, Bd0, (Bd0.kx + c).astype(int), (Bd0.ky + c).astype(int), VC, cal, sid=sid, iid=iid)
    ctfabs = a.ctf == 'abs'
    if a.ctf == 'auto':   # signed vs |CTF| decided at 12 A (band crosses the first zero) by the correlation with the consensus projection
        kmax12 = int(np.floor(box * SMPD / 12.)); Bd12 = Band(box, kmax12); ix12 = (Bd12.kx + c).astype(int); iy12 = (Bd12.ky + c).astype(int); cc = {}
        for mode in (False, True):
            D12 = Data(stk, tab, box, Bd12, ix12, iy12, sign, sid=sid, iid=iid, ctfabs=mode); F, C, R = D12.chunk(cal); P = C * slice_batch(VC[0], VC[1], R, Bd12.kx, Bd12.ky, conv) * shift_phase(Bd12, tab[cal, 3], tab[cal, 4], box, ss, fs)
            m = Bd12.shell >= int(np.floor(box * SMPD / 16.)); cc[mode] = float(np.mean(np.real(np.sum(np.conj(F[:, m]) * P[:, m], 1)) / np.sqrt(np.sum(np.abs(F[:, m]) ** 2, 1) * np.sum(np.abs(P[:, m]) ** 2, 1))))
        ctfabs = cc[True] > cc[False]; print('%s CTF model at 16-12 A: signed ncc %.4f, |CTF| ncc %.4f -> %s' % (T, cc[False], cc[True], 'abs' if ctfabs else 'signed'), flush=True)
    # ---- initial responsibilities
    _, i1c, i2c, i3c = grid_kernel(a.n1, a.n2, 1., a.n3); i1 = np.tile(i1c, a.nlayer); i2 = np.tile(i2c, a.nlayer); i3 = np.tile(i3c, a.nlayer); layer_of = np.repeat(np.arange(a.nlayer), ncell)
    r = np.zeros((n, K))
    if a.init == 'npz' and MZ is not None:
        zs = (MZ - MZ.mean(0)) / np.maximum(MZ.std(0), 1e-12); zz = zs[rows][:, :3]
        print('%s init: npz latent %s, PCs 1-2 (std per dim %s)' % (T, a.zkey, np.array2string(MZ.std(0)[:6], precision=3)), flush=True)
        g1, g2, g3 = _ginit(zz, a.n1, a.n2, a.n3)
    elif a.init in ('pca', 'pca12', 'spectral') and a.flexdir:
        from pool_place import read_embedding
        _, pinds, z, _ = read_embedding(a.flexdir + '/flex_pca_embedding.bin'); zs = (z - z.mean(0)) / np.maximum(z.std(0), 1e-12)
        ctr = {}
        try:
            for l in open(a.flexdir + '/flex_pca_contrast.txt'):
                if l.strip() and not l.startswith('#'): w = l.split(); ctr[int(w[0])] = float(w[1])
        except Exception: pass
        cvec = np.array([ctr.get(int(p), 0.) for p in pinds]); cc = np.array([abs(np.corrcoef(zs[:, q], cvec)[0, 1]) if cvec.std() > 0 else 0. for q in range(zs.shape[1])])
        keepq = [q for q in range(zs.shape[1]) if cc[q] < 0.3] if a.init != 'pca12' else list(range(zs.shape[1])); keep = keepq[:3]   # pca12: raw PCs 1-2 (a deleted domain lowers the contrast, so the contrast filter can drop the compositional axis)
        if a.init == 'spectral':
            # diffusion-map coordinates: symmetric kNN graph on the contrast-free latent of a 40k subsample, normalised
            # Laplacian, eigenvectors 2-3; the rest of the particles placed by their nearest subsample neighbours
            from scipy.sparse import csr_matrix, diags
            from scipy.sparse.linalg import eigsh
            Zq = zs[:, keepq]; isub = rng.choice(len(Zq), min(40000, len(Zq)), replace=False); Zs = Zq[isub]; kn = 20
            dd, nn = cKDTree(Zs).query(Zs, k=kn + 1); sig_g = np.median(dd[:, 1:]) ** 2
            rows = np.repeat(np.arange(len(Zs)), kn); cols = nn[:, 1:].ravel(); wts = np.exp(-dd[:, 1:].ravel() ** 2 / (2 * sig_g))
            Wg = csr_matrix((wts, (rows, cols)), shape=(len(Zs), len(Zs))); Wg = Wg.maximum(Wg.T); dg = np.asarray(Wg.sum(1)).ravel()
            Dm = diags(1. / np.sqrt(dg)); Lsym = Dm @ Wg @ Dm; vals, vecs = eigsh(Lsym, k=4, which='LA'); order = np.argsort(vals)[::-1]; vecs = vecs[:, order] / np.sqrt(dg)[:, None]
            emb = vecs[:, 1:4]; print('%s init: spectral (kNN %d on %d contrast-free PCs, eigenvalues %s)' % (T, kn, len(keepq), np.array2string(vals[order][:4], precision=4)), flush=True)
            _, nn_all = cKDTree(Zs).query(Zq, k=5); zz_all = emb[nn_all].mean(1); pos = {int(p): zz_all[j] for j, p in enumerate(pinds)}
        else:
            print('%s init: contrast correlation per PC %s -> using PCs %s' % (T, np.array2string(cc[:8], precision=2), [q + 1 for q in keep]), flush=True)
            pos = {int(p): zs[j, keep] for j, p in enumerate(pinds)}
        zz = np.array([pos.get(i + 1, np.zeros(3)) for i in sel])
        g1, g2, g3 = _ginit(zz, a.n1, a.n2, a.n3)
    else:
        g1 = rng.uniform(0, a.n1 - 1, n); g2 = rng.uniform(0, a.n2 - 1, n); g3 = rng.uniform(0, a.n3 - 1, n)
    r0 = np.exp(-0.5 * ((g1[:, None] - i1c[None, :]) ** 2 + (g2[:, None] - i2c[None, :]) ** 2 + (g3[:, None] - i3c[None, :]) ** 2) / 0.8 ** 2)
    if a.nlayer > 1:
        pref = rng.dirichlet(np.ones(a.nlayer) * 3., size=n)   # random per-particle layer preference breaks the symmetry
        r[:, :ngrid] = np.concatenate([r0 * pref[:, [l]] for l in range(a.nlayer)], axis=1)
    else: r[:, :ngrid] = r0
    r[:, ngrid:ngrid + a.nfree] = 0.02 * r[:, :ngrid].sum(1, keepdims=True) / max(a.nfree, 1)
    if a.nkseed > 0 and a.flexdir:
        from pool_place import read_embedding as _re, kmeans as _km
        _, pk, zk, _ = _re(a.flexdir + '/flex_pca_embedding.bin'); zk = (zk - zk.mean(0)) / np.maximum(zk.std(0), 1e-12); nrm = np.linalg.norm(zk, axis=1); core = np.where(nrm <= np.percentile(nrm, 95))[0]
        lab_c, cent, _ = _km(zk[rng.choice(core, min(60000, len(core)), replace=False)], a.nkseed, seeds=10, rng=rng); lab_all = ((zk[:, None, :] - cent[None]) ** 2).sum(2).argmin(1)
        labmap = dict(zip(pk, lab_all)); lab_sel = np.array([labmap.get(i + 1, -1) for i in sel])
        for k in range(a.nkseed):
            m = lab_sel == k; r[m, ngrid + a.nfree + k] = 0.5 * r[m, :ngrid].sum(1)   # half of the mass starts on the k-means state, half on the grid
        print('%s init: %d k-means seeded free components, sizes %s' % (T, a.nkseed, np.bincount(lab_sel[lab_sel >= 0], minlength=a.nkseed)), flush=True)
    r /= r.sum(1, keepdims=True)
    readout(T + ' init', r, i1, i2, ngrid, conf, closed, axis, rng, layer_of=layer_of, i3=i3)
    split_at = len(lps) if a.split_rounds > 0 and a.nlayer == 1 else -1
    if split_at >= 0: lps = lps + [lps[-1]] * a.split_rounds; sigmas = sigmas + [sigmas[-1]] * a.split_rounds
    # ---- EM over the schedule
    G.update(D=None, tab=tab, ss=ss, fs=fs, box=box, sel=sel, conv=conv); cur_lp = None; ctx = mp.get_context('fork'); hist = []; FLIB = fort_lib(a.nproc)
    chunks = [(s0, min(s0 + a.chunk, n)) for s0 in range(0, n, a.chunk)]
    avg_mode = bool(a.avgs) and bag is not None
    if avg_mode:   # chunk boundaries only where the bag changes (unbagged particles may split anywhere)
        cuts = [0]; last = 0
        for i in range(1, n):
            if i - last >= a.chunk and (bag[i] <= 0 or bag[i] != bag[i - 1]): cuts.append(i); last = i
        cuts.append(n); chunks = [(cuts[j], cuts[j + 1]) for j in range(len(cuts) - 1)]
        G['bag'] = bag; G['AVGMODE'] = a.avg_mode
        if a.avg_mode == 'unmixc':
            assert par is not None, 'unmixc needs --bag_parents'; G['PAROF'] = {int(b): int(par[bag == b][0]) for b in np.unique(bag[bag > 0]) if (par[bag == b] > 0).any()}
        h_ = np.fromfile(a.avgs, np.int32, 256); navg = int(h_[2]); AV = {}
        for suf, key in (('_even', 'E'), ('_odd', 'O')):
            fn = a.avgs.replace('.mrc', suf + '.mrc'); fn = fn if os.path.exists(fn) else a.avgs
            hh_ = np.fromfile(fn, np.int32, 256); AV[key] = np.memmap(fn, np.float32, 'r', offset=1024 + int(hh_[23]), shape=(int(hh_[2]), box, box))
        ub = [int(b) for b in np.unique(bag[bag > 0]) if b <= navg]; G['AROW'] = {b: j for j, b in enumerate(ub)}
        G['INPL'] = None; G['INPLC'] = None; G['INPL_ALL'] = None
        # pose samples per bag and half from ALL members with state > 0 (the EM subset may hold only a few members of each bag)
        bl_all = []
        for l in open(a.bags):
            if not l.strip() or l.startswith('#'): continue
            kv = dict(t.split('=', 1) for t in l.split() if '=' in t); bl_all.append(int(float(kv['class'])) if 'class' in kv else int(float(l.split()[0])))
        bl_all = np.array(bl_all); G['BAGROWS'] = {}
        for b in ub:
            rows_b = np.where((bl_all == b) & (tab[:len(bl_all), 9] > 0))[0]
            if len(rows_b) > a.avg_npose: rows_b = np.sort(rng.choice(rows_b, a.avg_npose, replace=False))
            G['BAGROWS'][b] = rows_b
        print('%s: pose samples for the renders: up to %d members per bag from all %d labelled particles (median %d per bag)' % (T, a.avg_npose, int((bl_all > 0).sum()), int(np.median([len(v) for v in G['BAGROWS'].values()]))), flush=True)
        with open(a.bags) as fh: first = next((l for l in fh if l.strip() and not l.startswith('#')), '')
        if 'e3=' in first and 'x=' in first:   # averages in the 2D class frame: members enter rotated/shifted by their ptcl2D alignment
            E3, X2, Y2 = [], [], []
            for l in open(a.bags):
                if not l.strip() or l.startswith('#'): continue
                kv = dict(t.split('=', 1) for t in l.split() if '=' in t); E3.append(float(kv.get('e3', 0))); X2.append(float(kv.get('x', 0))); Y2.append(float(kv.get('y', 0)))
            G['INPL_ALL'] = (np.array(E3), np.array(X2), np.array(Y2)); E3, X2, Y2 = (np.array(v)[sel] for v in (E3, X2, Y2)); G['INPL'] = (E3, X2, Y2)
            print('%s: the bags file carries 2D in-plane parameters (e3 sd %.0f deg, |shift| median %.1f px): class-frame renders, convention calibrated at the first band' % (T, E3.std(), np.median(np.hypot(X2, Y2))), flush=True)
        print('%s: sub-class E-step on %d averages (%s; even/odd %s), %d chunks aligned to bag boundaries, avg_mix %.2f' % (T, len(ub), a.avgs, 'yes' if os.path.exists(a.avgs.replace('.mrc', '_even.mrc')) else 'NO (same stack)', len(chunks), a.avg_mix), flush=True)
    for it, (lp, sig) in enumerate(zip(lps, sigmas)):
        if it == split_at:   # ---- split: every grid cell gets a twin in a second layer, each particle's grid mass divided at random between the layers
            a.nlayer = 2; ngrid = 2 * ncell; K = ngrid + nfree_tot; i1 = np.tile(i1c, 2); i2 = np.tile(i2c, 2); i3 = np.tile(i3c, 2); layer_of = np.repeat(np.arange(2), ncell)
            sp = rng.random(n); r = np.c_[r[:, :ncell] * (1 - sp)[:, None], r[:, :ncell] * sp[:, None], r[:, ncell:]]
            print('%s split: grid duplicated into 2 layers with a random per-particle split, %d more rounds at lp %g (a layer difference that reproduces across halves = a discrete factor)' % (T, a.split_rounds, lp), flush=True)
        if lp != cur_lp:
            kmax = int(np.floor(box * SMPD / lp)); Bd = Band(box, kmax); ix = (Bd.kx + c).astype(int); iy = (Bd.ky + c).astype(int)
            G.update(D=Data(stk, tab, box, Bd, ix, iy, sign, sid=sid, iid=iid, ctfabs=ctfabs), Bd=Bd, kmax=kmax)
            with ctx.Pool(a.nproc) as pool: parts = pool.map(_load_chunk, chunks)
            G['F'] = np.concatenate([p[0] for p in parts]); G['C'] = np.concatenate([p[1] for p in parts]); G['R'] = np.concatenate([p[2] for p in parts]); G['eo'] = eo; del parts
            # consensus of each half (all particles, weight 1) -> shrinkage target, noise model, reference scale
            numc = denc = None
            with ctx.Pool(a.nproc) as pool:
                if FLIB is not None: numc, denc = fort_mstep(FLIB, chunks, lambda s0, s1: np.ascontiguousarray(pw[s0:s1, None]), 1)
                else:
                    for p in pool.imap_unordered(_mstep_chunk, [(s0, s1, pw[s0:s1, None].copy()) for s0, s1 in chunks]): numc = p[0] if numc is None else numc + p[0]; denc = p[1] if denc is None else denc + p[1]
            n3 = 2 * kmax + 1
            eps0 = 0.05 * np.median(denc[denc > 0]); Gc = [(numc[hh, 0] / (denc[hh, 0] + eps0)).reshape(n3, n3, n3) for hh in (0, 1)]; eps = a.shrink * np.median(denc[denc > 0])
            G['REF'] = [[(np.ascontiguousarray(Gc[hh].real, np.float32), np.ascontiguousarray(Gc[hh].imag, np.float32))] for hh in (0, 1)]; G['SC'] = np.ones(Bd.n); G['sig2'] = np.ones(Bd.n)
            sub = rng.choice(n, min(2000, n), replace=False); Fs, Cs, Rs = G['F'][sub].astype(complex), G['C'][sub].astype(float), G['R'][sub].astype(float)
            P = np.zeros_like(Fs)
            for hh in (0, 1):
                m = eo[sub] == hh
                if not m.any(): continue
                re, im = G['REF'][1 - hh][0]; P[m] = Cs[m] * slice_batch(re, im, Rs[m], Bd.kx, Bd.ky, conv)
            shells = [Bd.shell == s for s in range(kmax + 1)]
            num = np.array([np.sum(Bd.wq[m] * np.real(np.conj(P[:, m]) * Fs[:, m])) if m.any() else 0. for m in shells]); den = np.array([np.sum(Bd.wq[m] * np.abs(P[:, m]) ** 2) if m.any() else 0. for m in shells])
            sc = np.where(den > 0, num / np.maximum(den, 1e-30), 0.); sc = np.convolve(np.pad(sc, 1, mode='edge'), np.ones(3) / 3., mode='valid'); sc = np.maximum(sc, 0.05 * max(sc[1:5].mean(), 1e-30)); SC = sc[Bd.shell]; P *= SC
            res2 = np.abs(Fs - P) ** 2; sig2 = np.maximum(np.array([res2[:, m].mean() if m.any() else 1. for m in shells]), 1e-12)[Bd.shell]; G['SC'] = SC; G['sig2'] = sig2; G['sig2_base'] = sig2.copy()
            G['pfac'] = None
            if a.pscale or a.pnoise:
                G['GC'] = [(np.ascontiguousarray(Gc[hh].real, np.float32), np.ascontiguousarray(Gc[hh].imag, np.float32)) for hh in (0, 1)]
                pn, pd, pr = _consensus_stats(chunks); sc_i = np.clip(pn / np.maximum(pd, 1e-30), 0.3, 3.)
                if a.pscale: G['C'] = (G['C'] * sc_i[:, None]).astype(G['C'].dtype); print('%s   per-particle contrast scale: median %.2f, 5-95 pct %.2f-%.2f (folded into the CTF weights)' % (T, np.median(sc_i), np.percentile(sc_i, 5), np.percentile(sc_i, 95)), flush=True)
                if a.pnoise: pf = np.clip(pr / max(np.median(pr), 1e-30), 0.3, 3.); G['pfac'] = pf; print('%s   per-particle noise factor: 5-95 pct %.2f-%.2f' % (T, np.percentile(pf, 5), np.percentile(pf, 95)), flush=True)
            snr = [float(np.mean(np.abs(P[:, m]) ** 2) / sig2[m][0]) for m in shells[1:]]
            print('%s band %g A: %d shells, %d coefficients, grid %d^3; reference scale %s; signal fraction per shell %s (%.0f s)' % (T, lp, kmax, Bd.n, n3, np.array2string(sc[1::max(1, kmax // 6)] / sc[1], precision=2), np.array2string(np.array(snr[::max(1, kmax // 6)]), precision=3), time.time() - t0), flush=True)
            if avg_mode:
                ub = sorted(G['AROW'], key=G['AROW'].get); nsh = Bd.shell.max() + 1; AE = np.zeros((len(ub), Bd.n), np.complex64); AO = np.zeros_like(AE); ASIG = np.zeros((len(ub), nsh))
                for j, b in enumerate(ub):
                    fe = fft2c(np.array(AV['E'][b - 1], float))[iy, ix] * sign; fo = fft2c(np.array(AV['O'][b - 1], float))[iy, ix] * sign; AE[j] = fe; AO[j] = fo
                    d2 = np.abs(fe - fo) ** 2 / 2.; ASIG[j] = np.bincount(Bd.shell, Bd.wq * d2, minlength=nsh) / np.maximum(np.bincount(Bd.shell, Bd.wq, minlength=nsh), 1)
                    ASIG[j] = np.maximum(np.convolve(np.pad(ASIG[j], 1, mode='edge'), np.ones(3) / 3., mode='valid'), 1e-6 * ASIG[j].max())
                G['AE'], G['AO'], G['ASIG'] = AE, AO, ASIG; G['GC'] = [(np.ascontiguousarray(Gc[hh].real, np.float32), np.ascontiguousarray(Gc[hh].imag, np.float32)) for hh in (0, 1)]
                if G['INPL'] is not None and G['INPLC'] is None:   # 12 candidate conventions, scored by the correlation of the consensus render with the averages
                    cnt = {b: 0 for b in ub}
                    for b in bag[bag > 0]: cnt[int(b)] = cnt.get(int(b), 0) + 1
                    cal_b = sorted(ub, key=lambda b: -cnt[b])[:30]; best = None; table = []
                    for sr in (1, -1):
                        for s2 in (1, -1):
                            for rot3 in (0, 1, -1):
                                cc = []
                                for b in cal_b:
                                    m = np.where(bag == b)[0]; Rm, phm = _class_frame(G['R'][m].astype(float), tab[sel[m], 3], tab[sel[m], 4], G['INPL'][0][m], G['INPL'][1][m], G['INPL'][2][m], (sr, s2, rot3), Bd, box, ss, fs)
                                    Pc = np.mean(slice_batch(Gc[0].real.astype(np.float32), Gc[0].imag.astype(np.float32), Rm, Bd.kx, Bd.ky, conv) * phm, 0); A = 0.5 * (AE[G['AROW'][b]] + AO[G['AROW'][b]]).astype(complex); q = Bd.shell >= 2
                                    cc.append(np.sum(Bd.wq[q] * np.real(np.conj(Pc[q]) * A[q])) / np.sqrt(np.sum(Bd.wq[q] * np.abs(Pc[q]) ** 2) * np.sum(Bd.wq[q] * np.abs(A[q]) ** 2)))
                                table.append(((sr, s2, rot3), float(np.mean(cc))))
                                if best is None or np.mean(cc) > best[1]: best = ((sr, s2, rot3), float(np.mean(cc)))
                    G['INPLC'] = best[0]; print('%s   class-frame convention: best (rot sign %+d, 2D-shift sign %+d, 3D shift rotated %+d) ncc %.3f; all: %s' % (T, *best[0], best[1], ' '.join('%+d%+d%+d:%.2f' % (*c, v) for c, v in table)), flush=True)
                BP = {}
                for b in ub:
                    rows_b = G['BAGROWS'][b]; per = []
                    for h in (0, 1):
                        rb = rows_b[eo_all[rows_b] == h]
                        if len(rb) == 0: per.append((np.zeros((0, 3, 3)), np.zeros(0), np.zeros(0))); continue
                        Rb = np.array([euler2m([-tab[i, 0], tab[i, 1], -tab[i, 2]]) for i in rb]); xb, yb = tab[rb, 3].copy(), tab[rb, 4].copy()
                        if G['INPLC'] is not None:
                            sr, s2, rot3 = G['INPLC']; E3a, X2a, Y2a = (v[rb] for v in G['INPL_ALL']); Rz = np.array([rotmat(sr * t, 3) for t in E3a]); Rb = np.einsum('mij,mjk->mik', Rz, Rb)
                            if rot3 != 0: c_, s_ = np.cos(np.radians(rot3 * E3a)), np.sin(np.radians(rot3 * E3a)); xb, yb = c_ * xb - s_ * yb, s_ * xb + c_ * yb
                            xb, yb = xb + s2 * X2a, yb + s2 * Y2a
                        per.append((Rb, xb, yb))
                    BP[b] = per
                G['BAGPOSE'] = BP
                if a.avg_mode == 'unmixc':   # parent-centred: the parent's average = member-weighted mean of its sub-classes' averages; its pose sample = union of theirs
                    par_of = G.get('PAROF', {}); PA = {}
                    for pp in set(par_of.values()):
                        bs = [b for b in ub if par_of.get(b) == pp]
                        if not bs: continue
                        wts = np.array([len(G['BAGROWS'][b]) for b in bs], float); j_ = [G['AROW'][b] for b in bs]
                        Ap = [(np.sum(wts[:, None] * AE[j_], 0) / wts.sum()), (np.sum(wts[:, None] * AO[j_], 0) / wts.sum())]
                        per = []
                        for h in (0, 1):
                            Rs = np.concatenate([BP[b][h][0] for b in bs]); xs = np.concatenate([BP[b][h][1] for b in bs]); ys = np.concatenate([BP[b][h][2] for b in bs])
                            if len(Rs) > a.avg_npose: sel_ = rng.choice(len(Rs), a.avg_npose, replace=False); Rs, xs, ys = Rs[sel_], xs[sel_], ys[sel_]
                            per.append((Rs, xs, ys))
                        PA[pp] = (Ap, per)
                    G['PARAVG'] = PA; print('%s   parent-centred prior: %d parents with averages and pose samples' % (T, len(PA)), flush=True)
                snr_avg = float(np.median([np.sum(np.abs(AE[j]) ** 2) / max(np.sum(ASIG[j][Bd.shell]), 1e-30) for j in range(len(ub))]))
                print('%s   averages sampled at the band: %d x %d coefficients, median power / even-odd noise %.1f' % (T, len(ub), Bd.n, snr_avg), flush=True)
            cur_lp = lp
        # ---- M-step: smoothed responsibilities -> soft cross-half backprojection of every component
        Ksm1, _, _, _ = grid_kernel(a.n1, a.n2, sig, a.n3)
        if a.kernel == 'graph' and it >= a.graph_after and 'REF' in G and len(G['REF'][0]) >= ncell:
            # ---- state graph: similarity = cross-half correlation of the state deviations (noise in one half cannot create a neighbour)
            n3_ = G['REF'][0][0][0].shape[0]; kk_ = np.arange(n3_) - n3_ // 2; rr_ = np.sqrt(sum(q ** 2 for q in np.meshgrid(kk_, kk_, kk_, indexing='ij'))); bm_ = (rr_ <= n3_ // 2) & (rr_ >= 1)
            dv = [np.array([((G['REF'][hh][m_][0] - Gc[hh].real) + 1j * (G['REF'][hh][m_][1] - Gc[hh].imag))[bm_] for m_ in range(ncell)]) for hh in (0, 1)]
            nrm = [np.sqrt(np.sum(np.abs(d) ** 2, 1)) for d in dv]; C01 = np.real(np.conj(dv[0]) @ dv[1].T) / np.maximum(np.outer(nrm[0], nrm[1]), 1e-30); Sg = 0.5 * (C01 + C01.T)
            d2g = np.maximum(2 * (1 - Sg), 0.); np.fill_diagonal(d2g, 0.)
            nn_eff = a.graph_nn * sig / sigmas[-1]   # wide neighbourhood early (as the lattice sigma), graph_nn states at the end
            # per-target bandwidth by bisection so that the effective number of neighbours (1 / sum p^2 over the kernel column) equals nn_eff
            Kg = np.zeros((ncell, ncell)); sgs = np.zeros(ncell)
            for m_ in range(ncell):
                lo_, hi_ = 1e-3, 10.
                for _ in range(40):
                    sg_ = np.sqrt(lo_ * hi_); pk = np.exp(-0.5 * d2g[:, m_] / sg_ ** 2); pk /= pk.sum(); ne = 1. / np.sum(pk ** 2)
                    if ne < nn_eff: lo_ = sg_
                    else: hi_ = sg_
                sg_ = np.sqrt(lo_ * hi_); Kg[:, m_] = np.exp(-0.5 * d2g[:, m_] / sg_ ** 2); sgs[m_] = sg_
            sg = float(np.median(sgs)); Ksm1 = Kg / Kg.sum(0, keepdims=True)
            # display coordinates: 2-D spectral embedding of the graph (normalised Laplacian), scaled onto the board's range
            dg = Kg.sum(1); Lsym = (Kg / np.sqrt(np.outer(dg, dg))); w_, v_ = np.linalg.eigh(Lsym); emb = v_[:, -3:-1][:, ::-1] / np.sqrt(dg)[:, None]
            for q in range(2): emb[:, q] = (emb[:, q] - emb[:, q].min()) / max(emb[:, q].max() - emb[:, q].min(), 1e-12) * ((a.n1 if q == 0 else a.n2) - 1)
            i1 = np.concatenate([emb[:, 0]] * a.nlayer); i2 = np.concatenate([emb[:, 1]] * a.nlayer)
            print('%s   state graph: cross-half self-reproducibility median %.2f; nearest-neighbour similarity median %.2f; bandwidth %.2f (nn_eff %.1f); effective neighbours per state median %.1f' % (T, np.median(np.diag(C01)), np.median(np.sort(Sg + np.eye(ncell) * -9, 1)[:, -1]), sg, nn_eff, np.median(1. / np.sum(Ksm1 ** 2, 0))), flush=True)
        if a.adapt > 0 and it > 0 and 'REF' in G and len(G['REF'][0]) >= ncell:   # per-cell width from the cross-half reproducibility of the cell deviation
            n3_ = G['REF'][0][0][0].shape[0]; kk_ = np.arange(n3_) - n3_ // 2; rr_ = np.sqrt(sum(q ** 2 for q in np.meshgrid(kk_, kk_, kk_, indexing='ij'))); bm_ = (rr_ <= n3_ // 2) & (rr_ >= 1)
            crep = np.zeros(ncell)
            for m_ in range(ncell):
                dv = [((G['REF'][hh][m_][0] - Gc[hh].real) + 1j * (G['REF'][hh][m_][1] - Gc[hh].imag))[bm_] for hh in (0, 1)]
                crep[m_] = float(np.real(np.sum(np.conj(dv[0]) * dv[1])) / max(np.sqrt(np.sum(np.abs(dv[0]) ** 2) * np.sum(np.abs(dv[1]) ** 2)), 1e-30))
            fac = 1 + a.adapt * (1 - 2 * np.clip((crep - 0.3) / 0.5, 0, 1)); sig_m = sig * fac
            d2_ = (i1c[:, None] - i1c[None, :]) ** 2 + (i2c[:, None] - i2c[None, :]) ** 2 + (i3c[:, None] - i3c[None, :]) ** 2
            Kc = np.exp(-0.5 * d2_ / np.maximum(sig_m[None, :], 1e-6) ** 2); Ksm1 = Kc / Kc.sum(0, keepdims=True)   # column m = the kernel of target cell m, width sigma_m
            print('%s   adaptive kernel: cell-deviation cross-half corr median %.2f (min %.2f, max %.2f); sigma range %.2f-%.2f cells' % (T, np.median(crep), crep.min(), crep.max(), sig_m.min(), sig_m.max()), flush=True)
        Ksm = np.kron(np.eye(a.nlayer), Ksm1); W = np.c_[r[:, :ngrid] @ Ksm, r[:, ngrid:]] * pw[:, None]
        num = den = None
        with ctx.Pool(a.nproc) as pool:
            if FLIB is not None and G.get('LLS') is not None:   # pose-marginalised M-step: every particle at its neighbourhood poses, weights = joint pose-state posterior (kernel-smoothed over states)
                lls_ = G['LLS']; mx_ = lls_.max(axis=(0, 2), keepdims=True); jp = np.exp(lls_ - mx_); jp /= jp.sum(axis=(0, 2), keepdims=True)   # (npose, n, K) joint posterior
                if G.get('pfac') is not None and a.pnoise_m: jp = jp / G['pfac'][None, :, None]
                R0 = G['R'].copy(); F0 = G['F']; num = den = None
                for jpose, (d_, dx_, dy_) in enumerate(G['PM']):
                    G['R'] = np.ascontiguousarray(np.einsum('mij,jk->mik', R0.astype(float), d_), R0.dtype)
                    if dx_ != 0. or dy_ != 0.: G['F'] = (F0 * np.exp(-2j * np.pi * (Bd.kx[None, :] * (ss * fs * dx_) + Bd.ky[None, :] * (ss * fs * dy_)) / box)).astype(F0.dtype)
                    Wj = np.c_[jp[jpose][:, :ngrid] @ Ksm, jp[jpose][:, ngrid:]] * pw[:, None]; nj, dj = fort_mstep(FLIB, chunks, lambda s0, s1: Wj[s0:s1], Wj.shape[1]); num = nj if num is None else num + nj; den = dj if den is None else den + dj
                    G['F'] = F0
                G['R'] = R0
            elif FLIB is not None:
                Wm = W if not (a.pnoise_m and G.get('pfac') is not None) else W / G['pfac'][:, None]
                num, den = fort_mstep(FLIB, chunks, lambda s0, s1: Wm[s0:s1], Wm.shape[1])
                if os.environ.get('KLEM_CHECK', '') and it == 0:   # correctness gate: numpy vs Fortran M-step on the first chunk
                    s0_, s1_ = chunks[0]; n0, d0 = _mstep_chunk((s0_, s1_, W[s0_:s1_])); n1, d1 = fort_mstep(FLIB, [chunks[0]], lambda s0, s1: W[s0:s1], W.shape[1])
                    print('%s   KLEM_CHECK M-step: max |numpy - fortran| / max: num %.2e, den %.2e' % (T, np.abs(n0 - n1).max() / max(np.abs(n0).max(), 1e-30), np.abs(d0 - d1).max() / max(np.abs(d0).max(), 1e-30)), flush=True)
            else:
                for p in pool.imap_unordered(_mstep_chunk, [(s0, s1, W[s0:s1]) for s0, s1 in chunks]): num = p[0] if num is None else num + p[0]; den = p[1] if den is None else den + p[1]
        if a.tie_layers and a.nlayer == 2:
            # occupancy model: base map of cell m from ALL its particles (layer-2 particles with the component removed), one shared component Delta
            Dl = []
            for hh in (0, 1):
                V1 = num[hh, :ncell] / (den[hh, :ncell] + eps); V2 = num[hh, ncell:2 * ncell] / (den[hh, ncell:2 * ncell] + eps)
                wgt = np.minimum(den[hh, :ncell], den[hh, ncell:2 * ncell]).sum(1); wgt = wgt / max(wgt.sum(), 1e-30)
                Dl.append(np.sum(wgt[:, None] * (V2 - V1), 0))
            if a.tie_support > 0 or a.tie_support_mask:   # positive, compact component: clip negatives in real space; support = reproducible positive part (largest components) or a given site mask
                ifft_t = lambda X: np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(X)))); fft_t = lambda v: np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v)))
                from scipy.ndimage import gaussian_filter as _gf, label as _lab, binary_dilation as _bd
                dr = [ifft_t(Dl[hh].reshape(n3, n3, n3)) for hh in (0, 1)]; S_ = np.ones((n3, n3, n3))
                if a.tie_support_mask:
                    if G.get('TMASK_SRC') is None: G['TMASK_SRC'] = read_mrc(a.tie_support_mask)
                    mk_ = G['TMASK_SRC']; nb_ = mk_.shape[0]; Fm_ = fft_t(mk_); c_ = nb_ // 2; h_ = n3 // 2
                    S_ = np.clip(ifft_t(Fm_[c_ - h_:c_ - h_ + n3, c_ - h_:c_ - h_ + n3, c_ - h_:c_ - h_ + n3]) * (n3 / nb_) ** 3, 0, 1)
                elif it >= a.tie_support_after:
                    p_ = _gf(np.maximum(dr[0], 0) * np.maximum(dr[1], 0), 1.0); m_ = p_ > a.tie_support * p_.max(); lb_, nl_ = _lab(m_)
                    if nl_ > 0:
                        ms_ = np.array([p_[lb_ == i].sum() for i in range(1, nl_ + 1)]); keep_ = [i + 1 for i in np.argsort(-ms_)[:a.tie_ncomp]]; m_ = np.isin(lb_, keep_)
                    S_ = np.clip(_gf(_bd(m_, iterations=1).astype(float), 1.0), 0, 1)
                ds = [S_ * np.maximum(dr[hh], 0) for hh in (0, 1)]; Dl = [fft_t(ds[hh]).ravel() for hh in (0, 1)]
                cc_t = float(np.sum(ds[0] * ds[1]) / max(np.sqrt(np.sum(ds[0] ** 2) * np.sum(ds[1] ** 2)), 1e-30)); zz_t, yy_t, xx_t = np.indices((n3, n3, n3)); rad_t = np.hypot(np.hypot(xx_t - n3 / 2, yy_t - n3 / 2), zz_t - n3 / 2) * (box * SMPD / n3)
                print('%s   tied component: support %.1f%% of the coarse voxels (%s), centroid radius %.0f A, cross-half corr of the positive component %.2f' % (T, 100 * np.mean(S_ > 0.5), 'given mask' if a.tie_support_mask else ('reproducible positive part' if it >= a.tie_support_after else 'unrestricted'), float(np.sum(rad_t * S_) / max(S_.sum(), 1e-30)), cc_t), flush=True)
            for hh in (0, 1):
                Delta = Dl[hh]
                base = (num[hh, :ncell] + num[hh, ncell:2 * ncell] - den[hh, ncell:2 * ncell] * Delta[None, :]) / (den[hh, :ncell] + den[hh, ncell:2 * ncell] + eps)
                num[hh, :ncell] = base * (den[hh, :ncell] + eps); den[hh, :ncell] = den[hh, :ncell] + eps - eps   # so that the ratio below reproduces base
                num[hh, :ncell] = base * (den[hh, :ncell] + eps) - eps * Gc[hh].ravel()[None, :]
                num[hh, ncell:2 * ncell] = (base + Delta[None, :]) * (den[hh, ncell:2 * ncell] + eps) - eps * Gc[hh].ravel()[None, :]
            if it % 3 == 2: print('%s   tied layers: component power / base power (shell-summed) %.3f' % (T, float(np.sum(np.abs(Delta) ** 2) / max(np.sum(np.abs(base) ** 2) / ncell, 1e-30))), flush=True)
        REF = [[], []]
        if it == 0 and a.nfree > 0:   # junk-like initial free components: a 35 A blob and a half-contrast consensus (no particle responsibilities yet)
            kk = np.arange(-kmax, kmax + 1); rr = np.sqrt(sum(q ** 2 for q in np.meshgrid(kk, kk, kk, indexing='ij'))); lp35 = (rr <= box * SMPD / 35.).ravel()
            for hh in (0, 1):
                num[hh, ngrid] = Gc[hh].ravel() * lp35 * (den[hh, ngrid].max() + eps); den[hh, ngrid] = den[hh, ngrid].max() + eps
                if a.nfree > 1: num[hh, ngrid + 1] = 0.5 * Gc[hh].ravel() * (den[hh, ngrid + 1].max() + eps); den[hh, ngrid + 1] = den[hh, ngrid + 1].max() + eps
        FB = flex_basis(a.flex_prior, n3, G) if a.flex_prior else None
        for hh in (0, 1):
            for k in range(K):
                tgt = Gc[hh].ravel()
                if FB is not None and k < ngrid:   # PPCA prior: the target is the consensus + the cell's own projection onto the flex subspace
                    d0 = (num[hh, k] + eps * Gc[hh].ravel()) / (den[hh, k] + eps) - Gc[hh].ravel(); tgt = Gc[hh].ravel() + FB @ (FB.conj().T @ d0)
                Gk = ((num[hh, k] + eps * tgt) / (den[hh, k] + eps)).reshape(n3, n3, n3)
                REF[hh].append((np.ascontiguousarray(Gk.real, np.float32), np.ascontiguousarray(Gk.imag, np.float32)))
        if a.node_basis > 0 and it + 1 >= a.node_basis_after:   # ---- node basis: PPCA-EM M-step on the per-image residuals (own half), cross-half validated
            q = a.node_basis; n3_ = REF[0][0][0].shape[0]; nv_ = n3_ ** 3; G['NB_MARG'] = a.node_basis_marg; init_now = False
            cells = [[(REF[hh][k][0] + 1j * REF[hh][k][1]) for k in range(ngrid)] for hh in (0, 1)]
            devs = [[c - Gc[hh] for c in cells[hh]] for hh in (0, 1)]
            if G.get('NB') is None or G['NB'][0][0][0].shape[0] != n3_ or G.get('ZPOST') is None:   # (re)initialise from the mass-weighted PCA of this half's cell deviations (no deflation, no evidence term until the first residual fit)
                NBn = []; init_now = True
                for hh in (0, 1):
                    Dm = np.array([d.ravel() for d in devs[hh]]); wm = np.maximum(r[:, :ngrid].sum(0), 1e-6); wm = wm / wm.sum(); mu = (wm[:, None] * Dm).sum(0); Dm = Dm - mu
                    Gm = np.real((np.sqrt(wm)[:, None] * Dm) @ np.conj(np.sqrt(wm)[:, None] * Dm).T); ev, U = np.linalg.eigh(Gm); ev, U = ev[::-1], U[:, ::-1]
                    comps = []
                    for j in range(q):
                        bj = (U[:, j] * np.sqrt(wm)) @ Dm; bj = bj / max(np.linalg.norm(bj), 1e-30) * np.sqrt(max(ev[j], 0.)); comps.append(bj.reshape(n3_, n3_, n3_))
                    NBn.append(comps)
                print('%s   node basis (q %d): initialised from the cell PCA (relative eigenvalues %s)' % (T, q, np.array2string(np.maximum(ev[:q], 0) / max(ev[0], 1e-30), precision=2)), flush=True)
            else:
                ZPn = G['ZPOST'][:, :ngrid, :]; W1 = np.einsum('ik,ikj->ikj', r[:, :ngrid], ZPn).reshape(n, ngrid * q); W2 = np.einsum('ik,ikj->ij', r[:, :ngrid], ZPn ** 2); Wb = np.ascontiguousarray(np.c_[W1, W2], np.float32)
                if FLIB is not None: nb_, db_ = fort_mstep(FLIB, chunks, lambda s0, s1: Wb[s0:s1], Wb.shape[1])
                else:
                    nb_ = db_ = None
                    with ctx.Pool(a.nproc) as pool:
                        for p in pool.imap_unordered(_mstep_chunk, [(s0, s1, Wb[s0:s1]) for s0, s1 in chunks]): nb_ = p[0] if nb_ is None else nb_ + p[0]; db_ = p[1] if db_ is None else db_ + p[1]
                NBn = []
                for hh in (0, 1):
                    comps = []
                    for j in range(q):
                        numB = np.zeros(nv_, complex)
                        for k in range(ngrid): numB += nb_[hh, k * q + j] - db_[hh, k * q + j] * cells[hh][k].ravel()
                        comps.append((numB / (db_[hh, ngrid * q + j] + eps)).reshape(n3_, n3_, n3_))
                    NBn.append(comps)
            if a.node_basis_deflate and not init_now:   # remove the span of this half's cell deviations: the basis models within-cell residual motion only
                for hh in (0, 1):
                    Uo = []
                    for d in devs[hh]:
                        v = d.ravel().copy()
                        for u in Uo: v -= np.vdot(u, v) * u
                        nv = np.linalg.norm(v)
                        if nv > 1e-6 * max(np.linalg.norm(d), 1e-30): Uo.append(v / nv)
                    for j in range(q):
                        v = NBn[hh][j].ravel().copy()
                        for u in Uo: v -= np.vdot(u, v) * u
                        NBn[hh][j] = v.reshape(n3_, n3_, n3_)
            for hh in (0, 1):
                for j in range(1, q):
                    for l in range(j): NBn[hh][j] = NBn[hh][j] - np.vdot(NBn[hh][l].ravel(), NBn[hh][j].ravel()) / max(np.vdot(NBn[hh][l].ravel(), NBn[hh][l].ravel()).real, 1e-30) * NBn[hh][l]
            for j in range(q):   # the halves' components are defined up to a sign: align them
                if np.real(np.vdot(NBn[0][j].ravel(), NBn[1][j].ravel())) < 0: NBn[1][j] = -NBn[1][j]
            reps = [np.real(np.vdot(NBn[0][j].ravel(), NBn[1][j].ravel())) / max(np.linalg.norm(NBn[0][j]) * np.linalg.norm(NBn[1][j]), 1e-30) for j in range(q)]
            amp = [0.5 * (np.linalg.norm(NBn[0][j]) + np.linalg.norm(NBn[1][j])) / max(np.sqrt(np.mean([np.linalg.norm(d) ** 2 for d in devs[0]])), 1e-30) for j in range(q)]
            for j in range(q):
                if reps[j] < a.node_basis_rep: NBn[0][j][:] = 0.; NBn[1][j][:] = 0.
            G['NB'] = [[(np.ascontiguousarray(c.real, np.float32), np.ascontiguousarray(c.imag, np.float32)) for c in NBn[hh]] for hh in (0, 1)]; G['NB_INIT'] = init_now
            print('%s   node basis: cross-half reproducibility per component %s (kept if >= %.2f); amplitude relative to the RMS cell deviation %s' % (T, np.array2string(np.array(reps), precision=2), a.node_basis_rep, np.array2string(np.array(amp), precision=2)), flush=True)
        if (a.devw > 0 or a.focus > 0 or a.focus_mask) and it > 0:
            n3_ = REF[0][0][0].shape[0]; kk_ = np.arange(n3_) - n3_ // 2; rr_ = np.sqrt(sum(q ** 2 for q in np.meshgrid(kk_, kk_, kk_, indexing='ij'))); sh_ = np.rint(rr_).astype(int).ravel(); nsh_ = sh_.max() + 1
            D0 = np.array([(REF[0][k_][0] - Gc[0].real) + 1j * (REF[0][k_][1] - Gc[0].imag) for k_ in range(ncell)]); D1 = np.array([(REF[1][k_][0] - Gc[1].real) + 1j * (REF[1][k_][1] - Gc[1].imag) for k_ in range(ncell)])
            if a.devw > 0:   # pooled cross-half FSC of the cell deviations per shell -> likelihood shell weights 2F/(1+F)
                num_ = np.bincount(sh_, np.sum(np.real(np.conj(D0) * D1), 0).ravel(), minlength=nsh_); d0_ = np.bincount(sh_, np.sum(np.abs(D0) ** 2, 0).ravel(), minlength=nsh_); d1_ = np.bincount(sh_, np.sum(np.abs(D1) ** 2, 0).ravel(), minlength=nsh_)
                fsc_ = num_ / np.maximum(np.sqrt(d0_ * d1_), 1e-30); wsh = np.clip(2 * fsc_ / (1 + np.maximum(fsc_, 1e-6)), a.devw, 1.); wsh[:2] = 1.
                G['sig2'] = G['sig2_base'] / wsh[np.minimum(Bd.shell, nsh_ - 1)]
                print('%s   deviation FSC per shell %s -> likelihood weights %s' % (T, np.array2string(fsc_[1:kmax + 1:max(1, kmax // 6)], precision=2), np.array2string(wsh[1:kmax + 1:max(1, kmax // 6)], precision=2)), flush=True)
            if a.focus > 0 or a.focus_mask:   # focus mask in real space: data-derived (cross-half product of the deviations) or a given focal mask, applied to the E-step deviations
                ifft = lambda X: np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(X)))); fft_ = lambda v: np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v)))
                from scipy.ndimage import gaussian_filter
                if a.focus_mask:   # Fourier-crop the full-box mask to the coarse grid of this band
                    if G.get('FMASK_SRC') is None: G['FMASK_SRC'] = read_mrc(a.focus_mask)
                    mk_ = G['FMASK_SRC']; nb_ = mk_.shape[0]; Fm_ = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(mk_))); c_ = nb_ // 2; h_ = n3_ // 2
                    fm = np.clip(np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(Fm_[c_ - h_:c_ - h_ + n3_, c_ - h_:c_ - h_ + n3_, c_ - h_:c_ - h_ + n3_])))) * (n3_ / nb_) ** 3, 0, 1)
                else:
                    v_ = np.zeros((n3_, n3_, n3_))
                    for k_ in range(ncell): v_ += ifft(D0[k_]) * ifft(D1[k_])
                    v_ = gaussian_filter(np.maximum(v_, 0), 1.5); fm = np.clip(v_ / max(np.percentile(v_, 99) * a.focus, 1e-30), 0, 1) ** 0.5
                for hh in (0, 1):
                    for k_ in range(K):
                        Dk = (REF[hh][k_][0] - Gc[hh].real) + 1j * (REF[hh][k_][1] - Gc[hh].imag); Fk = Gc[hh] + fft_(fm * ifft(Dk))
                        REF[hh][k_] = (np.ascontiguousarray(Fk.real, np.float32), np.ascontiguousarray(Fk.imag, np.float32))
                print('%s   focus mask: %.1f%% of the coarse voxels above 0.5 (%s)' % (T, 100 * np.mean(fm > 0.5), 'focal mask ' + a.focus_mask if a.focus_mask else 'reproducible variability of the cells'), flush=True)
        if a.deconv > 0:   # one-step deconvolution of the grid kernel on the cell volumes: V + deconv x (V - K V), per layer (the consensus part passes through)
            for hh in (0, 1):
                for lay in range(a.nlayer):
                    Vr = np.array([REF[hh][lay * ncell + k][0] for k in range(ncell)]); Vi = np.array([REF[hh][lay * ncell + k][1] for k in range(ncell)])
                    Vr = Vr + a.deconv * (Vr - np.tensordot(Ksm1, Vr, 1)); Vi = Vi + a.deconv * (Vi - np.tensordot(Ksm1, Vi, 1))
                    for k in range(ncell): REF[hh][lay * ncell + k] = (np.ascontiguousarray(Vr[k], np.float32), np.ascontiguousarray(Vi[k], np.float32))
        if split_at >= 0 and it >= split_at:   # reproducibility of the layer difference per cell: cross-half correlation (band = the whole coarse grid)
            n3_ = REF[0][0][0].shape[0]; kk_ = np.arange(n3_) - n3_ // 2; rr_ = np.sqrt(sum(q ** 2 for q in np.meshgrid(kk_, kk_, kk_, indexing='ij'))); bm_ = (rr_ <= n3_ // 2) & (rr_ >= 1)
            rep_, pw_ = [], []
            for k in range(ncell):
                dd = [((REF[hh][ncell + k][0] - REF[hh][k][0]) + 1j * (REF[hh][ncell + k][1] - REF[hh][k][1]))[bm_] for hh in (0, 1)]; v0 = (REF[0][k][0] + 1j * REF[0][k][1])[bm_]
                rep_.append(float(np.real(np.sum(np.conj(dd[0]) * dd[1])) / max(np.sqrt(np.sum(np.abs(dd[0]) ** 2) * np.sum(np.abs(dd[1]) ** 2)), 1e-30))); pw_.append(float(np.mean(np.abs(dd[0]) ** 2) / max(np.mean(np.abs(v0) ** 2), 1e-30)))
            rep_, pw_ = np.array(rep_), np.array(pw_); kb = int(np.argmax(rep_))
            print('%s   split reproducibility: cross-half corr of the layer difference per cell median %.2f, max %.2f at cell (%d,%d) (difference power / cell power there %.3f, median %.3f); cells above 0.5: %d of %d' % (T, np.median(rep_), rep_[kb], i1c[kb], i2c[kb], pw_[kb], np.median(pw_), int((rep_ > 0.5).sum()), ncell), flush=True)
        G['REF'] = REF
        # ---- E-step
        if avg_mode:
            if FLIB is not None: ll = fort_estep(FLIB, chunks)
            else:
                with ctx.Pool(a.nproc) as pool: parts = pool.map(_estep_chunk, chunks)
                ll = np.concatenate(parts); del parts
            ubl = list(G['BAGPOSE'].keys()); nb_ = max(1, len(ubl) // (4 * a.nproc)); groups = [ubl[j:j + nb_] for j in range(0, len(ubl), nb_)]
            with ctx.Pool(a.nproc) as pool: parts = pool.map(_estep_bag_chunk, groups)
            bagll = {}
            for q in parts: bagll.update(q)
            del parts
            # members inherit their bag's half-specific log-likelihood (a half without members in the bag falls back to the other half)
            llp = ll.copy(); nb = 0
            for b, res in bagll.items():
                for h in (0, 1):
                    m = (bag == b) & (eo == h); v = res[h] if res[h] is not None else res[1 - h]
                    if not (m.any() and v is not None): continue
                    if a.avg_mode in ('unmix', 'unmixc'): ll[m] = llp[m] + a.avg_prior * v[None, :]
                    else: ll[m] = (1 - a.avg_mix) * v[None, :] + a.avg_mix * llp[m]
                    nb += int(m.sum())
            Lb = ll - ll.max(1, keepdims=True); rb = np.exp(Lb); rb /= rb.sum(1, keepdims=True)
            print('%s   sub-class E-step (%s): %d particles got their sub-class %s (%d bags), mean max r %.2f' % (T, a.avg_mode, nb, 'mixture prior' if a.avg_mode == 'unmix' else 'responsibilities', len(bagll), rb.max(1).mean()), flush=True)
        elif a.pose_marg > 0:   # marginalise over a pose neighbourhood (and shifts), logsumexp; Fortran or numpy
            from scipy.spatial.transform import Rotation as Rot
            if a.pose_marg_n >= 27: gv = np.array([[i, j, k_] for i in (-1, 0, 1) for j in (-1, 0, 1) for k_ in (-1, 0, 1)], float)
            else: gv = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float)
            dR = Rot.from_rotvec(gv * np.radians(a.pose_marg)).as_matrix(); shifts = [(0., 0.)] + ([(a.shift_marg, 0.), (-a.shift_marg, 0.), (0., a.shift_marg), (0., -a.shift_marg)] if a.shift_marg > 0 else [])
            R0 = G['R'].copy(); F0 = G['F']; lls = []; PM = []
            for d_ in dR:
                G['R'] = np.ascontiguousarray(np.einsum('mij,jk->mik', R0.astype(float), d_), R0.dtype)
                for (dx_, dy_) in shifts:
                    if dx_ != 0. or dy_ != 0.: G['F'] = (F0 * np.exp(-2j * np.pi * (Bd.kx[None, :] * (ss * fs * dx_) + Bd.ky[None, :] * (ss * fs * dy_)) / box)).astype(F0.dtype)
                    if FLIB is not None: lls.append(fort_estep(FLIB, chunks))
                    else:
                        with ctx.Pool(a.nproc) as pool: parts = pool.map(_estep_chunk, chunks)
                        lls.append(np.concatenate(parts)); del parts
                    PM.append((d_, dx_, dy_)); G['F'] = F0
            G['R'] = R0; lls = np.array(lls); mx = lls.max(0); ll = mx + np.log(np.mean(np.exp(lls - mx[None]), 0))
            if a.pose_marg_mstep: G['LLS'] = lls; G['PM'] = PM
            else: G.pop('LLS', None)
        elif FLIB is not None and G.get('NB') is None:
            ll = fort_estep(FLIB, chunks)
            if os.environ.get('KLEM_CHECK', '') and it == 0:   # correctness gate: numpy vs Fortran on the first chunk
                ll0 = _estep_chunk(chunks[0]); d_ = np.abs(ll0 - ll[chunks[0][0]:chunks[0][1]]).max() / max(np.abs(ll0).max(), 1e-30); print('%s   KLEM_CHECK E-step: max |numpy - fortran| / max = %.2e' % (T, d_), flush=True)
        else:
            with ctx.Pool(a.nproc) as pool: parts = pool.map(_estep_chunk, chunks)
            if G.get('NB') is not None: ll = np.concatenate([p[0] for p in parts]); G['ZPOST'] = np.concatenate([p[1] for p in parts])
            else: ll = np.concatenate(parts); G.pop('ZPOST', None)
            del parts
        if bag is not None and a.bag_beta > 0 and not avg_mode:   # bag prior: add bag_beta x the bag-mean log-likelihood (members share evidence, keep their own)
            bu, bi = np.unique(bag, return_inverse=True); cnt = np.bincount(bi); S_ = np.zeros((len(bu), ll.shape[1]))
            np.add.at(S_, bi, ll); mean_bag = S_ / cnt[:, None]; prior = mean_bag[bi]
            if par is not None:   # centred on the parent class: a view's own preference among the cells cancels, only the within-parent part is shared
                pu, pi = np.unique(par, return_inverse=True); pc = np.bincount(pi); P_ = np.zeros((len(pu), ll.shape[1])); np.add.at(P_, pi, ll); prior = prior - (P_ / pc[:, None])[pi]; prior[par <= 0] = 0.
            prior = a.bag_beta * prior; prior[bag <= 0] = 0.; ll = ll + prior
        if G.get('pfac') is not None: ll = ll / G['pfac'][:, None]
        if a.mixprop and 'mixp' in G: ll = ll + np.log(G['mixp'])[None, :]
        if a.temp_final < 1 and it >= len(lps) - 3: ll = ll / a.temp_final
        if a.temp_start > 1:
            nt_ = max(1, int(round(2 * len(lps) / 3))); tcur = a.temp_start + (1. - a.temp_start) * min(it, nt_) / nt_; ll = ll / tcur
        L = ll - ll.max(1, keepdims=True); r = np.exp(L); r /= r.sum(1, keepdims=True)
        if a.balance > 0 and (a.balance_until <= 0 or it < a.balance_until):   # ---- OT-balanced E-step (Sinkhorn): grid columns -> equal shares of the grid mass, rows -> 1; free components unconstrained
            from scipy.special import logsumexp
            gmass = float(r[:, :ngrid].sum()); tgt = np.log(max(gmass, 1e-30) / ngrid); u = np.zeros(n); v = np.zeros(K); Lb = L.astype(float)
            if a.nlayer > 1:   # layer-aware: cells balanced within each layer, the layer split itself stays free (the previous round's layer masses are the targets)
                lmass = np.array([float(r[:, :ngrid][:, layer_of == l].sum()) for l in range(a.nlayer)]); tgt = np.log(np.maximum(lmass[layer_of], 1e-30) / ncell)
                print('%s   balanced E-step: layer masses %s (targets per cell within each layer)' % (T, np.array2string(lmass / max(gmass, 1e-30), precision=3)), flush=True)
            for _ in range(30):
                v[:ngrid] = tgt - logsumexp(Lb[:, :ngrid] + u[:, None], axis=0); u = -logsumexp(Lb + v[None, :], axis=1)
            rb = np.exp(Lb + u[:, None] + a.balance * v[None, :]); rb /= rb.sum(1, keepdims=True); occ0 = r[:, :ngrid].sum(0) / n; occ1 = rb[:, :ngrid].sum(0) / n
            print('%s   balanced E-step (lambda %.2f): grid occupancy min/max %.3f/%.3f -> %.3f/%.3f, free %s -> %s; mean total variation |r - r_bal| %.3f; argmax changed for %.2f of particles' % (T, a.balance, occ0.min(), occ0.max(), occ1.min(), occ1.max(), np.array2string(r[:, ngrid:].sum(0) / n, precision=3), np.array2string(rb[:, ngrid:].sum(0) / n, precision=3), float(np.abs(rb - r).sum(1).mean() / 2), float(np.mean(rb.argmax(1) != r.argmax(1)))), flush=True)
            r = rb
        if G.get('ZPOST') is not None:   # node coordinate per image: responsibility-weighted posterior mean over the cells
            G['zb'] = np.einsum('ik,ikj->ij', r, G['ZPOST']); zb_ = G['zb']; gq_ = (conf > 0)
            try:
                from scipy.stats import spearmanr as _sp; cl_ = np.array([c in closed for c in conf]); ax_ = np.array([axis.get(c, 0.) for c in conf])
                aucs = [abs(0.5 - np.argsort(np.argsort(zb_[gq_, j]))[cl_[gq_]].mean() / gq_.sum()) * 2 for j in range(zb_.shape[1])]
                print('%s   node basis coordinate (per image): Spearman(zb_j, truth axis) %s | closed-vs-rest |AUC-0.5|x2 %s | corr with lattice z1 %+.2f, z2 %+.2f' % (T, np.array2string(np.array([_sp(zb_[gq_, j], ax_[gq_])[0] for j in range(zb_.shape[1])]), precision=2), np.array2string(np.array(aucs), precision=2), np.corrcoef(zb_[:, 0], r[:, :ngrid] @ i1)[0, 1], np.corrcoef(zb_[:, 0], r[:, :ngrid] @ i2)[0, 1]), flush=True)
            except Exception as e_: print('%s   node basis readout failed: %s' % (T, e_))
        if a.mixprop: G['mixp'] = np.maximum(r.mean(0), 1e-4 / r.shape[1]); G['mixp'] /= G['mixp'].sum()
        if a.pose_every > 0 and (it + 1) % a.pose_every == 0:
            am_now = r.argmax(1); tp = time.time(); REF_keep = None
            if a.pose_mask or a.pose_ref == 'consensus':   # node alignment policy: the pose step against masked / consensus references (core vs region vs everything)
                RF_ = G['REF']; n3_ = RF_[0][0][0].shape[0]; K_ = len(RF_[0]); ifft_ = lambda X: np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(X)))); fft__ = lambda v: np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v)))
                pm = np.ones((n3_, n3_, n3_))
                if a.pose_mask in ('core', 'region'):
                    from scipy.ndimage import gaussian_filter as gf_
                    v_ = np.zeros((n3_, n3_, n3_))
                    for k_ in range(ncell): v_ += ifft_((RF_[0][k_][0] - Gc[0].real) + 1j * (RF_[0][k_][1] - Gc[0].imag)) * ifft_((RF_[1][k_][0] - Gc[1].real) + 1j * (RF_[1][k_][1] - Gc[1].imag))
                    v_ = gf_(np.maximum(v_, 0), 1.5); fm_ = np.clip(v_ / max(np.percentile(v_, 99) * (a.focus if a.focus > 0 else 0.5), 1e-30), 0, 1) ** 0.5
                    d_ = gf_(ifft_(0.5 * (Gc[0] + Gc[1])), 1.0); d_ = np.clip((d_ - d_.mean()) / max(d_.std(), 1e-30) - 1.0, 0, 1)
                    pm = d_ * (1 - fm_) if a.pose_mask == 'core' else fm_
                elif a.pose_mask:
                    if G.get('PMASK_SRC') is None: G['PMASK_SRC'] = read_mrc(a.pose_mask)
                    mk_ = G['PMASK_SRC']; nb_ = mk_.shape[0]; Fm_ = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(mk_))); c_ = nb_ // 2; h_ = n3_ // 2
                    pm = np.clip(np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(Fm_[c_ - h_:c_ - h_ + n3_, c_ - h_:c_ - h_ + n3_, c_ - h_:c_ - h_ + n3_])))) * (n3_ / nb_) ** 3, 0, 1)
                REF_keep = RF_; REFm = [[None] * K_ for _ in (0, 1)]
                for hh in (0, 1):
                    for k_ in range(K_):
                        src_ = Gc[hh] if a.pose_ref == 'consensus' else (RF_[hh][k_][0] + 1j * RF_[hh][k_][1]); Fk_ = fft__(pm * ifft_(src_))
                        REFm[hh][k_] = (np.ascontiguousarray(Fk_.real, np.float32), np.ascontiguousarray(Fk_.imag, np.float32))
                G['REF'] = REFm
                print('%s   pose step references: %s x mask %s (%.1f%% of the coarse voxels above 0.5)' % (T, a.pose_ref, a.pose_mask or 'none', 100 * np.mean(pm > 0.5)), flush=True)
            with ctx.Pool(a.nproc) as pool: parts = pool.map(_pose_chunk, [(s0, s1, am_now[s0:s1], a.pose_step, a.pose_prior) for s0, s1 in chunks])
            if REF_keep is not None: G['REF'] = REF_keep
            Rn = np.concatenate([q[0] for q in parts]); gain = np.concatenate([q[1] for q in parts]); del parts
            # rotation of the particle frame: R = euler2m([-e1, e2, -e3]); back to SIMPLE Euler angles for the table
            moved = np.mean(np.linalg.norm(Rn - G['R'], axis=(1, 2)) > 1e-6); G['R'] = Rn
            for j, i in enumerate(sel):
                e = m2euler_simple(Rn[j].astype(float)); tab[i, 0], tab[i, 1], tab[i, 2] = (-e[0]) % 360., e[1], (-e[2]) % 360.
            G['tab'] = tab; G['D'].tab = tab
            with open(a.out + '/ptcl3D_pose_it%d.txt' % (it + 1), 'w') as fo_:   # the moved poses (SIMPLE euler convention), one row per EM particle
                for j, i in enumerate(sel): fo_.write('e1=%.4f e2=%.4f e3=%.4f x=%.3f y=%.3f state=1 pind=%d\n' % (tab[i, 0], tab[i, 1], tab[i, 2], tab[i, 3], tab[i, 4], i + 1))
            msg = '%s   pose step (it %d, step %.1f deg): %.2f of particles moved, mean log-likelihood gain %.2f (%.0f s)' % (T, it + 1, a.pose_step, moved, gain.mean(), time.time() - tp)
            if a.pose_truth:
                tt = read_table(a.pose_truth, ('e1', 'e2', 'e3')); ang = []
                for j, i in enumerate(sel[:20000]):
                    Rt = euler2m([-tt[i, 0], tt[i, 1], -tt[i, 2]]); ang.append(np.degrees(np.arccos(np.clip((np.trace(Rn[j].astype(float) @ Rt.T) - 1) / 2, -1, 1))))
                msg += ' | pose error vs truth: median %.2f deg (p90 %.2f) on %d particles' % (np.median(ang), np.percentile(ang, 90), len(ang))
            print(msg, flush=True)
        if split_at >= 0 and it >= split_at:
            l2 = r[:, ncell:ngrid].sum(1) / np.maximum(r[:, :ngrid].sum(1), 1e-12); cl_ = [c_ for c_ in np.unique(conf) if c_ > 0][:14]
            print('%s   layer-2 share: all %.2f | per conformer %s | 10 pct of particles most layer-polarised: %.2f of them beyond 0.9 or below 0.1' % (T, l2.mean(), ' '.join('%d:%.2f' % (c_, l2[conf == c_].mean()) for c_ in cl_), np.mean(np.abs(l2 - 0.5) > 0.4)), flush=True)
        occ = r.sum(0) / n
        print('%s it %d (lp %g, sigma %.1f): occupancy grid min/median/max %.3f/%.3f/%.3f, free %s (%.0f s)' % (T, it + 1, lp, sig, occ[:ngrid].min(), np.median(occ[:ngrid]), occ[:ngrid].max(), np.array2string(occ[ngrid:], precision=3), time.time() - t0), flush=True)
        res = readout('%s it %d' % (T, it + 1), r, i1, i2, ngrid, conf, closed, axis, rng, layer_of=layer_of, i3=i3); hist.append({k: v for k, v in res.items() if k not in ('z1', 'z2', 'z3', 'am')} | dict(it=it + 1, lp=lp, sigma=sig))
    # ---- optional final E-step over all particles (subset EM -> latent for everyone)
    if a.final_all and len(sel_full) > len(sel):   # subset EM (by --nsub or --subset) -> one E-step over ALL particles against this board's references
        sel = sel_full; n = len(sel); eo = eo_all[sel]; G['sel'] = sel; chunks = [(s0, min(s0 + a.chunk, n)) for s0 in range(0, n, a.chunk)]
        if a.meta: rows = np.array([pos_in_full[int(i)] for i in sel]); conf = (lab_ext[rows] if lab_ext is not None else np.ones(n, int)); closed = np.zeros(n, bool)
        else: conf = np.array([gt.get(i + 1, (0, 0))[0] for i in sel]); closed = np.array([gt.get(i + 1, (0, 0))[1] for i in sel]) == 1
        bag = None
        with ctx.Pool(a.nproc) as pool: parts = pool.map(_load_chunk, chunks)
        G['F'] = np.concatenate([p[0] for p in parts]); G['C'] = np.concatenate([p[1] for p in parts]); G['R'] = np.concatenate([p[2] for p in parts]); G['eo'] = eo; del parts
        if FLIB is not None: ll = fort_estep(FLIB, chunks)
        else:
            with ctx.Pool(a.nproc) as pool: parts = pool.map(_estep_chunk, chunks)
            ll = np.concatenate(parts); del parts
        L = ll - ll.max(1, keepdims=True); r = np.exp(L); r /= r.sum(1, keepdims=True)
        print('%s final E-step over all %d particles (%.0f s)' % (T, n, time.time() - t0), flush=True)
    # ---- outputs
    res = readout(T + ' final', r, i1, i2, ngrid, conf, closed, axis, rng, layer_of=layer_of, i3=i3)
    if os.environ.get('KLEM_DUMP', ''):   # Fortran port reference: the last E-step's inputs and outputs (band, noise, scale, per-particle cell log-likelihoods) and the half cell volumes
        try:
            Bd_ = G['Bd']; n3_ = G['REF'][0][0][0].shape[0]; ifft3 = lambda X: np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(X))))
            np.savez_compressed(a.out + '/klem_estep_ref.npz', ll=ll.astype(np.float32), pinds=sel + 1, eo=G['eo'], kx=Bd_.kx, ky=Bd_.ky, wq=Bd_.wq, shell=Bd_.shell, kmax=G['kmax'], n3=n3_, sig2=G['sig2'], SC=G['SC'], conv=int(G['conv']), ss=G['ss'], fs=G['fs'], box=G['box'],
                                tab=G['tab'][sel][:, :10], R=G['R'], lp=lp, ncell=ncell, ngrid=ngrid, K=len(G['REF'][0]), sign=float(getattr(G['D'], 'sign', 1.)), smpd=SMPD, F200=G['F'][:200], C200=G['C'][:200], R200=G['R'][:200])
            os.makedirs(a.out + '/cells_half', exist_ok=True)
            for hh in (0, 1):
                for k in range(len(G['REF'][hh])):
                    Xc = G['REF'][hh][k][0] + 1j * G['REF'][hh][k][1]; Xf = np.zeros((G['box'],) * 3, complex); cb = G['box'] // 2; hb = n3_ // 2; Xf[cb - hb:cb - hb + n3_, cb - hb:cb - hb + n3_, cb - hb:cb - hb + n3_] = Xc   # the coarse grid zero-padded to the full (even) box: same Fourier lattice, SIMPLE-readable
                    V3 = ifft3(Xf); write_mrc('%s/cells_half/cell_%02d_%s.mrc' % (a.out, k + 1, 'even' if hh == 0 else 'odd'), V3.astype(np.float32), SMPD)
            print('%s   KLEM_DUMP: klem_estep_ref.npz (ll %s) and %d half cell volumes (coarse grid %d^3 zero-padded to the full box)' % (T, ll.shape, 2 * len(G['REF'][0]), n3_), flush=True)
        except Exception as e: print('KLEM_DUMP failed:', e)
    with open(a.out + '/klem_latent.txt', 'w') as f:
        f.write('# pind z1 z2 argmax maxr gridmass layer2mass%s%s\n' % (' z3' if res['z3'] is not None else '', ''.join(' zb%d' % (j_ + 1) for j_ in range(G['zb'].shape[1])) if G.get('zb') is not None else ''))
        for j in range(n): f.write('%d %.4f %.4f %d %.4f %.4f %.4f%s%s\n' % (sel[j] + 1, res['z1'][j], res['z2'][j], res['am'][j] + 1, r[j].max(), r[j, :ngrid].sum(), r[j, :ngrid][layer_of >= 1].sum() / max(r[j, :ngrid].sum(), 1e-12), (' %.4f' % res['z3'][j]) if res['z3'] is not None else '', (''.join(' %.4f' % v_ for v_ in G['zb'][j]) if G.get('zb') is not None else '')))
    np.savez_compressed(a.out + '/klem_resp.npz', pinds=sel + 1, r=r.astype(np.float32), i1=i1, i2=i2); json.dump(hist, open(a.out + '/klem_hist.json', 'w'), indent=1)
    # per-state survival readout: cross-half reproducibility of every state's deviation from the consensus, population (argmax)
    try:
        n3_ = G['REF'][0][0][0].shape[0]; kk_ = np.arange(n3_) - n3_ // 2; rr_ = np.sqrt(sum(q ** 2 for q in np.meshgrid(kk_, kk_, kk_, indexing='ij'))); bm_ = (rr_ <= min(n3_ // 2, box * SMPD / 12.)) & (rr_ >= 1); pop = np.bincount(res['am'], minlength=K)   # survival judged at 12 A and coarser, whatever the band
        with open(a.out + '/klem_states.txt', 'w') as f:
            f.write('# state population cross_half_corr_of_deviation grid(1)/free(0)\n')
            crep = []
            for k in range(K):
                dv = [((G['REF'][hh][k][0] - Gc[hh].real) + 1j * (G['REF'][hh][k][1] - Gc[hh].imag))[bm_] for hh in (0, 1)]
                c_ = float(np.real(np.sum(np.conj(dv[0]) * dv[1])) / max(np.sqrt(np.sum(np.abs(dv[0]) ** 2) * np.sum(np.abs(dv[1]) ** 2)), 1e-30)); crep.append(c_); f.write('%d %d %.4f %d\n' % (k + 1, pop[k], c_, int(k < ngrid)))
        np.savez_compressed(a.out + '/klem_halves.npz', D0=np.array([((G['REF'][0][k][0] - Gc[0].real) + 1j * (G['REF'][0][k][1] - Gc[0].imag))[bm_] for k in range(ngrid)]).astype(np.complex64),
                            D1=np.array([((G['REF'][1][k][0] - Gc[1].real) + 1j * (G['REF'][1][k][1] - Gc[1].imag))[bm_] for k in range(ngrid)]).astype(np.complex64), w=r[:, :ngrid].sum(0), i1=i1[:ngrid], i2=i2[:ngrid], n3=n3_, band=np.where(bm_.ravel())[0])   # cell deviations per half at <= 12 A: the annealing split criterion reads these
        crep = np.array(crep[:ngrid]); print('%s states: cross-half reproducibility of the deviation median %.2f, %d of %d above 0.5 (survive), populations min/median/max %d/%d/%d' % (T, np.median(crep), int((crep > 0.5).sum()), ngrid, pop[:ngrid].min(), int(np.median(pop[:ngrid])), pop[:ngrid].max()), flush=True)
    except Exception as e: print('state readout skipped:', e)
    # cell maps on the coarse lattice (full map = mean of halves), as MRC for inspection
    smpd_c = box * SMPD / n3; d = a.out + '/cells'; os.makedirs(d, exist_ok=True)
    for k in range(K):
        Gk = 0.5 * sum((G['REF'][hh][k][0] + 1j * G['REF'][hh][k][1]) for hh in (0, 1)); v = np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(Gk))))
        name = (('state_%03d.mrc' % k) if a.kernel == 'graph' else ('cell_%02d_%02d.mrc' % (i1[k], i2[k]) if a.n3 == 1 else 'cell_%02d_%02d_%02d.mrc' % (i1[k], i2[k], i3[k])) if a.nlayer == 1 else 'cell_L%d_%02d_%02d.mrc' % (layer_of[k] + 1, i1[k], i2[k])) if k < ngrid else 'free_%d.mrc' % (k - ngrid + 1); write_mrc('%s/%s' % (d, name), v, smpd_c)
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        good = conf > 0; fig, ax = plt.subplots(1, 2, figsize=(11, 5))
        H, _, _ = np.histogram2d(res['z1'][good], res['z2'][good], bins=[np.linspace(-0.5, a.n1 - 0.5, 4 * a.n1 + 1), np.linspace(-0.5, a.n2 - 0.5, 4 * a.n2 + 1)]); ax[0].imshow(H.T, origin='lower', extent=(-0.5, a.n1 - 0.5, -0.5, a.n2 - 0.5), cmap='Greys', aspect='auto'); ax[0].set_title('posterior-mean positions, conformer particles')
        cls = np.unique(conf[good]); cmap = plt.get_cmap('tab10')
        for j, cc_ in enumerate(cls):
            m = good & (conf == cc_); ax[1].scatter(res['z1'][m][:2000], res['z2'][m][:2000], s=2, alpha=0.3, color=cmap(j % 10), label='%d (axis %+.2f)' % (cc_, axis.get(cc_, 0)))
            ax[1].text(res['z1'][m].mean(), res['z2'][m].mean(), str(cc_), fontsize=9, weight='bold')
        ax[1].set_xlim(-0.5, a.n1 - 0.5); ax[1].set_ylim(-0.5, a.n2 - 0.5); ax[1].legend(fontsize=6, markerscale=4); ax[1].set_title('per conformer (truth colours)'); plt.tight_layout(); plt.savefig(a.out + '/fig_klem_latent.png', dpi=110)
    except Exception as e: print('figure skipped:', e)
    print('%s DONE (%.0f s) -> %s' % (T, time.time() - t0, a.out))


if __name__ == '__main__':
    main()
