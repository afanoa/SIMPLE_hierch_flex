#!/usr/bin/env python3
"""Deformation-latent prototype on the window maps of the latent-windows route (toy with truth):
  1. every window map registered to the canonical window map (middle window) by multi-scale demons -> flow fields u_w
  2. PCA of the flows across windows = deformation basis; mode-1 score per window vs window index / truth axis
  3. flow-corrected consensus: every window map warped INTO the canonical frame and population-averaged, even and odd
     halves separately -> FSC gain over the plain average; correlation with the canonical window's truth conformer
  4. per-particle deformation E-step on a subset: references = canonical (even/odd) warped along mode 1 at a grid of
     latent values, scored cross-half with the latent_windows E-step -> Spearman with the truth axis and conformer
     purity, against the window E-step on the same subset and band.
usage: deform_proto.py MAPDIR OUT --nw 8 --info windows.txt --cons toy_consensus.mrc --truth f1,.. --states .. --axis ..
       --stk all.mrcs --table ptcl3D.txt --gt gt_particles.txt [--nsub 20000] [--lp 12] [--ngrid 9] [--sigma 8] [--nproc 8]"""
import sys, os, argparse, time, types
import numpy as np
os.environ.setdefault('TOY_SMPD', '1.464')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model_class import read_mrc, write_mrc, Shells, fftc, BANDS, demons, warp, flow_pca
from latent_windows import cmd_estep, parse_axis
from scipy.ndimage import map_coordinates
from scipy.stats import spearmanr
from multiprocessing import Pool

SMPD = float(os.environ['TOY_SMPD'])


def inverse_field(u, iters=6):
    v = -u.copy(); B = u.shape[1]; zz, yy, xx = np.mgrid[:B, :B, :B]
    for _ in range(iters):
        v = -np.array([map_coordinates(u[k], [zz + v[0], yy + v[1], xx + v[2]], order=1, mode='nearest') for k in range(3)])
    return v


def reg_one(args):  # sigma = elastic regularisation (voxels)
    w, fn_w, fn_can, smpd, sigma = args
    S = Shells(160, smpd); Mw = read_mrc(fn_w); Mc = read_mrc(fn_can); S = Shells(Mw.shape[0], smpd); msk = S.soft_mask(100. / smpd)
    u = demons(Mw * msk, Mc * msk, S, sigma_elastic=sigma); return w, u.astype(np.float32)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('mapdir'); ap.add_argument('out'); ap.add_argument('--nw', type=int, default=8); ap.add_argument('--info', required=True)
    ap.add_argument('--cons', required=True); ap.add_argument('--truth', required=True); ap.add_argument('--states', required=True); ap.add_argument('--axis', required=True)
    ap.add_argument('--stk', default=''); ap.add_argument('--table', default=''); ap.add_argument('--gt', default=''); ap.add_argument('--nsub', type=int, default=20000)
    ap.add_argument('--lp', type=float, default=12.); ap.add_argument('--ngrid', type=int, default=9); ap.add_argument('--sigma', type=float, default=3.); ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--mskdiam', type=float, default=200.); ap.add_argument('--tag', default='deform')
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); t0 = time.time(); T = a.tag
    axis = parse_axis(a); states = [int(s) for s in a.states.split(',')]; truth = dict(zip(states, a.truth.split(',')))
    info = np.loadtxt(a.info, comments='#', ndmin=2); n = info[:, 1]; dom = info[:, 2].astype(int); axm = info[:, 5]
    can = int(np.argmax(n * (np.arange(a.nw) == a.nw // 2))) + 1 if False else a.nw // 2 + 1
    fn = lambda w, suf='': '%s/recvol_state%02d%s.mrc' % (a.mapdir, w, suf)
    cons = read_mrc(a.cons); B = cons.shape[0]; S = Shells(B, SMPD); msk = S.soft_mask(a.mskdiam / 2 / SMPD)
    print('%s: canonical = window %d (n %d, dominant conformer %d); registering %d windows (sigma %.0f px)' % (T, can, n[can - 1], dom[can - 1], a.nw, a.sigma), flush=True)
    with Pool(min(a.nproc, a.nw)) as pool: res = pool.map(reg_one, [(w, fn(w), fn(can), SMPD, a.sigma) for w in range(1, a.nw + 1)])
    U = np.zeros((a.nw, 3, B, B, B), np.float32)
    for w, u in res: U[w - 1] = u
    dens = S.lowpass(cons, 15); inside = dens > dens.mean() + 1.5 * dens.std()
    amp = [np.percentile(np.sqrt(np.sum(U[w] ** 2, 0))[inside] * SMPD, 90) for w in range(a.nw)]
    print('%s displacement p90 per window (A): %s (%.0f s)' % (T, ' '.join('%.1f' % v for v in amp), time.time() - t0), flush=True)
    # ---- flow PCA across windows (decimated, inside the density)
    Ud = U[:, :, ::2, ::2, ::2]; lam, scores, pr = flow_pca(Ud, inside[::2, ::2, ::2])
    s1 = scores[:, 0]; s1 = s1 * np.sign(spearmanr(s1, np.arange(a.nw))[0] or 1)
    print('%s flow PCA: effective dimension %.1f, variance shares %s; mode-1 score per window %s; Spearman(mode 1, window index) %+.2f, (mode 1, truth axis of windows) %+.2f' %
          (T, pr, np.array2string(lam[:4] / lam.sum(), precision=2), ' '.join('%+.2f' % v for v in s1 / np.abs(s1).max()), spearmanr(s1, np.arange(a.nw))[0], spearmanr(s1, axm)[0]), flush=True)
    # ---- flow-corrected consensus (windows warped INTO the canonical frame), even / odd
    def fsc_res(x, y):
        a_, b_ = fftc(x * msk), fftc(y * msk); sh = S.sh; num = np.bincount(sh.ravel(), np.real(np.conj(a_) * b_).ravel()); d1 = np.bincount(sh.ravel(), np.abs(a_.ravel()) ** 2); d2 = np.bincount(sh.ravel(), np.abs(b_.ravel()) ** 2)
        f = num / np.maximum(np.sqrt(d1 * d2), 1e-30); s = 1
        while s < B // 2 and f[s] > 0.143: s += 1
        return B * SMPD / s, f
    def bc(x, y, lo, hi):
        X, Y = fftc(x * msk), fftc(y * msk); m = S.band(lo, hi); return float(np.real(np.sum(np.conj(X[m]) * Y[m])) / max(np.sqrt(np.sum(np.abs(X[m]) ** 2) * np.sum(np.abs(Y[m]) ** 2)), 1e-30))
    V = {}
    for suf in ('_even', '_odd'):
        plain = np.zeros((B, B, B)); corr = np.zeros((B, B, B))
        for w in range(1, a.nw + 1):
            M = read_mrc(fn(w, suf)); plain += n[w - 1] * M; v = inverse_field(U[w - 1].astype(float)); corr += n[w - 1] * warp(M, v)
        V['plain' + suf] = plain / n.sum(); V['flow' + suf] = corr / n.sum()
    rp, fp = fsc_res(V['plain_even'], V['plain_odd']); rf, ff = fsc_res(V['flow_even'], V['flow_odd'])
    Tc = read_mrc(truth[dom[can - 1]]) - cons; Pf = 0.5 * (V['plain_even'] + V['plain_odd']) - cons; Ff = 0.5 * (V['flow_even'] + V['flow_odd']) - cons
    print('%s flow-corrected consensus: FSC 0.143 plain average %.1f A -> flow-corrected %.1f A; corr with (truth conformer %d - consensus) per band: plain %s | flow-corrected %s' %
          (T, rp, rf, dom[can - 1], ' '.join('%.2f' % bc(Pf, Tc, lo, hi) for lo, hi in BANDS), ' '.join('%.2f' % bc(Ff, Tc, lo, hi) for lo, hi in BANDS)), flush=True)
    write_mrc(a.out + '/flow_corrected.mrc', 0.5 * (V['flow_even'] + V['flow_odd']), SMPD); write_mrc(a.out + '/plain_average.mrc', 0.5 * (V['plain_even'] + V['plain_odd']), SMPD)
    np.save(a.out + '/flows_dec2.npy', Ud); np.save(a.out + '/flow_scores.npy', scores)
    # ---- per-particle deformation E-step: canonical warped along mode 1 at a grid of latent values (even/odd canonical, same fields)
    if a.stk:
        ubar = U.mean(0).astype(float); e1 = scores[:, 0]; lo, hi = e1.min(), e1.max(); grid = np.linspace(lo - 0.15 * (hi - lo), hi + 0.15 * (hi - lo), a.ngrid)
        Xd = Ud.reshape(a.nw, -1).astype(float); Xc = Xd - Xd.mean(0); mode = (Xc.T @ scores[:, 0]) / max(np.sum(scores[:, 0] ** 2), 1e-30)   # unit-score field (decimated)
        from scipy.ndimage import zoom
        mode = mode.reshape(3, B // 2, B // 2, B // 2); mode_full = np.array([zoom(c, 2, order=1) for c in mode]) * 2
        gdir = a.out + '/grid_refs'; os.makedirs(gdir, exist_ok=True)
        for suf in ('_even', '_odd'):
            Mc = read_mrc(fn(can, suf))
            for g, z in enumerate(grid): write_mrc('%s/recvol_state%02d%s.mrc' % (gdir, g + 1, suf), warp(Mc, ubar + z * mode_full), SMPD)
        for g, z in enumerate(grid): write_mrc('%s/recvol_state%02d.mrc' % (gdir, g + 1), warp(0.5 * (read_mrc(fn(can, '_even')) + read_mrc(fn(can, '_odd'))), ubar + z * mode_full), SMPD)
        print('%s grid references: %d latent values %s (window scores %s)' % (T, a.ngrid, np.array2string(grid, precision=2), np.array2string(e1, precision=2)), flush=True)
        base = dict(stk=a.stk, table=a.table, gt=a.gt, refs='', lp=a.lp, states=a.states, axis=a.axis, prev='', prior_sigma=1.5, nsub=a.nsub, mskdiam=a.mskdiam)
        cmd_estep(types.SimpleNamespace(mapdir=gdir, nw=a.ngrid, out=a.out + '/coord_deform.txt', tag=T + ' grid', **base))
        cmd_estep(types.SimpleNamespace(mapdir=a.mapdir, nw=a.nw, out=a.out + '/coord_windows.txt', tag=T + ' windows', **base))
    print('%s DONE (%.0f s)' % (T, time.time() - t0))


if __name__ == '__main__':
    main()
