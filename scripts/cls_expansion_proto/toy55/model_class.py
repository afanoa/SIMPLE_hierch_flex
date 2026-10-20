#!/usr/bin/env python3
"""Model-class test on the TRUTH maps of the RNase P toy (no particles): which latent decoder can represent the
conformer landscape, per band? For every conformer map V_i against the consensus V_c (d_i = V_i - V_c):
  linear K   : d_i in the span of the OTHER conformers' differences (leave-one-out PCA basis, K components)
  deform     : V_i ~ V_c(x + u(x)), u a smooth displacement field (multi-scale symmetric demons, Gaussian-regularised)
  deform+occ : V_i ~ o(x) V_c(x + u(x)), o a smooth occupancy map in [0, 1]
Readout per conformer and band (117-30, 30-20, 20-12, 12-8 A): unexplained fraction |resid|^2 / |d_i|^2 (1 = nothing
explained, 0 = exact). Also the displacement amplitude per conformer (A; 50th/90th percentile inside the density),
the PCA spectrum of the flow fields (participation ratio = effective dimension of the deformation space) against the
PCA spectrum of the volume differences per band (effective dimension a linear volume basis needs).
usage: model_class.py OUT CONSENSUS map1,map2,... [--names n1,..] [--smpd 1.464] [--mskdiam 200] [--sigma 8] [--nproc 8]"""
import sys, os, json, argparse, time
import numpy as np
from scipy.ndimage import map_coordinates, gaussian_filter, zoom

BANDS = [(117, 30), (30, 20), (20, 12), (12, 8)]


def read_mrc(fn):
    h = np.fromfile(fn, np.int32, 256); nx, ny, nz = int(h[0]), int(h[1]), int(h[2])
    return np.fromfile(fn, np.float32, offset=1024 + int(h[23]), count=nx * ny * nz).reshape(nz, ny, nx).astype(np.float64)


def write_mrc(fn, vol, smpd):
    v = np.ascontiguousarray(vol, np.float32); nz, ny, nx = v.shape
    h = np.zeros(256, np.int32); hf = h.view(np.float32)
    h[0:3] = (nx, ny, nz); h[3] = 2; h[7:10] = (nx, ny, nz); hf[10:13] = (nx * smpd, ny * smpd, nz * smpd); hf[13:16] = 90.; h[16:19] = (1, 2, 3)
    hf[19:22] = (v.min(), v.max(), v.mean()); h[52] = int.from_bytes(b'MAP ', 'little'); h[53] = int.from_bytes(bytes([68, 65, 0, 0]), 'little'); hf[54] = v.std()
    with open(fn, 'wb') as fo: fo.write(h.tobytes()); fo.write(v.tobytes())


def fftc(x): return np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(x)))
def ifftc(F): return np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(F))))


class Shells:
    def __init__(self, B, smpd):
        zz, yy, xx = np.mgrid[:B, :B, :B]; self.r = np.hypot(np.hypot(xx - B / 2, yy - B / 2), zz - B / 2); self.B = B; self.smpd = smpd
        self.sh = np.rint(self.r).astype(int)

    def kof(self, res): return self.B * self.smpd / res

    def band(self, lo, hi): return (self.sh >= self.kof(lo)) & (self.sh <= self.kof(hi))

    def lowpass(self, v, lp, width=2.):
        k = self.kof(lp); w = 0.5 * (1 + np.cos(np.pi * np.clip((self.r - (k - width)) / width, 0, 1)))
        return ifftc(fftc(v) * w)

    def soft_mask(self, radius_px, width=4.):
        return 0.5 * (1 + np.cos(np.pi * np.clip((self.r - (radius_px - width)) / width, 0, 1)))


def band_unexplained(d, resid, S):
    Fd, Fr = fftc(d), fftc(resid)
    return [float(np.sum(np.abs(Fr[m]) ** 2) / max(np.sum(np.abs(Fd[m]) ** 2), 1e-30)) for m in [S.band(lo, hi) for lo, hi in BANDS]]


def band_spectrum(D, S):
    """PCA spectrum of the volume differences per band -> participation ratio (effective dimension) and the
    variance share of the first component"""
    out = {}
    F = [fftc(d) for d in D]
    for lo, hi in BANDS:
        m = S.band(lo, hi); X = np.array([f[m] for f in F]); G = np.real(X @ X.conj().T); lam = np.maximum(np.linalg.eigvalsh(G)[::-1], 0)
        out['%d-%d' % (lo, hi)] = dict(pr=float(lam.sum() ** 2 / max(np.sum(lam ** 2), 1e-30)), share1=float(lam[0] / max(lam.sum(), 1e-30)),
                                      share2=float(lam[:2].sum() / max(lam.sum(), 1e-30)))
    return out


def warp(V, u):
    """V(x + u(x)); u = (3, B, B, B) in voxels, components (dz, dy, dx)"""
    B = V.shape[0]; zz, yy, xx = np.mgrid[:B, :B, :B]
    return map_coordinates(V, [zz + u[0], yy + u[1], xx + u[2]], order=1, mode='nearest')


def demons_level(fixed, moving, u, iters, sigma_fluid, sigma_elastic, alpha=2.):
    gF = np.array(np.gradient(fixed))
    for it in range(iters):
        Mw = warp(moving, u); gM = np.array(np.gradient(Mw)); g = 0.5 * (gF + gM); diff = Mw - fixed
        den = np.sum(g ** 2, 0) + diff ** 2 / alpha ** 2
        step = -diff * g / np.maximum(den, 1e-12)
        step = np.array([gaussian_filter(s, sigma_fluid) for s in step])
        u = u + step
        if sigma_elastic > 0: u = np.array([gaussian_filter(c, sigma_elastic) for c in u])
    return u


def demons(fixed, moving, S, sigma_elastic=3., scales=((30, 50, 2), (20, 40, 2), (14, 30, 1), (10, 20, 1)), sigma_fluid=5.):
    """u with moving(x + u) ~ fixed. scales: (low-pass A, iterations, decimation): the coarse levels run on a 2x
    decimated lattice (the field is up-sampled between levels). Returns u (3, B, B, B) in voxels of the full lattice."""
    B = fixed.shape[0]; u = None; cur_dec = None
    for lp, iters, dec in scales:
        Fx = S.lowpass(fixed, lp); Mv = S.lowpass(moving, lp)
        if dec > 1:
            Fx = Fx[::dec, ::dec, ::dec]; Mv = Mv[::dec, ::dec, ::dec]
        if u is None: u = np.zeros((3,) + Fx.shape)
        elif cur_dec != dec:
            f = cur_dec / dec; u = np.array([zoom(c, f, order=1) for c in u]) * f
        cur_dec = dec
        u = demons_level(Fx, Mv, u, iters, sigma_fluid / dec, sigma_elastic / dec)
    if cur_dec != 1: f = cur_dec; u = np.array([zoom(c, f, order=1) for c in u]) * f
    return u


def occupancy(W, T, sigma=6.):
    num = gaussian_filter(W * T, sigma); den = gaussian_filter(W * W, sigma)
    return np.clip(num / np.maximum(den, 1e-6 * den.max()), 0., 1.)


def flow_pca(U, inside):
    """PCA of the flow fields (n, 3, B, B, B) over the voxels 'inside' (bool): spectrum, participation ratio, scores"""
    X = np.array([u[:, inside].ravel() for u in U]); Xc = X - X.mean(0); G = Xc @ Xc.T; lam, E = np.linalg.eigh(G); lam, E = lam[::-1], E[:, ::-1]
    lam = np.maximum(lam, 0); scores = E * np.sqrt(lam)
    return lam, scores, float(lam.sum() ** 2 / max(np.sum(lam ** 2), 1e-30))


def one_conformer(args):
    i, name, fn, cons_fn, smpd, mskdiam, sigma, out = args
    t0 = time.time()
    Vc = read_mrc(cons_fn); Vi = read_mrc(fn); B = Vc.shape[0]; S = Shells(B, smpd); msk = S.soft_mask(mskdiam / 2 / smpd)
    Vc, Vi = Vc * msk, Vi * msk; d = Vi - Vc
    dens = S.lowpass(Vc, 15) > S.lowpass(Vc, 15).mean() + 1.5 * S.lowpass(Vc, 15).std()      # density support for amplitude stats
    res = dict(name=name, bands=['%d-%d' % b for b in BANDS])
    u = demons(Vi, Vc, S, sigma_elastic=sigma)
    Wv = warp(Vc, u) * msk; res['deform'] = band_unexplained(d, (Vi - Wv), S)
    o = occupancy(Wv, Vi); res['deform_occ'] = band_unexplained(d, Vi - o * Wv, S)
    amp = np.sqrt(np.sum(u ** 2, 0))[dens] * smpd; res['disp_p50_A'] = float(np.percentile(amp, 50)); res['disp_p90_A'] = float(np.percentile(amp, 90)); res['disp_max_A'] = float(amp.max())
    res['occ_min'] = float(np.percentile(o[dens], 5))
    # rigid-only reference: the best global shift (the trivial deformation), to show the field is not just a shift
    sh = np.array([np.mean(u[k][dens]) for k in range(3)]); us = np.broadcast_to(sh[:, None, None, None], u.shape).copy()
    res['shift_only'] = band_unexplained(d, Vi - warp(Vc, us) * msk, S); res['shift_px'] = sh.tolist()
    np.save(out + '/flow_%s.npy' % name, u[:, ::2, ::2, ::2].astype(np.float32))
    write_mrc(out + '/warped_%s.mrc' % name, Wv, smpd)
    res['time_s'] = time.time() - t0
    return res


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('cons'); ap.add_argument('maps')
    ap.add_argument('--names', default=''); ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--mskdiam', type=float, default=200.)
    ap.add_argument('--sigma', type=float, default=3., help='elastic regularisation of the field (voxels); updates are smoothed with 5 voxels'); ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--ks', default='1,2,3,5', help='linear leave-one-out ranks')
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    fns = a.maps.split(','); names = a.names.split(',') if a.names else [os.path.basename(f).split('_')[0] for f in fns]; n = len(fns)
    Vc = read_mrc(a.cons); B = Vc.shape[0]; S = Shells(B, a.smpd); msk = S.soft_mask(a.mskdiam / 2 / a.smpd); Vc *= msk
    V = [read_mrc(f) * msk for f in fns]; D = [v - Vc for v in V]
    print('model_class: %d maps, box %d, smpd %.3f; |d_i|/|V_c| = %s' % (n, B, a.smpd, ' '.join('%s:%.2f' % (nm, np.linalg.norm(d) / np.linalg.norm(Vc)) for nm, d in zip(names, D))), flush=True)
    # ---- linear leave-one-out
    Ks = [int(k) for k in a.ks.split(',')]; lin = {}
    X = np.array([d.ravel() for d in D])
    for i in range(n):
        Xo = np.delete(X, i, 0); G = Xo @ Xo.T; lam, E = np.linalg.eigh(G); order = np.argsort(lam)[::-1]; lam, E = lam[order], E[:, order]
        keep = lam > 1e-10 * lam[0]; Vt = (E[:, keep].T @ Xo) / np.sqrt(lam[keep])[:, None]        # orthonormal rows spanning the others
        lin[names[i]] = {}
        for K in Ks + ['all']:
            k = len(lam[keep]) if K == 'all' else min(K, len(lam[keep])); coef = Vt[:k] @ X[i]; resid = X[i] - Vt[:k].T @ coef
            lin[names[i]][str(K)] = band_unexplained(D[i], resid.reshape(Vc.shape), S)
        print('   linear LOO %s: ' % names[i] + '; '.join('K=%s %s' % (K, ' '.join('%.2f' % v for v in lin[names[i]][str(K)])) for K in Ks + ['all']), flush=True)
    spec = band_spectrum(D, S)
    print('volume-difference PCA per band (effective dimension, share of PC1, PC1+2): ' + '; '.join('%s A: %.1f %.2f %.2f' % (b, v['pr'], v['share1'], v['share2']) for b, v in spec.items()), flush=True)
    # ---- deformation (parallel over conformers)
    from multiprocessing import Pool
    jobs = [(i, names[i], fns[i], a.cons, a.smpd, a.mskdiam, a.sigma, a.out) for i in range(n)]
    with Pool(min(a.nproc, n)) as pool: res = pool.map(one_conformer, jobs)
    for r in res:
        print('   deform %s: shift-only %s | deform %s | deform+occ %s | displacement p50/p90/max %.1f/%.1f/%.1f A, occupancy p5 %.2f (%.0f s)' %
              (r['name'], ' '.join('%.2f' % v for v in r['shift_only']), ' '.join('%.2f' % v for v in r['deform']), ' '.join('%.2f' % v for v in r['deform_occ']),
               r['disp_p50_A'], r['disp_p90_A'], r['disp_max_A'], r['occ_min'], r['time_s']), flush=True)
    # ---- flow PCA (effective dimension of the deformation space) on the decimated fields inside the density
    U = np.array([np.load(a.out + '/flow_%s.npy' % nm) for nm in names]); dens = S.lowpass(Vc, 15)[::2, ::2, ::2]; inside = dens > dens.mean() + 1.5 * dens.std()
    lam, scores, pr = flow_pca(U, inside)
    print('flow-field PCA: effective dimension %.1f, variance shares %s; mode-1 scores %s' % (pr, np.array2string(lam[:5] / lam.sum(), precision=2), ' '.join('%s:%+.2f' % (nm, s) for nm, s in zip(names, scores[:, 0] / np.abs(scores[:, 0]).max()))), flush=True)
    # ---- summary table: mean over conformers per band
    def mean_rows(key, sub=None):
        rows = [lin[nm][sub] for nm in names] if key == 'lin' else [r[key] for r in res]; return np.mean(rows, 0)
    print('SUMMARY unexplained fraction per band (mean over conformers; bands %s A):' % ', '.join('%d-%d' % b for b in BANDS))
    for K in Ks: print('   linear K=%-3s %s' % (K, ' '.join('%.2f' % v for v in mean_rows('lin', str(K)))))
    print('   (K=all is exact by construction: the consensus is the weighted mean of the maps)')
    for key in ('shift_only', 'deform', 'deform_occ'): print('   %-11s  %s' % (key, ' '.join('%.2f' % v for v in mean_rows(key))))
    json.dump(dict(names=names, linear=lin, deform=res, spectrum=spec, flow_pr=pr, flow_lambda=lam.tolist(), flow_scores=scores[:, :3].tolist(), sigma=a.sigma), open(a.out + '/results.json', 'w'), indent=1)
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.2)); x = np.arange(len(BANDS)); w = 0.13
        series = [('linear K=1', mean_rows('lin', '1')), ('linear K=2', mean_rows('lin', '2')), ('linear K=5', mean_rows('lin', '5')), ('linear all others', mean_rows('lin', 'all')),
                  ('shift only', mean_rows('shift_only')), ('deform', mean_rows('deform')), ('deform+occ', mean_rows('deform_occ'))]
        for j, (lab, vals) in enumerate(series): ax[0].bar(x + (j - 3) * w, vals, w, label=lab)
        ax[0].set_xticks(x); ax[0].set_xticklabels(['%d-%d A' % b for b in BANDS]); ax[0].set_ylabel('unexplained fraction of V_i - V_c'); ax[0].set_ylim(0, 1.05); ax[0].legend(fontsize=7, ncol=2)
        ax[0].set_title('what each model class leaves unexplained (mean over %d conformers)' % n, fontsize=10)
        ax[1].bar(np.arange(n), [r['disp_p90_A'] for r in res]); ax[1].set_xticks(np.arange(n)); ax[1].set_xticklabels(names, rotation=60, fontsize=8); ax[1].set_ylabel('displacement p90 (A)')
        ax[1].set_title('deformation amplitude vs consensus; flow space effective dim %.1f' % pr, fontsize=10)
        plt.tight_layout(); plt.savefig(a.out + '/fig_model_class.png', dpi=110)
    except Exception as e: print('figure skipped:', e)
    print('DONE model_class ->', a.out)


if __name__ == '__main__':
    main()
