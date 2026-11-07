#!/usr/bin/env python3
"""Prototypes of audit 6.1 (null-calibrated mask core) and 6.2 (cross-validated SVD site factors) on a node's saved half maps.

mask: recompute falloff_mask's cross-half covariance v (15 A low-pass, least-squares scaling to the weighted mean inside the
      support, Gaussian 1.5 px, sphere), check it against mask_variance.mrc, then compare the 0.2-of-peak core with a
      null-calibrated core z = v / sd_null > Z where sd_null comes from noise surrogates built from the half differences
      through the same pipeline; cores are 6-connected components of at least DMINVOX voxels; both cores get the 3 px
      dilation and 6 px Gaussian falloff; molecule fractions and Dice are reported.
svd:  for every component of either core, the voxel-by-state deviation matrices D_e, D_o (columns weighted by sqrt(w)), the
      SVD of D_e and the odd half projected on its right singular vectors: per factor the cross-half correlation of its
      spatial map against the noise ceiling (1 / sqrt(n_eff), n_eff from surrogate noise pairs on the same voxels), beside
      the 2-means verdict of the tree's discovery log when present.
usage: mask_disc_proto.py NODE_DIR [--maps flex|cheap|recB0] [--z 3] [--nsur 24] [--mskdiam A]
"""
import sys, os, re, glob, argparse, math
import numpy as np
from scipy import ndimage

MLP, MSIG, MFRAC, MDIL, MFALL, MSUP, DMINVOX = 15., 1.5, 0.2, 3., 6., 0.1, 200

def read_mrc(fn):
    with open(fn, 'rb') as f:
        h = f.read(1024)
    nx, ny, nz, mode = np.frombuffer(h[:16], dtype='<i4')
    cella = np.frombuffer(h[40:52], dtype='<f4')
    smpd = float(cella[0]) / nx if cella[0] > 0 else 1.0
    dt = {0: '<i1', 1: '<i2', 2: '<f4', 6: '<u2'}[int(mode)]
    nsym = int(np.frombuffer(h[92:96], dtype='<i4')[0])
    d = np.fromfile(fn, dtype=dt, offset=1024 + nsym).astype(np.float32)[:nx * ny * nz]
    return d.reshape((nz, ny, nx)).transpose(2, 1, 0).copy(), smpd   # (x, y, z) like the Fortran rmat

def fourier_filter(v, smpd, mode):
    n = v.shape[0]; f = np.fft.fftfreq(n)
    fx, fy, fz = np.meshgrid(f, f, f, indexing='ij'); f2 = fx**2 + fy**2 + fz**2
    fac = np.exp(-f2 / smpd**2 * MLP**2 * math.log(2.)) if mode == 1 else np.exp(-2. * math.pi**2 * MSIG**2 * f2)
    return np.real(np.fft.ifftn(np.fft.fftn(v) * fac)).astype(np.float32)

def sphere(n, rad):
    c = n // 2; ax = np.arange(n) - c
    x, y, z = np.meshgrid(ax, ax, ax, indexing='ij')
    return (x**2 + y**2 + z**2) <= rad**2

def dist_px(mask):
    return ndimage.distance_transform_edt(~mask)

def falloff(core, sph):
    if not core.any(): return sph.astype(np.float32)
    core2 = dist_px(core) <= MDIL
    m = np.exp(-dist_px(core2)**2 / (2. * MFALL**2)).astype(np.float32); m[~sph] = 0.
    return m

def components(core):
    lab, nl = ndimage.label(core, structure=ndimage.generate_binary_structure(3, 1))
    sizes = ndimage.sum(core, lab, range(1, nl + 1)).astype(int) if nl else np.array([], int)
    keep = np.zeros_like(core)
    for k, s in enumerate(sizes, 1):
        if s >= DMINVOX: keep |= lab == k
    return keep, sorted([int(s) for s in sizes if s >= DMINVOX], reverse=True), lab

def load_maps(d, which):
    sub = {'flex': 'flex', 'cheap': 'cheap', 'recB0': 'recB0'}[which]
    evens = sorted(glob.glob(os.path.join(d, sub, 'recvol_state??_even_unfil.mrc')))
    odds = [e.replace('_even_', '_odd_') for e in evens]
    states = [int(re.search(r'state(\d\d)', e).group(1)) for e in evens]
    w = np.ones(len(evens))
    coord = os.path.join(d, sub, 'flex_pca_coordinates.txt')
    if which == 'flex' and os.path.exists(coord):
        lab = np.loadtxt(coord, comments='#', usecols=2).astype(int)
        w = np.array([max((lab == s).sum(), 1) for s in states], float)
    return evens, odds, states, w

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('node'); ap.add_argument('--maps', default='flex')
    ap.add_argument('--z', type=float, default=3.); ap.add_argument('--nsur', type=int, default=24); ap.add_argument('--mskdiam', type=float, default=0.)
    a = ap.parse_args(); rng = np.random.default_rng(0)
    evens, odds, states, w = load_maps(a.node, a.maps)
    if len(evens) < 2: sys.exit('no half maps under %s/%s' % (a.node, a.maps))
    E = []; O = []
    for e, o in zip(evens, odds):
        ve, smpd = read_mrc(e); vo, _ = read_mrc(o)
        E.append(fourier_filter(ve, smpd, 1)); O.append(fourier_filter(vo, smpd, 1))
    E = np.array(E); O = np.array(O); n = E.shape[1]; wn = w / w.sum()
    mskdiam = a.mskdiam
    if mskdiam <= 0:
        m = re.search(r'mskdiam=([0-9.]+)', open(os.path.join(os.path.dirname(a.node.rstrip('/')), '..', 'run_manifest.txt')).read()) if os.path.exists(os.path.join(os.path.dirname(a.node.rstrip('/')), '..', 'run_manifest.txt')) else None
        mskdiam = float(m.group(1)) if m else 0.8 * n * smpd
    sph = sphere(n, mskdiam / 2. / smpd)
    em = np.tensordot(wn, E, 1); om = np.tensordot(wn, O, 1)
    mean = 0.5 * (em + om); sup = (mean > MSUP * mean[sph].max()) & sph
    for j in range(len(E)):
        E[j] *= (E[j][sup] * em[sup]).sum() / max((E[j][sup]**2).sum(), 1e-30)
        O[j] *= (O[j][sup] * om[sup]).sum() / max((O[j][sup]**2).sum(), 1e-30)
    em = np.tensordot(wn, E, 1); om = np.tensordot(wn, O, 1)
    DE = E - em; DO = O - om
    v = np.tensordot(wn, DE * DO, 1); vs = fourier_filter(v, smpd, 2); vs[~sph] = 0.
    print('node %s: %d %s maps, box %d, %.3f A/px, mskdiam %.0f A, weights %s' % (a.node, len(E), a.maps, n, smpd, mskdiam, np.round(wn, 3).tolist()))
    mv = os.path.join(a.node, 'mask_variance.mrc')
    if os.path.exists(mv):
        ref, _ = read_mrc(mv); c = np.corrcoef(ref[sph], vs[sph])[0, 1]
        print('  recomputed v vs mask_variance.mrc: corr %.4f, peak ratio %.3f' % (c, vs.max() / max(ref.max(), 1e-30)))
    # ---- 6.1 null surrogates from the half differences
    D = E - O; Dm = np.tensordot(wn, D, 1); Dc = D - Dm; K = len(E)
    sur = []
    for r in range(a.nsur):
        perm = rng.permutation(K)
        while (perm == np.arange(K)).any(): perm = rng.permutation(K)
        vr = np.zeros_like(v)
        for j in range(K): vr += wn[j] * Dc[j] * Dc[perm[j]] * 0.5 * math.sqrt(w[perm[j]] / w[j])
        vr = fourier_filter(vr, smpd, 2); sur.append(vr[sph])
    sur = np.array(sur)
    sd_pooled = 1.4826 * np.median(np.abs(sur - np.median(sur)))
    z = np.zeros_like(vs); z[sph] = vs[sph] / max(sd_pooled, 1e-30)
    core_old = vs > MFRAC * vs.max()
    core_new_raw = z > a.z
    core_new, sizes_new, lab_new = components(core_new_raw)
    core_old_c, sizes_old, lab_old = components(core_old)
    m_old = falloff(core_old, sph); m_new = falloff(core_new, sph)
    frac = lambda m: (m > 0.5)[sup].sum() / max(sup.sum(), 1)
    dice = 2 * (core_old & core_new).sum() / max(core_old.sum() + core_new.sum(), 1)
    print('  null sd of the smoothed covariance (pooled MAD over %d surrogates) %.4e; peak v %.4e = z %.1f; z > %.0f on %.2f %% of the sphere (%.3f %% expected by chance)'
          % (a.nsur, sd_pooled, vs.max(), vs.max() / sd_pooled, a.z, 100. * core_new_raw[sph].mean(), 100. * 0.5 * math.erfc(a.z / math.sqrt(2))))
    print('  OLD core (0.2 x peak): %d voxels, components >= %d: %s; mask covers %.0f %% of the molecule support' % (core_old.sum(), DMINVOX, sizes_old[:8], 100 * frac(m_old)))
    print('  NEW core (z > %.0f, components >= %d): %d voxels, components: %s; mask covers %.0f %% of the molecule support; Dice(old core, new core) %.3f'
          % (a.z, DMINVOX, core_new.sum(), sizes_new[:8], 100 * frac(m_new), dice))
    # ---- 6.2 cross-validated SVD per component
    sqw = np.sqrt(wn)
    # the noise ceiling: correlation between two independent surrogate noise maps restricted to the component's voxels
    def neff_of(idx):
        rs = []
        for r in range(0, len(sur) - 1, 2):
            rs.append(np.corrcoef(sur[r][idx_in_sph[idx]], sur[r + 1][idx_in_sph[idx]])[0, 1])
        return 1. / max(np.var(rs), 1e-12) if len(rs) > 2 else float('nan')
    idx_in_sph = np.full(vs.shape, -1, int); idx_in_sph[sph] = np.arange(sph.sum())
    disc_log = os.path.join(a.node, 'sites', 'hetero_step.log'); verdicts = []
    if os.path.exists(disc_log):
        for l in open(disc_log, errors='replace'):
            m = re.search(r'component (\d+) (splits into two factors|is one factor)(.*)', l)
            if m: verdicts.append((int(m.group(1)), m.group(2), m.group(3).strip()))
    for tag, lab, sizes in (('old', lab_old, sizes_old), ('new', lab_new, sizes_new)):
        print('  SVD factors on the %s core components:' % tag)
        comp_ids = [k for k in range(1, lab.max() + 1) if (lab == k).sum() >= DMINVOX]
        comp_ids.sort(key=lambda k: -(lab == k).sum())
        for ci, k in enumerate(comp_ids[:8], 1):
            idx = lab == k; nvox = int(idx.sum())
            De = (DE[:, idx] * sqw[:, None]).T; Do = (DO[:, idx] * sqw[:, None]).T
            U, S, Vt = np.linalg.svd(De, full_matrices=False)
            ne = neff_of(idx)
            rs = []; share = S**2 / max((S**2).sum(), 1e-30)
            for kk in range(min(len(S), 4)):
                me = De @ Vt[kk]; mo = Do @ Vt[kk]
                rs.append(np.corrcoef(me, mo)[0, 1])
            zs = [math.atanh(min(max(r, -0.999999), 0.999999)) * math.sqrt(max(ne - 3, 1)) for r in rs]
            nfac = sum(1 for zz in zs if zz >= a.z)
            # singular-value noise floor: pure-noise matrices from the half differences on the same voxels
            smax = []
            for r in range(len(sur)):
                perm = rng.permutation(K)
                Nn = (Dc[perm][:, idx] * sqw[:, None]).T / math.sqrt(2.)
                smax.append(np.linalg.svd(Nn, compute_uv=False)[0])
            floor = np.mean(smax) + 2 * np.std(smax)
            nmag = int((S > floor).sum())
            vd = next((v for v in verdicts if v[0] == ci), None)
            print('    component %d (%d voxels, n_eff %.0f): factor cross-half r %s, Fisher z %s, variance shares %s -> %d reproducible factor(s); singular values %s vs noise floor %.3g -> %d above%s'
                  % (ci, nvox, ne, ' '.join('%.3f' % r for r in rs), ' '.join('%.1f' % zz for zz in zs), ' '.join('%.2f' % s for s in share[:4]), nfac,
                     ' '.join('%.3g' % x for x in S[:4]), floor, nmag,
                     ('; 2-means: ' + vd[1] + ' ' + vd[2]) if (vd and tag == 'old') else ''))

if __name__ == '__main__':
    main()
