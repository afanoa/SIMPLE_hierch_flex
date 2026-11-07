#!/usr/bin/env python3
"""shrinkB.py SRC DST MODE --counts N1,N2,... [--mskdiam D]
Stage-B references between single-iteration refine3D calls. Reads the K state maps the iteration wrote in SRC
(recvol_stateNN{,_even,_odd}.mrc: the shipped pair the next iteration aligns against) and writes them to DST under the same
names, every other file of the state's stem copied unchanged:
  copy  the maps unchanged (the control)
  dev1  each state's deviation from the count-weighted node mean shrunk per Fourier shell by the multivariate Wiener of the
        deviations: signal covariance = their cross-half covariance (PSD), noise = A diag(s2_half) A^T with A = I - 1 w^T
  dev2  the same in two parts, inside / outside tree_wiener's focus mask (20 A between-state cross-half variance), as its emp2
  dev2h dev2 at half dose: each deviation moved half-way to its Wiener estimate
The node mean is kept; each half is shrunk against half-map noise (tree_wiener's W2); merged = (even + odd) / 2."""
import argparse, glob, os, shutil, struct
import numpy as np
from scipy.ndimage import distance_transform_edt, gaussian_filter

LPMASK, MTHR, MSIG, MDIL, EDGE = 20., 0.15, 2., 3., 6.
BANDS = [(1e9, 20.), (20., 14.), (14., 10.), (10., 0.)]


def rdmrc(p):
    raw = open(p, 'rb').read()
    h = np.frombuffer(raw[:1024], np.int32)
    n, mode = int(h[0]), int(h[3])
    if mode != 2 or int(h[1]) != n or int(h[2]) != n:
        raise SystemExit('%s: not a cubic float32 map' % p)
    off = 1024 + int(h[23])
    x = np.frombuffer(raw[off:off + 4 * n ** 3], np.float32).reshape(n, n, n).astype(np.float64)
    smpd = float(np.frombuffer(raw[40:44], np.float32)[0]) / int(h[7])
    return x, raw[:off], smpd


def wrmrc(p, hdr, x):
    h = bytearray(hdr); x = x.astype(np.float32)
    for word, v in ((19, x.min()), (20, x.max()), (21, x.mean()), (54, x.std())):
        struct.pack_into('<f', h, 4 * word, float(v))
    with open(p, 'wb') as f:
        f.write(bytes(h)); f.write(x.tobytes())


def psd(c):
    lam, v = np.linalg.eigh(c)
    return (v * np.maximum(lam, 0.)[None]) @ v.T


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('src'); ap.add_argument('dst'); ap.add_argument('mode', choices=['copy', 'dev1', 'dev2', 'dev2h'])
    ap.add_argument('--counts', required=True); ap.add_argument('--mskdiam', type=float, default=0.)
    a = ap.parse_args()
    cnt = np.array([float(c) for c in a.counts.split(',')]); K = len(cnt)
    os.makedirs(a.dst, exist_ok=True)
    stems = ['recvol_state%02d' % (k + 1) for k in range(K)]
    maps = lambda k: [stems[k] + '.mrc', stems[k] + '_even.mrc', stems[k] + '_odd.mrc']
    for k in range(K):
        for f in glob.glob(os.path.join(a.src, stems[k] + '*')):
            if os.path.isfile(f) and (a.mode == 'copy' or os.path.basename(f) not in maps(k)):
                shutil.copy2(f, os.path.join(a.dst, os.path.basename(f)))
    if a.mode == 'copy':
        print('copy: %d states, counts %s' % (K, cnt.astype(int).tolist())); return
    E, O, HE, HO, HM = [], [], [], [], []
    for k in range(K):
        e, he, smpd = rdmrc(os.path.join(a.src, maps(k)[1])); o, ho, _ = rdmrc(os.path.join(a.src, maps(k)[2]))
        _, hm, _ = rdmrc(os.path.join(a.src, maps(k)[0]))
        E.append(e); O.append(o); HE.append(he); HO.append(ho); HM.append(hm)
    n = E[0].shape[0]; nsh = n // 2; w = cnt / cnt.sum()
    f = np.fft.fftfreq(n) * n; kz = np.arange(n // 2 + 1)
    kx, ky, kzz = np.meshgrid(f, f, kz, indexing='ij'); rad = np.sqrt(kx ** 2 + ky ** 2 + kzz ** 2)
    s3 = np.rint(rad).astype(int); wg3 = np.where((kzz == 0) | (kzz == n // 2), 1., 2.)
    sel = (s3 >= 1) & (s3 < nsh); shell = s3[sel]; wg = wg3[sel]
    order = np.argsort(shell, kind='stable'); so = shell[order]
    lo = np.searchsorted(so, np.arange(nsh)); hi = np.searchsorted(so, np.arange(nsh), 'right')
    rsh = np.where(np.arange(nsh) > 0, n * smpd / np.maximum(np.arange(nsh), 1), np.inf)
    dose = 0.5 if a.mode == 'dev2h' else 1.
    if a.mode.startswith('dev2'):
        f2 = (kx ** 2 + ky ** 2 + kzz ** 2) / float(n * n)
        lpg = lambda v: np.fft.irfftn(np.fft.rfftn(v) * np.exp(-f2 / smpd ** 2 * LPMASK ** 2 * np.log(2.)), s=v.shape)
        c = n // 2; zz, yy, xx = np.mgrid[:n, :n, :n]; rr = np.sqrt((xx - c) ** 2 + (yy - c) ** 2 + (zz - c) ** 2)
        r0 = a.mskdiam / 2. / smpd if a.mskdiam > 0 else n / 2. - EDGE
        msk = np.where(rr <= r0, 1., np.where(rr <= r0 + EDGE, 0.5 + 0.5 * np.cos(np.pi * (rr - r0) / EDGE), 0.)); inm = msk > 0.5
        le = [lpg(x) for x in E]; lo_ = [lpg(x) for x in O]
        eb = sum(w[k] * le[k] for k in range(K)); ob = sum(w[k] * lo_[k] for k in range(K))
        vv = gaussian_filter(sum(w[k] * le[k] * lo_[k] for k in range(K)) - eb * ob, MSIG)
        core = (vv > MTHR * vv[inm].max()) & inm
        mv = np.clip(gaussian_filter((distance_transform_edt(~core) <= MDIL).astype(float), MSIG), 0., 1.) * msk
        parts = [mv, 1. - mv]
        print('focus mask %.1f %% of the sphere (mskdiam %.1f A, r0 %.1f px)' % (100. * mv.sum() / msk.sum(), a.mskdiam, r0))
    else:
        parts = [None]
    A = np.eye(K) - np.outer(np.ones(K), w)
    ce = np.stack([np.fft.rfftn(x)[sel] for x in E]); co = np.stack([np.fft.rfftn(x)[sel] for x in O])
    oute, outo = ce.copy(), co.copy()
    num = np.zeros((len(parts), K, len(BANDS))); den = np.zeros_like(num); xh = np.zeros_like(num); ae = np.zeros_like(num); ao = np.zeros_like(num)
    bandof = lambda s: next(i for i, (hi_, lo2) in enumerate(BANDS) if rsh[s] <= hi_ and rsh[s] > lo2)
    for ip, pm in enumerate(parts):
        if pm is None:
            pe, po = ce, co
        else:
            pe = np.stack([np.fft.rfftn(x * pm)[sel] for x in E]); po = np.stack([np.fft.rfftn(x * pm)[sel] for x in O])
        de = pe - w @ pe; do = po - w @ po
        for s in range(1, nsh):
            idx = order[lo[s]:hi[s]]; g = wg[idx]; sg = g.sum()
            dE, dO = de[:, idx], do[:, idx]
            G = np.real((dE.conj() * g) @ dO.T) / sg
            C = psd(0.5 * (G + G.T))
            s2 = (g * np.abs(pe[:, idx] - po[:, idx]) ** 2).sum(1) / (2. * sg)
            W = C @ np.linalg.pinv(C + A @ np.diag(s2) @ A.T, rcond=1e-8, hermitian=True)
            he, ho = dE + dose * (W @ dE - dE), dO + dose * (W @ dO - dO)
            oute[:, idx] += he - dE; outo[:, idx] += ho - dO
            b = bandof(s)
            num[ip, :, b] += (g * (np.abs(he) ** 2 + np.abs(ho) ** 2)).sum(1); den[ip, :, b] += (g * (np.abs(dE) ** 2 + np.abs(dO) ** 2)).sum(1)
            xh[ip, :, b] += (g * np.real(dE.conj() * dO)).sum(1); ae[ip, :, b] += (g * np.abs(dE) ** 2).sum(1); ao[ip, :, b] += (g * np.abs(dO) ** 2).sum(1)
    hdr = '%-8s %6s' % ('part', 'state') + ''.join('   %s' % ('%g-%gA' % (min(h_, 999), l_)).rjust(14) for h_, l_ in BANDS)
    print('%s: K %d, counts %s, box %d at %.3f A/px; per band: kept power of the deviation / deviation cross-half FSC' % (a.mode, K, cnt.astype(int).tolist(), n, smpd))
    print(hdr)
    for ip in range(len(parts)):
        for k in range(K):
            cells = ''.join('   %6.3f/%6.3f' % (num[ip, k, b] / max(den[ip, k, b], 1e-300), xh[ip, k, b] / max(np.sqrt(ae[ip, k, b] * ao[ip, k, b]), 1e-300)) for b in range(len(BANDS)))
            print('%-8s %6d%s' % (['all', 'inside', 'outside'][0 if len(parts) == 1 else ip + 1], k + 1, cells))
    for k in range(K):
        Fe = np.zeros((n, n, n // 2 + 1), complex); Fo = np.zeros_like(Fe)
        Fe0 = np.fft.rfftn(E[k]); Fo0 = np.fft.rfftn(O[k])
        Fe[:] = Fe0; Fo[:] = Fo0; Fe[sel] = oute[k]; Fo[sel] = outo[k]
        e = np.fft.irfftn(Fe, s=(n, n, n)); o = np.fft.irfftn(Fo, s=(n, n, n))
        wrmrc(os.path.join(a.dst, maps(k)[1]), HE[k], e); wrmrc(os.path.join(a.dst, maps(k)[2]), HO[k], o)
        wrmrc(os.path.join(a.dst, maps(k)[0]), HM[k], 0.5 * (e + o))
        print('state %d: rms change even %.4g of rms %.4g' % (k + 1, np.sqrt(((e - E[k]) ** 2).mean()), np.sqrt((E[k] ** 2).mean())))


if __name__ == '__main__':
    main()
