#!/usr/bin/env python3
"""Map-level truth readout of K child maps (a reconstruct3D dir with recvol_stateNN_{even,odd}_unfil.mrc, states = children in
coordinate order) on a synthetic benchmark: per child the half-map FSC 0.143 resolution and its members' truth (n, swing share,
mean swing); for the two END children (1 and K) the difference map: split-half reproducibility per band (corr of the even and the
odd difference) and agreement with the TRUE difference (members' conformer-weighted truth maps, junk left out), per band; the same
for every pair of neighbouring children (median). Soft sphere of MSKDIAM.
usage: subcls_mapscore.py RECDIR LABELS GT TRUTHMAPDIR [--mskdiam 220] [--smpd 1.464] [--tag NAME]"""
import sys, glob, os, argparse, numpy as np
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/routes/toy55'); from model_class import read_mrc
ap = argparse.ArgumentParser(); ap.add_argument('rec'); ap.add_argument('labels'); ap.add_argument('gt'); ap.add_argument('tmaps')
ap.add_argument('--mskdiam', type=float, default=220.); ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--tag', default='')
a = ap.parse_args()
lab = np.loadtxt(a.labels).astype(int); G = np.loadtxt(a.gt, comments='#', ndmin=2); conf = np.zeros(len(lab), int); conf[G[:, 0].astype(int) - 1] = G[:, 1].astype(int)
K = lab.max(); kids = list(range(1, K + 1))
def rd(fn): return np.asarray(read_mrc(fn), np.float64)
v0 = rd(os.path.join(a.rec, 'recvol_state01_even_unfil.mrc')); n = v0.shape[0]; c = n // 2
zz, yy, xx = np.mgrid[:n, :n, :n]; rr = np.sqrt((xx - c) ** 2. + (yy - c) ** 2. + (zz - c) ** 2.); R = a.mskdiam / 2. / a.smpd; msk = np.clip((R + 3 - rr) / 6., 0, 1)
k = np.fft.fftfreq(n) * n; KZ, KY, KX = np.meshgrid(k, k, k, indexing='ij'); sh = np.sqrt(KX ** 2 + KY ** 2 + KZ ** 2).round().astype(int).ravel(); nsh = n // 2
def ft(v): return np.fft.fftn(v * msk).ravel()
def fsc(A, B):
    m = sh < nsh; num = np.bincount(sh[m], (A[m] * B[m].conj()).real, nsh); da = np.bincount(sh[m], np.abs(A[m]) ** 2, nsh); db = np.bincount(sh[m], np.abs(B[m]) ** 2, nsh)
    return num / np.sqrt(np.maximum(da * db, 1e-30))
def res143(f):
    below = np.where(f[1:] < 0.143)[0]; s = below[0] + 1 if len(below) else nsh - 1; return n * a.smpd / max(s, 1)
bands = [(117., 30.), (30., 20.), (20., 12.), (12., 8.)]
def bandcorr(A, B):
    out = []
    for lo, hi in bands:
        m = (sh >= n * a.smpd / lo) & (sh <= n * a.smpd / hi); out.append((A[m] * B[m].conj()).real.sum() / np.sqrt((np.abs(A[m]) ** 2).sum() * (np.abs(B[m]) ** 2).sum()))
    return np.array(out)
tm = {}
def tmap(cid):
    if cid not in tm:
        f = glob.glob(os.path.join(a.tmaps, 'conf%02d_*_m160.mrc' % cid)) if cid < 100 else glob.glob(os.path.join(a.tmaps, 'conf%d_*_m160.mrc' % cid))
        tm[cid] = ft(rd(f[0])) if f else None
    return tm[cid]
E, H = {}, {}
print('%schild: members, swing share, mean swing (A), half-map FSC 0.143 (A)' % ((a.tag + ' ') if a.tag else ''))
for kk in kids:
    m = lab == kk; e = ft(rd(os.path.join(a.rec, 'recvol_state%02d_even_unfil.mrc' % kk))); o = ft(rd(os.path.join(a.rec, 'recvol_state%02d_odd_unfil.mrc' % kk))); H[kk] = (e, o)
    cs, cn = np.unique(conf[m & (conf > 0)], return_counts=True); acc = 0; wsum = 0
    for ci, ni in zip(cs, cn):
        t = tmap(int(ci))
        if t is not None: acc = acc + ni * t; wsum += ni
    E[kk] = acc / max(wsum, 1); sw = (conf[m] >= 200) & (conf[m] <= 230)
    print('  %2d: %6d  %.2f  %5.1f  %5.1f' % (kk, m.sum(), sw.mean(), (conf[m][sw] - 200).mean() if sw.any() else np.nan, res143(fsc(e, o))))
def diff(i, j):
    De = H[j][0] - H[i][0]; Do = H[j][1] - H[i][1]; return bandcorr(De, Do), bandcorr(De + Do, E[j] - E[i])
rep, fid = diff(kids[0], kids[-1]); nb = [diff(kids[i], kids[i + 1]) for i in range(len(kids) - 1)]
fmt = lambda v: ' '.join('%.2f' % x for x in v)
print('%send children %d vs %d, difference map per band 117-30 / 30-20 / 20-12 / 12-8 A: split-half reproducibility %s | agreement with the true difference %s'
      % ((a.tag + ' ') if a.tag else '', kids[0], kids[-1], fmt(rep), fmt(fid)))
print('%sneighbouring children (median over %d pairs): reproducibility %s | agreement with truth %s'
      % ((a.tag + ' ') if a.tag else '', len(nb), fmt(np.median([x[0] for x in nb], 0)), fmt(np.median([x[1] for x in nb], 0))))
