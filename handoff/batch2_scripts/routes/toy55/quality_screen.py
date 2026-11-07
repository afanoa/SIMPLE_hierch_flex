#!/usr/bin/env python3
"""usage: quality_screen.py CASES
A knob-free quality screen for a motion candidate or a split: is the lower-corr group worse particles, or another state?
CASES lines: name truth mapA evenA oddA mapB evenB oddB mask dump.npz groups.txt
(groups.txt: make_groups.py's 1 = A / 2 = B / 0 per project row; dump.npz: dump_ptcl.py's corr / dfx per project row).
L is the group with the lower median corr, H the other.
(a) corr AUC of L vs H inside dfx deciles (matched defocus); its null sd from 20 within-decile shuffles of the labels; z_a.
(b) the map signature of D = M_L - M_H (soft mask, shells to 15 A): its signal cosine to the two groups' mean map, and the share of the
    mask where the 15 A low-passed D exceeds 3 sd of its half-map noise.
(c) the stable region: the mask where the low-passed abs(D) stays under 3 sd of its noise (the part of the molecule that does not
    change; the 30 % of the mask with the smallest abs(D) if that is under 10 %). Per shell from 15 to 6 A, each group's full-map SSNR
    from its halves' FSC there (2 FSC / (1 - FSC), FSC in 0.02-0.98) per particle; dq = the mean over shells of
    log(SSNR_L / n_L) - log(SSNR_H / n_H), SE from the shell-to-shell scatter, z_c = dq / SE.
Worse particles are worse everywhere, so the low-corr group's stable region loses SSNR per particle (z_c << 0). A real state differs
from the consensus only where it changes, so its stable region is as good (z_c ~ 0). One QSCREEN line per case."""
import sys, numpy as np
def read_mrc(f):
    raw = open(f, 'rb').read(); h = np.frombuffer(raw[:1024], np.int32); hf = h.view(np.float32); n = int(h[0])
    return np.frombuffer(raw[1024 + int(h[23]):], np.float32, count=n ** 3).astype(np.float64).reshape(n, n, n), float(hf[10]) / n
def fresample(v, n):
    m = v.shape[0]
    if m == n: return v
    F = np.fft.fftshift(np.fft.fftn(v)); c = m // 2
    if n < m: s = c - n // 2; G = F[s:s + n, s:s + n, s:s + n]
    else: G = np.zeros((n, n, n), complex); s = n // 2 - c; G[s:s + m, s:s + m, s:s + m] = F
    return np.real(np.fft.ifftn(np.fft.ifftshift(G))) * (n / m) ** 3
def kr(n): f = np.fft.fftfreq(n) * n; return np.sqrt(f[:, None, None] ** 2 + f[None, :, None] ** 2 + f[None, None, :] ** 2)
def cdot(x, y): return float((x.conj() * y).real.sum())
def soft(m, n): return np.real(np.fft.ifftn(np.fft.fftn(m.astype(float)) * np.exp(-0.5 * (2. * 2 * np.pi * kr(n) / n) ** 2)))
def auc(x, y):
    x, y = x[np.isfinite(x)], y[np.isfinite(y)]
    if len(x) < 20 or len(y) < 20: return np.nan
    r = np.argsort(np.argsort(np.r_[x, y])); return (r[:len(x)].sum() - len(x) * (len(x) - 1) / 2.) / (len(x) * len(y))
DC = {}
for l in open(sys.argv[1]):
    w = l.split()
    if len(w) < 11 or w[0].startswith('#'): continue
    nm, truth = w[0], w[1]; MA, MB = [read_mrc(f)[0] for f in w[2:5]], [read_mrc(f)[0] for f in w[5:8]]; smpd = read_mrc(w[2])[1]; n = MA[0].shape[0]
    if w[9] not in DC: DC[w[9]] = np.load(w[9])
    D = DC[w[9]]; corr, dfx = D['corr'], D['dfx']; grp = np.loadtxt(w[10], dtype=int); ok = np.isfinite(dfx)
    dd = np.full(len(dfx), -1); dd[ok] = np.searchsorted(np.quantile(dfx[ok], np.linspace(0.1, 0.9, 9)), dfx[ok])
    ia, ib = grp == 1, grp == 2; Lb = np.nanmedian(corr[ib]) < np.nanmedian(corr[ia]); iL, iH = (ib, ia) if Lb else (ia, ib); ML, MH = (MB, MA) if Lb else (MA, MB)
    a_ = np.nanmean([auc(corr[iL & (dd == j)], corr[iH & (dd == j)]) for j in range(10)]); rng = np.random.default_rng(3); nul = []
    for _ in range(20):
        lab = iL.copy(); both = iL | iH
        for j in range(10):
            q = np.nonzero(both & (dd == j))[0]; lab[q] = rng.permutation(lab[q])
        nul.append(np.nanmean([auc(corr[lab & (dd == j)], corr[both & ~lab & (dd == j)]) for j in range(10)]))
    za = (a_ - 0.5) / max(np.std(nul), 1e-6)
    mk = fresample(read_mrc(w[8])[0], n) > 0.5; sh = np.rint(kr(n)).astype(int); res = 15.
    lpw = np.clip((1.1 * smpd / res - kr(n) / n) / (0.2 * smpd / res), 0., 1.); lpw = 0.5 - 0.5 * np.cos(np.pi * lpw)
    FL, FH = [np.fft.fftn(x) for x in ML], [np.fft.fftn(x) for x in MH]
    Dr = np.real(np.fft.ifftn((FL[0] - FH[0]) * lpw)); Nr = np.real(np.fft.ifftn(((FL[1] - FL[2]) - (FH[1] - FH[2])) * lpw)) / 2; sd = Nr[mk].std()
    share = (np.abs(Dr[mk]) > 3 * sd).mean(); stab = mk & (np.abs(Dr) <= 3 * sd)
    if stab.sum() < 0.1 * mk.sum(): stab = mk & (np.abs(Dr) <= np.percentile(np.abs(Dr[mk]), 30))
    sm = soft(mk, n); sel = (sh >= 1) & (sh <= int(n * smpd / res))
    QL, QH = [np.fft.fftn(x * sm)[sel] for x in ML], [np.fft.fftn(x * sm)[sel] for x in MH]
    dE, dO = QL[1] - QH[1], QL[2] - QH[2]; pE, pO = (QL[1] + QH[1]) / 2, (QL[2] + QH[2]) / 2; rep = cdot(dE, dO)
    cs = cdot(dE, pO) / np.sqrt(max(rep, 1e-30) * max(cdot(pE, pO), 1e-30)) if rep > 0 else np.nan
    ss = soft(stab, n); s1, s2 = int(np.ceil(n * smpd / 15.)), int(n * smpd / 6.); lr = []
    EL, OL, EH, OH = (np.fft.fftn(x * ss) for x in (ML[1], ML[2], MH[1], MH[2])); nL, nH = iL.sum(), iH.sum()
    for k in range(s1, s2 + 1):
        q = sh == k; f = []
        for E_, O_ in ((EL, OL), (EH, OH)): f.append(float((E_[q] * O_[q].conj()).real.sum() / np.sqrt((np.abs(E_[q]) ** 2).sum() * (np.abs(O_[q]) ** 2).sum())))
        if all(0.02 < x < 0.98 for x in f): lr.append(np.log(2 * f[0] / (1 - f[0]) / nL) - np.log(2 * f[1] / (1 - f[1]) / nH))
    lr = np.array(lr); dq = lr.mean() if len(lr) else np.nan; se = lr.std(ddof=1) / np.sqrt(len(lr)) if len(lr) > 2 else np.nan; zc = dq / se if se > 0 else np.nan
    print('QSCREEN %-24s %-8s | L = %s (%d) vs H (%d) | (a) AUC %.3f, null sd %.4f, z_a %+.1f | (b) cos %+.3f, share %.3f | (c) stable %.2f of mask, %d shells, dq %+.3f, SE %.3f, z_c %+.1f'
          % (nm, truth, 'B' if Lb else 'A', nL, nH, a_, np.std(nul), za, cs, share, stab.sum() / mk.sum(), len(lr), dq, se, zc)); sys.stdout.flush()
