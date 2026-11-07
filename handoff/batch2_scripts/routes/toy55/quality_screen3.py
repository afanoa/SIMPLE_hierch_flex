#!/usr/bin/env python3
"""usage: quality_screen3.py CASES   (quality_screen.py's case list; L is read from its QSCREEN lines on stdin)
quality_screen2.py's blur and quality_screen.py's SSNR per particle from the UNFILTERED halves' cross-power, so a map's own filter or
regulariser cannot enter (the full maps and _even / _odd halves can be FSC-weighted or regularised per state):
- halves: X_even_unfil.mrc / X_odd_unfil.mrc next to each listed half (a series' halves/binK_even.mrc maps to work/recvol_state0K_even_unfil.mrc);
- stable region as before (mask voxels with the 15 A low-passed abs(D) under 3 sd of its half-map noise, D from the listed full maps);
- per shell from 15 A to max(6 A, 2.5 px): the signal power Re<E, O> of each group in the soft stable region;
  dBx = -4 x the slope of 0.5 log(S_L / S_H) against s^2 (A^2, > 0 = L blurrier), SE from the fit;
  dqx = the mean over shells of log((S_L / N_L^2) ... no: of log(SSNR_L / n_L) - log(SSNR_H / n_H) with SSNR = S / noise power
  (noise power |E - O|^2 / 4 of the same shell), SE from the shell scatter.
One QS3 line per case."""
import os, re, sys, numpy as np
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
def soft(m, n): return np.real(np.fft.ifftn(np.fft.fftn(m.astype(float)) * np.exp(-0.5 * (2. * 2 * np.pi * kr(n) / n) ** 2)))
def unfil(p):
    m = re.search(r'^(.*)/halves/bin(\d)_(even|odd)\.mrc$', p)
    q = '%s/work/recvol_state%02d_%s_unfil.mrc' % (m.group(1), int(m.group(2)), m.group(3)) if m else (p if p.endswith('_unfil.mrc') else p[:-4] + '_unfil.mrc')
    return q if os.path.exists(q) else None
LB, NN = {}, {}
for l in sys.stdin:
    m = re.match(r'QSCREEN (\S+)\s+\S+\s+\| L = ([AB]) \((\d+)\) vs H \((\d+)\)', l)
    if m: LB[m.group(1)] = m.group(2) == 'B'; NN[m.group(1)] = (int(m.group(3)), int(m.group(4)))
for l in open(sys.argv[1]):
    w = l.split()
    if len(w) < 11 or w[0].startswith('#') or w[0] not in LB: continue
    hv = [unfil(p) for p in (w[3], w[4], w[6], w[7])]
    if None in hv: print('QS3 %-36s %-8s | no unfiltered halves' % (w[0], w[1])); continue
    A, B = read_mrc(w[2])[0], read_mrc(w[5])[0]; smpd = read_mrc(w[2])[1]; n = A.shape[0]; EA, OA, EB, OB = (read_mrc(p)[0] for p in hv)
    (ML, EL, OL, MH, EH, OH) = (B, EB, OB, A, EA, OA) if LB[w[0]] else (A, EA, OA, B, EB, OB); nL, nH = NN[w[0]]
    mk = fresample(read_mrc(w[8])[0], n) > 0.5; K = kr(n); sh = np.rint(K).astype(int); res = 15.
    lpw = np.clip((1.1 * smpd / res - K / n) / (0.2 * smpd / res), 0., 1.); lpw = 0.5 - 0.5 * np.cos(np.pi * lpw)
    Dr = np.real(np.fft.ifftn((np.fft.fftn(ML) - np.fft.fftn(MH)) * lpw)); Nr = np.real(np.fft.ifftn((np.fft.fftn(EL - OL) - np.fft.fftn(EH - OH)) * lpw)) / 2
    stab = mk & (np.abs(Dr) <= 3 * Nr[mk].std())
    if stab.sum() < 0.1 * mk.sum(): stab = mk & (np.abs(Dr) <= np.percentile(np.abs(Dr[mk]), 30))
    ss = soft(stab, n); FEL, FOL, FEH, FOH = (np.fft.fftn(x * ss) for x in (EL, OL, EH, OH)); x, y, lq = [], [], []
    for k in range(int(np.ceil(n * smpd / 15.)), int(n * smpd / max(6., 2.5 * smpd)) + 1):
        q = sh == k; SL = (FEL[q] * FOL[q].conj()).real.sum(); SH = (FEH[q] * FOH[q].conj()).real.sum()
        NL = (np.abs(FEL[q] - FOL[q]) ** 2).sum() / 4; NH = (np.abs(FEH[q] - FOH[q]) ** 2).sum() / 4
        if SL > 0 and SH > 0: x.append((k / (n * smpd)) ** 2); y.append(0.5 * np.log(SL / SH)); lq.append(np.log(SL / NL / nL) - np.log(SH / NH / nH))
    x, y, lq = np.array(x), np.array(y), np.array(lq)
    if len(x) < 4: print('QS3 %-36s %-8s | too few shells' % (w[0], w[1])); continue
    c = np.polyfit(x, y, 1); r = y - np.polyval(c, x); se = np.sqrt(r @ r / (len(x) - 2) / ((x - x.mean()) ** 2).sum())
    print('QS3 %-36s %-8s | %d shells | dBx %+.1f A^2 (SE %.1f) | dqx %+.3f (SE %.3f, z %+.1f)' % (w[0], w[1], len(x), -4 * c[0], 4 * se, lq.mean(), lq.std(ddof=1) / np.sqrt(len(lq)),
          lq.mean() / (lq.std(ddof=1) / np.sqrt(len(lq))))); sys.stdout.flush()
