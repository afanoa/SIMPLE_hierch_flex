#!/usr/bin/env python3
"""usage: bshare.py CASES < QSCREEN   (quality_screen3.py's case list; L from its QSCREEN lines)
Reference-free quality test: is the end-bin difference D = M_L - M_H a radially filtered copy of the two bins' mean map M (the same
molecule seen blurrier or dimmer: particle quality, dilution) or a localised change (a real state)? From the unfiltered halves in the
case's soft mask, per shell from 40 A to max(6 A, 2.5 px): P_DD = Re<D_e, D_o> (D's reproducible power), P_MM = Re<M_e, M_o> and
P_DM = (Re<D_e, M_o> + Re<D_o, M_e>) / 2.
- share = sum P_DM^2 / P_MM over sum P_DD: the share of D's reproducible power that a per-shell rescaling of M explains (1 = a pure
  filter of the mean map; 0 = orthogonal to it);
- c(s) = P_DM / P_MM fitted against s^2 (weights P_MM): c0 = the scale difference, dBm = -4 x slope = the blur of L vs H (A^2, > 0 = L
  blurrier) under D / M = -dB s^2 / 4;
- rel = sum P_DD / sum P_MM: the size of the change.
One BS line per case."""
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
LB = {}
for l in sys.stdin:
    m = re.match(r'QSCREEN (\S+)\s+\S+\s+\| L = ([AB])', l)
    if m: LB[m.group(1)] = m.group(2) == 'B'
for l in open(sys.argv[1]):
    w = l.split()
    if len(w) < 11 or w[0].startswith('#') or w[0] not in LB: continue
    hv = [unfil(p) for p in (w[3], w[4], w[6], w[7])]
    if None in hv: print('BS %-36s | no unfiltered halves' % w[0]); continue
    smpd = read_mrc(w[2])[1]; EA, OA, EB, OB = (read_mrc(p)[0] for p in hv); n = EA.shape[0]
    EL, OL, EH, OH = (EB, OB, EA, OA) if LB[w[0]] else (EA, OA, EB, OB); sm = soft(fresample(read_mrc(w[8])[0], n) > 0.5, n)
    FEL, FOL, FEH, FOH = (np.fft.fftn(x * sm) for x in (EL, OL, EH, OH)); De, Do, Me, Mo = FEL - FEH, FOL - FOH, (FEL + FEH) / 2, (FOL + FOH) / 2
    sh = np.rint(kr(n)).astype(int); k0, k1 = max(1, int(np.ceil(n * smpd / 40.))), int(n * smpd / max(6., 2.5 * smpd)); sel = (sh >= k0) & (sh <= k1); s = sh[sel]
    P = {nm: np.bincount(s, v, k1 + 1)[k0:] for nm, v in (('DD', (De[sel] * Do[sel].conj()).real), ('MM', (Me[sel] * Mo[sel].conj()).real),
         ('DM', ((De[sel] * Mo[sel].conj()).real + (Do[sel] * Me[sel].conj()).real) / 2))}
    ok = P['MM'] > 0; s2 = (np.arange(k0, k1 + 1) / (n * smpd)) ** 2; c = P['DM'][ok] / P['MM'][ok]
    share = (P['DM'][ok] ** 2 / P['MM'][ok]).sum() / P['DD'].sum(); c1, c0 = np.polyfit(s2[ok], c, 1, w=np.sqrt(P['MM'][ok]))
    print('BS %-36s %-8s | %d shells | share %.3f | c0 %+.4f dBm %+.1f A^2 | rel %.4f' % (w[0], w[1], ok.sum(), share, c0, -4 * c1, P['DD'].sum() / P['MM'].sum())); sys.stdout.flush()
