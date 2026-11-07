#!/usr/bin/env python3
"""After the masked split test of one board: which original sub-states ended in the same final state, the range of motion
kept (largest masked distance between surviving states over the largest between all tested sub-states; maps = the f692e
state maps, population-weighted within a surviving state, 12 A low-pass, standardised inside the mask) and the closed-arm
occupancy of each sub-state (arm_readout output). usage: prune_post.py OUTDIR STATES_DIR MASK ARM_TXT NAME"""
import sys, os, numpy as np
out, sdir, maskf, armf, name = sys.argv[1:6]
lab = np.loadtxt(out + '/labels.txt', comments='#', ndmin=2).astype(int); small = np.loadtxt(out + '/small.txt', comments='#', ndmin=2).astype(int)
def rd(f):
    h = np.fromfile(f, np.int32, 256); nx, ny, nz = h[:3]; return np.fromfile(f, np.float32, offset=1024 + int(h[23]), count=int(nx * ny * nz)).reshape(nz, ny, nx).astype(np.float64)
groups = {}
if os.path.isfile(out + '/children_A.txt') and len(lab):
    c, ch = np.loadtxt(out + '/cand_A.txt').astype(int), np.loadtxt(out + '/children_A.txt').astype(int)
    for cand, orig in lab[:, :2]:
        m = c == cand; g = int(np.bincount(ch[m]).argmax()) if m.any() else 0; groups.setdefault(g, []).append(int(orig))
arm = {}
if os.path.isfile(armf):
    for l in open(armf):
        if l.startswith('#') or not l.strip(): continue
        p = l.split(); b = os.path.basename(p[0])
        if b.startswith('state_'): arm[int(b[6:8])] = float(p[1])
M = rd(maskf) if maskf != 'none' and os.path.isfile(maskf) else None; B = None; maps = {}
for orig in lab[:, 1] if len(lab) else []:
    f = '%s/state_%02d.mrc' % (sdir, orig)
    if os.path.isfile(f): maps[int(orig)] = rd(f)
rng = float('nan')
if len(maps) >= 2:
    n = next(iter(maps.values())).shape[0]; fr = np.fft.fftfreq(n, d=1.464); kx, ky, kz = np.meshgrid(fr, fr, fr, indexing='ij'); G = np.exp(-0.5 * (kx**2 + ky**2 + kz**2) * (12. / 1.2) ** 2)
    m = (M > 0.5) if M is not None else np.ones((n, n, n), bool)
    def prep(v): v = np.real(np.fft.ifftn(np.fft.fftn(v) * G))[m]; v = v - v.mean(); return v / max(v.std(), 1e-12)
    P = {k: prep(v) for k, v in maps.items()}; w = {int(o): int(a + b) for _, o, a, b in lab}
    def dmax(vs): return max(np.sqrt(np.mean((vs[i] - vs[j]) ** 2)) for i in range(len(vs)) for j in range(i + 1, len(vs))) if len(vs) > 1 else 0.
    allv = list(P.values()); surv = []
    for g, mem in groups.items():
        mem = [k for k in mem if k in P]
        if mem: v = sum(w[k] * P[k] for k in mem) / sum(w[k] for k in mem); surv.append(v / max(v.std(), 1e-12))
    rng = dmax(surv) / max(dmax(allv), 1e-12)
closed = [k for k, v in arm.items() if v >= 0.8]
def gtxt(): return ' '.join('[%s]' % ','.join('%d%s' % (k, '*' if k in closed else '') for k in sorted(mem)) for g, mem in sorted(groups.items()))
kst = len(groups); ntest = len(lab); nsmall = len(small)
csurv = [k for k in closed if any(k in mem and len([x for x in mem if x not in closed]) == 0 for mem in groups.values())]
print('PRUNE %-34s sub-states %3d | too small %3d | tested %3d -> supported %3d | range kept %.2f | closed-arm>=0.8: %s (kept as own state: %s) | groups %s'
      % (name, ntest + nsmall, nsmall, ntest, kst, rng, ','.join(map(str, sorted(closed))) or '-', ','.join(map(str, sorted(csurv))) or '-', gtxt()))
