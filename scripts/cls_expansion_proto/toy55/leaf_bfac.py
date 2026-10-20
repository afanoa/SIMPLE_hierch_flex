#!/usr/bin/env python3
"""Do the leaf states of a klem tree differ by sharpness rather than structure? Per leaf (fam?/<LEAFGLOB>/level1/maps/
recvol_state01{,_even,_odd}.mrc): particles, FSC 0.143 resolution of its halves, relative B-factor vs the consensus
(Guinier slope of the amplitude ratio over 20..8 A inside the envelope), and the correlation of the leaf's deviation from
the consensus BEFORE and AFTER per-shell amplitude matching with every other leaf's. usage: leaf_bfac.py RUN CONSENSUS"""
import sys, os, glob, re, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc
run, consf = sys.argv[1:3]; smpd = float(os.environ.get('SMPD', 1.464)); lo, hi = float(os.environ.get('BLO', 20.)), float(os.environ.get('BHI', 8.))
cons = read_mrc(consf).astype(np.float32); n = cons.shape[0]; x = np.arange(n) - n // 2; rr = np.sqrt(x[:, None, None] ** 2 + x[None, :, None] ** 2 + x[None, None, :] ** 2)
env = np.clip((100. / smpd + 6 - rr) / 12., 0, 1).astype(np.float32); fsh = np.round(np.sqrt(np.sum(np.array(np.meshgrid(*[np.fft.fftfreq(n) * n] * 3, indexing='ij')) ** 2, 0))).astype(int); nsh = n // 2 + 1
klo, khi = int(round(n * smpd / lo)), int(round(n * smpd / hi)); shells = np.arange(nsh); s2 = (shells / (n * smpd)) ** 2
def F(v): return np.fft.fftn(v * env)
cnt = np.bincount(fsh.ravel(), None, minlength=nsh)[:nsh]
def pshell(A): return np.bincount(fsh.ravel(), (np.abs(A) ** 2).ravel(), minlength=nsh)[:nsh] / np.maximum(cnt, 1)
def res143(E, O):
    num = np.bincount(fsh.ravel(), np.real(E * np.conj(O)).ravel(), minlength=nsh)[:nsh] / np.maximum(cnt, 1); fsc = num / np.sqrt(np.maximum(pshell(E) * pshell(O), 1e-30))
    k = np.where(fsc[1:] < 0.143)[0]; return n * smpd / (k[0] + 1) if len(k) else hi
CF = F(cons); PC = pshell(CF); leaves = []
for d in sorted(glob.glob(run + '/fam?/' + os.environ.get('LEAFGLOB', 'region*') + '/level1/maps')):
    m = re.search(r'fam(\d)/(\w+)/level1', d); fe, fo, ff = [d + '/recvol_state01%s.mrc' % s for s in ('_even', '_odd', '')]
    if not (os.path.isfile(fe) and os.path.isfile(fo)): continue
    E, O = F(read_mrc(fe).astype(np.float32)), F(read_mrc(fo).astype(np.float32)); A = 0.5 * (E + O); PA = pshell(A)
    ratio = np.sqrt(np.maximum(PA, 1e-30) / np.maximum(PC, 1e-30)); sel = (shells >= klo) & (shells <= khi); slope = np.polyfit(s2[sel], np.log(ratio[sel]), 1)[0]; B = -4. * slope
    cl = d.replace('/level1/maps', '/states_clusters.txt'); cl = cl.replace('/maps', ''); npart = 0
    for c in (d.replace('/maps', '/states_clusters.txt'),):
        if os.path.isfile(c): npart = int((np.loadtxt(c, comments='#', ndmin=2)[:, 0] == 1).sum())
    w = np.where(PA > 0, np.sqrt(np.maximum(PC, 1e-30) / np.maximum(PA, 1e-30)), 1.)[np.minimum(fsh, nsh - 1)]   # per-shell amplitude match to the consensus
    band = (fsh >= klo) & (fsh <= khi); dev = (A - CF)[band]; devm = (w * A - CF)[band]
    leaves.append(('f%s_%s' % (m.group(1), m.group(2)), npart, res143(E, O), B, dev, devm))
print('%-18s %6s %8s %8s' % ('leaf', 'n', 'res143', 'B_rel')); [print('%-18s %6d %8.2f %8.1f' % l[:4]) for l in leaves]
L = len(leaves); C = np.zeros((L, L)); Cm = np.zeros((L, L))
for i in range(L):
    for j in range(L):
        C[i, j] = np.real(np.vdot(leaves[i][4], leaves[j][4])) / (np.linalg.norm(leaves[i][4]) * np.linalg.norm(leaves[j][4]) + 1e-30); Cm[i, j] = np.real(np.vdot(leaves[i][5], leaves[j][5])) / (np.linalg.norm(leaves[i][5]) * np.linalg.norm(leaves[j][5]) + 1e-30)
off = ~np.eye(L, dtype=bool); print('deviation-from-consensus correlation between leaves (%.0f..%.0f A): raw median %.3f (min %.3f), after per-shell amplitude matching median %.3f (min %.3f)' % (lo, hi, np.median(C[off]), C[off].min(), np.median(Cm[off]), Cm[off].min()))
Bs = np.array([l[3] for l in leaves]); print('relative B-factors: median %.1f, spread (std) %.1f, range %.1f..%.1f' % (np.median(Bs), Bs.std(), Bs.min(), Bs.max()))
