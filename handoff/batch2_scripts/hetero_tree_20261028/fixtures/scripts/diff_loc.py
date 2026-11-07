#!/usr/bin/env python3
"""where two state maps differ: both maps made into 10 A views at box 160 (make_view), low-passed and standardised as in
arm_readout.py; the top 5 % of |B - A| over the density support are counted in arm_readout's closed-arm mask (the deposited closed
conformer's density the consensus lacks), its open-arm mask (the reverse) and the rest, against the masks' share of the support
(1.0 = no enrichment); the signed mean B - A in each mask says which state holds more of that arm.
usage: diff_loc.py OUT_DIR SPEC   (SPEC lines: NAME A_EVEN A_ODD B_EVEN B_ODD)"""
import sys, os, subprocess, numpy as np
from pathlib import Path
OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True); SC = Path(__file__).resolve().parent; PY = sys.executable
sys.path.insert(0, str(SC)); from model_class import read_mrc
R160 = Path('/mnt/beegfs/elmlund/afan/routes/toy55/maps160'); CL = sorted(R160.glob('docked/conf18_*_m160.mrc'))[0]
pairs = [l.split() for l in open(sys.argv[2]) if l.strip() and not l.startswith('#')]
with open(OUT / 'views_spec.txt', 'w') as f:
    for nm, ae, ao, be, bo in pairs: f.write('%s__A %s %s\n%s__B %s %s\n' % (nm, ae, ao, nm, be, bo))
env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
subprocess.run([PY, str(SC / 'make_view.py'), str(OUT / 'views'), str(OUT / 'views_spec.txt'), '--box', '160'], check=True, env=env, stdout=subprocess.DEVNULL)
smpd, LP = 1.464, 10.
cons0 = read_mrc(str(R160 / 'toy_consensus.mrc')); n = cons0.shape[0]; f = np.fft.fftfreq(n, d=smpd)
kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); G = np.exp(-0.5 * (kx**2 + ky**2 + kz**2) * (LP / 1.2)**2)
lp = lambda v: np.real(np.fft.ifftn(np.fft.fftn(v.astype(np.float32)) * G))
cons = lp(cons0); supp = cons > 0.2 * cons.max(); outside = cons < 0.02 * cons.max()
def std(v): v = v - v[outside].mean(); return v / max(np.sqrt(np.mean(v[supp] ** 2)), 1e-12)
cons = std(cons); closed = std(lp(read_mrc(str(CL)))); d = closed - cons
thr = 0.5 * np.percentile(np.abs(d[supp | (closed > 0.2 * closed.max())]), 99); mC = d > thr; mO = d < -thr
print('closed-arm mask %d voxels, open-arm mask %d voxels; enrichment = share of the top-5 %% |B - A| voxels / the mask share of the region'
      % (mC.sum(), mO.sum()))
for nm, *_ in pairs:
    A = std(lp(read_mrc(str(OUT / 'views' / (nm + '__A.mrc'))))); B = std(lp(read_mrc(str(OUT / 'views' / (nm + '__B.mrc')))))
    reg = supp | mC | mO; dd = np.abs(B - A); top = reg & (dd >= np.percentile(dd[reg], 95))
    sC, sO = (top & mC).sum() / top.sum(), (top & mO).sum() / top.sum(); bC, bO = mC[reg].mean(), mO[reg].mean()
    print('%-26s top-5 %%: closed-arm %.2f (x%.1f), open-arm %.2f (x%.1f), elsewhere %.2f | mean B - A: closed-arm %+.3f, open-arm %+.3f, rest %+.3f'
          % (nm, sC, sC / max(bC, 1e-12), sO, sO / max(bO, 1e-12), 1 - sC - sO, (B - A)[mC].mean(), (B - A)[mO].mean(), (B - A)[supp & ~mC & ~mO].mean()))
