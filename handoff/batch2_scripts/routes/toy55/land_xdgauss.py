#!/usr/bin/env python3
"""land_xdgauss.py: control for land_xdcompare.py. From its saved fits (--save NPZ, all mixtures already in the reference frame): the
Bhattacharyya coefficient of each run's deconvolved mixture against its own moment-matched Gaussian (how non-Gaussian the landscape
is) and between the runs' moment-matched Gaussians (agreement of the first two moments alone). Run agreement only says more than
second moments if BC(run, run') exceeds these.  usage: land_xdgauss.py NPZ name1,name2,... [--nmc 200000]"""
import sys, numpy as np
from scipy.special import logsumexp
Z = np.load(sys.argv[1]); names = sys.argv[2].split(','); nmc = int(sys.argv[sys.argv.index('--nmc') + 1]) if '--nmc' in sys.argv else 200000; r = np.random.default_rng(3); L2P = np.log(2 * np.pi)
def mix(nm, part): return dict(pi=Z['%s_%s_pi' % (nm, part)], mu=Z['%s_%s_mu' % (nm, part)], Sig=Z['%s_%s_Sig' % (nm, part)])
def gauss(f):
    m = f['pi'] @ f['mu']; C = np.einsum('k,kij->ij', f['pi'], f['Sig'] + np.einsum('ki,kj->kij', f['mu'] - m, f['mu'] - m)); return dict(pi=np.ones(1), mu=m[None], Sig=C[None])
def lp(X, f):
    Si = np.linalg.inv(f['Sig']); ld = np.linalg.slogdet(f['Sig'])[1]
    return logsumexp(np.stack([np.log(f['pi'][k]) - 0.5 * (np.einsum('ni,ij,nj->n', X - f['mu'][k], Si[k], X - f['mu'][k]) + ld[k] + X.shape[1] * L2P) for k in range(len(f['pi']))], 1), 1)
def smp(f, m):
    c = r.choice(len(f['pi']), m, p=f['pi'] / f['pi'].sum()); Ls = np.linalg.cholesky(f['Sig']); return f['mu'][c] + np.einsum('nij,nj->ni', Ls[c], r.standard_normal((m, f['mu'].shape[1])))
def bc(f, g): X = smp(f, nmc); return float(np.mean(np.exp(0.5 * (lp(X, g) - lp(X, f)))))
for nm in names: f = mix(nm, 'all'); print('%-6s mixture vs its moment-matched Gaussian: BC %.4f (K %d)' % (nm, bc(f, gauss(f)), len(f['pi'])), flush=True)
for i, x in enumerate(names):
    for y in names[i + 1:]: print('%s vs %s: Gaussians only BC %.4f | mixtures BC %.4f' % (x, y, bc(gauss(mix(x, 'all')), gauss(mix(y, 'all'))), bc(mix(x, 'all'), mix(y, 'all'))), flush=True)
