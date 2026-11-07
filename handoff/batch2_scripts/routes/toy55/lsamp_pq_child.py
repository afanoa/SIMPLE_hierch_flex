#!/usr/bin/env python3
"""lsamp_pq_child.py: does the closed child's own pole pair (site_rank at depth 1) carry the root's second landscape axis? On the child's
particles: |Spearman| of the child's pole posteriors (z1lo / z1hi) with root coordinates (NAME=coords.npz: axis-2 landscapes, pole-pair
coordinates), the child pole pair across runs, and on the root particles outside the child the same for the open-open pole pair.
usage: lsamp_pq_child.py --child RUN=POSTDIR,.. --root NAME=COORDS,.. """
import argparse, numpy as np
from scipy.stats import spearmanr
ap = argparse.ArgumentParser(); ap.add_argument('--child', required=True); ap.add_argument('--root', required=True); a = ap.parse_args()
CH = {}
for e in a.child.split(','):
    nm, d = e.split('='); CH[nm] = {t: np.loadtxt('%s/site_post_%s.txt' % (d, t)) for t in ('z1lo', 'z1hi')}
NM = len(next(iter(CH.values()))['z1lo']); RT = {}
for e in a.root.split(','):
    nm, f = e.split('='); C = np.load(f); t = np.full(NM, np.nan); t[C['rows']] = C['tfull'][:, 0]; RT[nm] = t
def rho(x, y):
    m = np.isfinite(x) & np.isfinite(y); return spearmanr(x[m], y[m])[0] if m.sum() > 100 else np.nan
for nm, P in CH.items():
    inc = P['z1lo'] >= 0; print('== child %s: %d particles; z1lo vs z1hi posterior Spearman %.3f' % (nm, inc.sum(), rho(np.where(inc, P['z1lo'], np.nan), np.where(inc, P['z1hi'], np.nan))))
    for rn, t in RT.items(): print('   %-14s |rho| with child z1lo %.3f, z1hi %.3f' % (rn, abs(rho(np.where(inc, P['z1lo'], np.nan), t)), abs(rho(np.where(inc, P['z1hi'], np.nan), t))))
nms = list(CH)
for i in range(len(nms)):
    for j in range(i + 1, len(nms)):
        x, y = CH[nms[i]], CH[nms[j]]; both = (x['z1lo'] >= 0) & (y['z1lo'] >= 0)
        print('child pole pair across runs (%s vs %s, %d shared particles): Spearman z1lo/z1lo %+.3f, z1lo/z1hi %+.3f' % (nms[i], nms[j], both.sum(), rho(np.where(both, x['z1lo'], np.nan), y['z1lo']), rho(np.where(both, x['z1lo'], np.nan), y['z1hi'])))
out = np.ones(NM, bool)
for P in CH.values(): out &= P['z1lo'] < 0
rn = [r for r in RT if r.startswith('oA_oB')]; ax = [r for r in RT if r.startswith('ax2')]
print('== root particles outside every closed child (%d with root coordinates): |rho| of the open-open pole pair with axis 2' % int((out & np.isfinite(RT[rn[0]])).sum()) if rn else '')
for r in rn: print('   %-14s %s' % (r, ' '.join('%s %.3f' % (x, abs(rho(np.where(out, RT[r], np.nan), RT[x]))) for x in ax)))
