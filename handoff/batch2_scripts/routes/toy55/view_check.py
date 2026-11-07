#!/usr/bin/env python3
"""usage: simple_exec prg=print_project_field projfile=P oritype=ptcl3D | view_check.py NAME NSYM
Do a 5-bin series' bins (the project's ptcl3D states 1..5) differ in view or defocus? Views are histogrammed on 10 cos(e2) bands x
12 e1 sectors over one asymmetric unit of cNSYM (equal-area cells), defocus (dfx) on its 10 quantiles over the binned rows. For the
end bins 1 and 5: the total-variation distance TV = 0.5 sum abs(p1 - p5) and its mean over 20 shuffles of the two bins' labels (the
finite-count baseline). Also per bin: the share of top views (abs cos e2 > 0.8) and the mean dfx. One VIEW line."""
import sys, numpy as np
NAME, NSYM = sys.argv[1], int(sys.argv[2]); st, e1, e2, df = [], [], [], []
for l in sys.stdin:
    d = dict(t.split('=', 1) for t in l.split() if '=' in t)
    if 'state' not in d or 'e1' not in d: continue
    try: s = int(float(d['state']))
    except ValueError: continue
    if 1 <= s <= 5: st.append(s); e1.append(float(d['e1'])); e2.append(float(d['e2'])); df.append(float(d.get('dfx', 'nan')))
st, e1, e2, df = np.array(st), np.array(e1), np.array(e2), np.array(df)
if len(st) == 0: print('VIEW %s | no binned rows' % NAME); sys.exit()
cz = np.cos(np.radians(e2)); ph = np.mod(e1, 360. / NSYM) / (360. / NSYM); vc = np.clip((cz + 1) / 2 * 10, 0, 9.999).astype(int) * 12 + np.clip(ph * 12, 0, 11.999).astype(int)
okd = np.isfinite(df); dq = np.searchsorted(np.quantile(df[okd], np.linspace(0.1, 0.9, 9)), df) if okd.any() else np.zeros(len(df), int)
def tv(cell, a, b, nc): p = np.bincount(cell[a], minlength=nc) / max(a.sum(), 1); q = np.bincount(cell[b], minlength=nc) / max(b.sum(), 1); return 0.5 * np.abs(p - q).sum()
rng = np.random.default_rng(5); e = (st == 1) | (st == 5); out = []
for cell, nc in ((vc, 120), (dq, 10)):
    t0 = tv(cell, st == 1, st == 5, nc); lab = st[e].copy(); c = cell[e]; tp = []
    for _ in range(20): rng.shuffle(lab); tp.append(tv(c, lab == 1, lab == 5, nc))
    out.append((t0, np.mean(tp)))
top = [np.mean(np.abs(cz[st == k]) > 0.8) for k in range(1, 6)]; mdf = [np.nanmean(df[st == k]) for k in range(1, 6)]
print('VIEW %s | rows/bin %s | views TV %.3f (shuffled %.3f) | dfx TV %.3f (shuffled %.3f) | top-view share %s | mean dfx %s' % (NAME, ' '.join(str((st == k).sum()) for k in range(1, 6)),
      out[0][0], out[0][1], out[1][0], out[1][1], ' '.join('%.3f' % x for x in top), ' '.join('%.3f' % x for x in mdf)))
