#!/usr/bin/env python3
"""usage: simple_exec prg=print_project_field projfile=P oritype=ptcl3D | quality_bins.py NAME NSYM [STATES_FILE] [REF_STATE]
The quality confound of a set of states or bins: per state (the project's ptcl3D state, or one value per row from STATES_FILE),
the particles, median corr and dfx, the top-view share (abs cos e2 > 0.8), and corr's AUC against REF_STATE (default: the
middle state of those present; 'max': the most populated) inside dfx deciles, so the AUC is at matched defocus (0.5 = no difference). For the two outermost
states: the view and dfx total-variation excess over 20 label shuffles (view_check.py's cells). One QUAL line per state, one QEND line."""
import sys, numpy as np
NAME, NSYM = sys.argv[1], int(sys.argv[2]); SF = sys.argv[3] if len(sys.argv) > 3 and sys.argv[3] != '-' else None
st, e1, e2, df, co = [], [], [], [], []
for l in sys.stdin:
    d = dict(t.split('=', 1) for t in l.split() if '=' in t)
    if 'e1' not in d: continue
    try: st.append(int(float(d.get('state', '0')))); e1.append(float(d['e1'])); e2.append(float(d['e2'])); df.append(float(d.get('dfx', 'nan'))); co.append(float(d.get('corr', 'nan')))
    except ValueError: continue
st, e1, e2, df, co = (np.array(x) for x in (st, e1, e2, df, co))
if SF: s2 = np.loadtxt(SF, dtype=int).ravel(); assert len(s2) == len(st), 'states file %d rows, project %d' % (len(s2), len(st)); st = s2
ks = [k for k in sorted(set(st.tolist())) if k > 0]
REF = (max(ks, key=lambda k: (st == k).sum()) if sys.argv[4] == 'max' else int(sys.argv[4])) if len(sys.argv) > 4 else ks[len(ks) // 2]
ok = np.isfinite(df); dd = np.full(len(df), -1); dd[ok] = np.searchsorted(np.quantile(df[ok], np.linspace(0.1, 0.9, 9)), df[ok])
def auc(a, b):
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 20 or len(b) < 20: return np.nan
    r = np.argsort(np.argsort(np.r_[a, b])); return (r[:len(a)].sum() - len(a) * (len(a) - 1) / 2.) / (len(a) * len(b))
cz = np.cos(np.radians(e2))
for k in ks:
    m = st == k; A = np.nanmean([auc(co[m & (dd == j)], co[(st == REF) & (dd == j)]) for j in range(10)]) if k != REF else 0.5
    print('QUAL %s state %d | particles %d | median corr %.4f | median dfx %.3f | top-view share %.3f | corr AUC vs state %d at matched defocus %.3f' % (NAME, k, m.sum(),
          np.nanmedian(co[m]), np.nanmedian(df[m]), np.mean(np.abs(cz[m]) > 0.8), REF, A))
lo, hi = ks[0], ks[-1]; vc = np.clip((cz + 1) / 2 * 10, 0, 9.999).astype(int) * 12 + np.clip(np.mod(e1, 360. / NSYM) / (360. / NSYM) * 12, 0, 11.999).astype(int)
def tv(cell, a, b, nc): p = np.bincount(cell[a], minlength=nc) / max(a.sum(), 1); q = np.bincount(cell[b], minlength=nc) / max(b.sum(), 1); return 0.5 * np.abs(p - q).sum()
rng = np.random.default_rng(5); e = (st == lo) | (st == hi); out = []
for cell, nc in ((vc, 120), (np.maximum(dd, 0), 10)):
    t0 = tv(cell, st == lo, st == hi, nc); lab = st[e].copy(); c = cell[e]; tp = []
    for _ in range(20): rng.shuffle(lab); tp.append(tv(c, lab == lo, lab == hi, nc))
    out.append(t0 - np.mean(tp))
print('QEND %s states %d vs %d | view TV excess %.3f | dfx TV excess %.3f | corr AUC %d vs %d at matched defocus %.3f' % (NAME, lo, hi, out[0], out[1], hi, lo,
      np.nanmean([auc(co[(st == hi) & (dd == j)], co[(st == lo) & (dd == j)]) for j in range(10)])))
