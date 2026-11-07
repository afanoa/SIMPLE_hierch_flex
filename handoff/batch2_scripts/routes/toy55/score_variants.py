#!/usr/bin/env python3
"""usage: score_variants.py NAME=TABLE [NAME=TABLE ..] --truth T --tmasks A.mrc,B.mrc,C.mrc [--tcols C_closed,A_domain,B_blob --null 0.1]
Gate variants of the factorial site layer on a toy with truth. TABLE = the port's factorial/table.txt ('name level labels post rows [mask]';
labels 0 / j = pole j; post one column per pole, < 0 outside the level); masks from the sites.txt beside it. Per kept site: Dice with each
true site mask, and per truth column the signal max over poles of |AUC - 0.5| of the pole posterior on the level's genuine particles; the
site's variable = the column with the largest signal, NULL (a false site) when every signal < NULL. --cont 1: continuous truth
columns (arm angles): signal = max over poles of |Spearman(posterior, column)|, plus eta^2 of the column over the site's labels. Per variant: the closed site's best pole
(precision / recall on genuine particles, AUC, top-2.5k purity), the same for domain A, the false sites, and copies of one variable."""
import os, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('variants', nargs='+'); ap.add_argument('--truth', required=True); ap.add_argument('--tmasks', required=True)
ap.add_argument('--tcols', default='C_closed,A_domain,B_blob'); ap.add_argument('--null', type=float, default=0.1); ap.add_argument('--cont', type=int, default=0); a = ap.parse_args()
VAR = a.variants; TCOLS = a.tcols.split(','); NULL = a.null
hd = [l.split()[1:] for l in open(a.truth) if l.startswith('#')][0]; T = np.genfromtxt(a.truth, comments='#'); tr = T[:, 0].astype(int) - 1; NM = tr.max() + 1
def col(c): v = np.zeros(NM); v[tr] = T[:, hd.index(c)]; return v
junk = col('junk') > 0; TC = {c: (col(c) > 0) & ~junk for c in TCOLS} if not a.cont else {c: col(c) for c in TCOLS}
def spr(x, y): return float(np.corrcoef(np.argsort(np.argsort(x)), np.argsort(np.argsort(y)))[0, 1])
def eta2(lab, y): mu = y.mean(); return sum(((lab == k).sum() * (y[lab == k].mean() - mu) ** 2) for k in np.unique(lab)) / max(((y - mu) ** 2).sum(), 1e-30)
def rd(p): h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n) > 0.5
TM = [(os.path.basename(p)[5:-4], rd(p)) for p in a.tmasks.split(',')]
def auc(s, y): rk = np.argsort(np.argsort(s)) + 1.; n1 = y.sum(); n0 = len(y) - n1; return (rk[y].sum() - n1 * (n1 + 1) / 2) / max(n1 * n0, 1)
def rowmask(f):
    r = np.loadtxt(f, dtype=int).ravel()
    if len(r) == NM and r.max() <= 1: return r > 0
    m = np.zeros(NM, bool); m[r - 1] = True; return m
summary = []
for v in VAR:
    vn, tab = v.split('=', 1); d = os.path.dirname(tab); masks = {}
    if os.path.exists(d + '/sites.txt'):
        for l in open(d + '/sites.txt'):
            w = l.split('|')[0].split()
            if w and not l.startswith('#') and len(w) >= 10: masks[w[0]] = w[9]
    print('==== %s (%s)' % (vn, tab)); best = {}; null = []; byvar = {}; bestauc = {}
    for w in [l.split() for l in open(tab) if l.strip() and not l.startswith('#')]:
        name, lev, labf, postf, rowf = w[:5]
        if len(w) > 5: masks[name] = w[5]
        L = np.loadtxt(labf, dtype=int); Q = np.loadtxt(postf, ndmin=2); Q = Q.T if Q.shape[0] != NM else Q; rows = rowmask(rowf); g = rows & ~junk
        dice = ''
        if name in masks and os.path.exists(masks[name]):
            m = rd(masks[name]); dice = ' '.join('%s %.2f' % (k, 2. * (m & t).sum() / max(m.sum() + t.sum(), 1)) for k, t in TM)
        sig = {}; poles = sorted(set(np.unique(L[L > 0]).tolist()) | set(range(1, Q.shape[1] + 1)))
        for c, y in TC.items():
            s_ = []
            for j in range(Q.shape[1]):
                P = Q[:, j]; ok = g & np.isfinite(P) & (P >= 0)
                if a.cont: r_ = spr(P[ok], y[ok]); s_.append((abs(r_), abs(r_), j + 1))
                else: a_ = auc(P[ok], y[ok]); s_.append((abs(a_ - 0.5), max(a_, 1 - a_), j + 1))
            sig[c] = max(s_)
        if a.cont:   # continuous truth (arm angles): |Spearman| of the poles' posteriors and eta^2 of the site's labels
            e2 = {c: eta2(L[g], y[g]) for c, y in TC.items()}; var = max(sig, key=lambda c: sig[c][0]); isnull = sig[var][0] < NULL
            print('  %-8s level %-4s poles %s | %s -> %s' % (name, lev, ' '.join('%d:%d' % (j, (L == j).sum()) for j in poles if (L == j).sum()),
                  ' '.join('%s %.2f/%.2f' % (c, sig[c][0], e2[c]) for c in TC), 'NULL (false site)' if isnull else var))
            if isnull: null.append(name); continue
            byvar.setdefault(var, []).append(name)
            for c in TC:
                if c not in best or sig[c][0] > best[c][0]: best[c] = (sig[c][0], e2[c], name)
            continue
        var = max(sig, key=lambda c: sig[c][0]); isnull = sig[var][0] < NULL
        line = '  %-8s level %-4s poles %s | mask Dice: %s | signal |AUC-0.5|: %s -> %s' % (name, lev, ' '.join('%d:%d' % (j, (L == j).sum()) for j in poles if (L == j).sum()), dice or '-',
               ' '.join('%s %.3f' % (c, sig[c][0]) for c in TC), 'NULL (false site)' if isnull else var)
        print(line)
        if isnull: null.append(name); continue
        byvar.setdefault(var, []).append(name); y = TC[var]
        for j in poles:
            m = (L == j) & rows
            if not m.any(): continue
            prec = (m & y).sum() / max((m & ~junk).sum(), 1); rec = (m & y).sum() / max((g & y).sum(), 1); P = Q[:, min(j, Q.shape[1]) - 1]; ok = g & np.isfinite(P) & (P >= 0)
            a_ = auc(P[ok], y[ok]); o = np.argsort(-np.where(ok, P if a_ >= 0.5 else -P, -np.inf)); top = y[o[:2500]].mean()
            print('      pole %d: %d labelled, %s precision %.3f recall %.3f | its posterior AUC %.3f, top-2.5k %.3f' % (j, m.sum(), var, prec, rec, max(a_, 1 - a_), top))
            if var not in best or prec > best[var][0]: best[var] = (prec, rec, max(a_, 1 - a_), top, '%s pole %d' % (name, j), (g & y).sum() / max(g.sum(), 1))
            bestauc[var] = max(bestauc.get(var, (0., ''))[0], max(a_, 1 - a_)), ('%s pole %d' % (name, j)) if max(a_, 1 - a_) >= bestauc.get(var, (0., ''))[0] else bestauc[var][1]
    for c in TC:
        b = best.get(c)
        if a.cont: print('  %s: %s' % (c, 'best site %s: |Spearman| %.3f, eta^2 of its labels %.3f' % (b[2], b[0], b[1]) if b else 'NOT recovered')); continue
        print('  %s: %s' % (c, 'best pole %s: precision %.3f (base %.3f) recall %.3f, AUC %.3f, top-2.5k %.3f' % (b[4], b[0], b[5], b[1], b[2], b[3]) if b else 'NOT recovered'))
    print('  false sites (no truth signal >= %.2f): %d %s; copies: %s' % (NULL, len(null), ' '.join(null), ', '.join('%s x%d' % (k, len(v_)) for k, v_ in byvar.items() if len(v_) > 1) or 'none'))
    for c, (v_, w_) in bestauc.items(): print('  %s: best posterior AUC %.3f (%s)' % (c, v_, w_))
    summary.append((vn, best, null, bestauc))
print('==== verdicts')
ref = summary[0][1].get('A_domain') if summary else None
for vn, best, null, bestauc in summary:
    if a.cont: print('%-6s %s | false sites %d' % (vn, ' '.join('%s %.2f' % (c, best[c][0]) if c in best else '%s -' % c for c in TC), len(null))); continue
    c = best.get('C_closed'); a_ = best.get('A_domain')
    print('%-6s closed %s | A %s | false sites %d' % (vn, 'best pole precision %.3f recall %.3f (AUC %.3f); best posterior AUC %.3f' % (c[0], c[1], c[2], bestauc['C_closed'][0]) if c else 'NOT recovered',
          'AUC %.3f%s' % (a_[2], '' if ref is None or abs(a_[2] - ref[2]) <= 0.01 else ' (changed vs %s)' % summary[0][0]) if a_ else 'lost', len(null)))
