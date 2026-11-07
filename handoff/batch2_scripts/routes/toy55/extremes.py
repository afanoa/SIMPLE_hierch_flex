#!/usr/bin/env python3
"""usage: extremes.py TRUTH OUT NAME:TREE_DIR:TABLE[:FLEX_EMB:FFLEX_EMB] ..
Do candidate particle sets reach the extremes of each truth motion coordinate (the IgM arm angles) at a usable size?
- Genuine particles only. For coordinate t and end (low / high), e(S) = (mean_t(S) - mean_t(G)) / sd_t(G), signed toward that end,
  over S's genuine members.
- e*(m) = the mean standardised t of the m most extreme genuine particles at that end; ratio = e(S) / e*(m), where m = S's genuine count.
- Per end and n_min (2k, 5k, 10k): the best ratio among a family's sets with m >= n_min, with that set's m, junk share and name.
Families per tree:
- leaves; merged units (hetero_tree_merged.txt);
- site poles (the layer's labels; and posterior > 0.5);
- site score tails (top / bottom n of each pole's posterior);
- root flex and fflex latent tails (top / bottom n of each raw MAP z component).
Tails use n = 2k, 5k, 10k, 20k. Writes OUT.tsv (every end x n_min x family) and prints the summary and the per-end bests."""
import sys, numpy as np
TF, OUT = sys.argv[1:3]; SPECS = [a.split(':') for a in sys.argv[3:]]
hd = [l.split()[1:] for l in open(TF) if l.startswith('#')][0]; T = np.genfromtxt(TF, comments='#'); tr = T[:, 0].astype(int) - 1; NM = tr.max() + 1
ang = hd[3:]; X = np.full((NM, len(ang)), np.nan); X[tr] = T[:, 3:3 + len(ang)]
junk = np.zeros(NM, bool); junk[tr] = T[:, hd.index('junk')] > 0; G = ~junk & ~np.isnan(X[:, 0])
Z = (X - X[G].mean(0)) / X[G].std(0); NMIN, NT = (2000, 5000, 10000), (2000, 5000, 10000, 20000)
ENDS = [(c, s) for c in range(len(ang)) for s in (-1, 1)]; EN = ['%s %s' % (ang[c], 'low' if s < 0 else 'high') for c, s in ENDS]
EST = [np.cumsum(np.sort(s * Z[G, c])[::-1]) / np.arange(1, G.sum() + 1) for c, s in ENDS]   # e*(m) = EST[e][m - 1]
def read_emb(f):
    b = open(f, 'rb').read(); ver, n, nc, box = (int(v) for v in np.frombuffer(b[8:24], '<i4'))
    assert b[:8] == b'SIMPLFXC' and ver == 5, 'flex_pca embedding version 5 expected'
    o = 28; pinds = np.frombuffer(b[o:o + 4 * n], '<i4').astype(int); o += 4 * n
    return pinds - 1, np.frombuffer(b[o:o + 8 * n * nc], '<f8').reshape(nc, n).T.copy()
def table(D, tab):
    p = tab if tab.startswith('/') else D + '/' + tab; out = []
    for w in [l.split() for l in open(p) if l.strip() and not l.startswith('#')]:
        Q = np.loadtxt(w[3], ndmin=2); Q = Q if Q.shape[0] == NM else Q.T; r = np.loadtxt(w[4], dtype=int).ravel()
        rows = (r > 0) if (len(r) == NM and r.max() <= 1) else np.isin(np.arange(NM), r - 1); out.append((w[0], np.loadtxt(w[2], dtype=int), Q, rows))
    return out
def sets(sp):
    nm, D, tab = sp[:3]; L = np.loadtxt(D + '/hetero_tree_labels.txt', dtype=int)[:NM]
    names = [l.split()[0] for l in open(D + '/hetero_tree_leaves.txt') if l.strip()]
    for k in np.unique(L[L > 0]): yield 'leaves', names[k - 1], np.where(L == k)[0]
    try: U = np.loadtxt(D + '/hetero_tree_merged.txt', dtype=int)[:NM]
    except OSError: U = None
    if U is not None:
        for u in np.unique(U[U > 0]): yield 'units', 'u%d' % u, np.where(U == u)[0]
    for s, lab, Q, rows in table(D, tab):
        for j in range(Q.shape[1]):
            yield 'site poles (labels)', '%s pole %d' % (s, j + 1), np.where((lab == j + 1) & rows)[0]
            yield 'site poles (post > 0.5)', '%s pole %d' % (s, j + 1), np.where(Q[:, j] > 0.5)[0]
            v = np.where(Q[:, j] >= 0)[0]; o = v[np.argsort(-Q[v, j], kind='stable')]
            for n in NT:
                if n < len(o): yield 'site score tails', '%s pole %d top %d' % (s, j + 1, n), o[:n]; yield 'site score tails', '%s pole %d bottom %d' % (s, j + 1, n), o[-n:]
    for fam, f in (('flex latent tails', sp[3] if len(sp) > 3 else ''), ('fflex latent tails', sp[4] if len(sp) > 4 else '')):
        if not f: continue
        rows, z = read_emb(f)
        for c in range(z.shape[1]):
            o = rows[np.argsort(z[:, c], kind='stable')]
            for n in NT: yield fam, 'z%d top %d' % (c + 1, n), o[-n:]; yield fam, 'z%d bottom %d' % (c + 1, n), o[:n]
best = {}   # (tree, family, end, nmin) -> (ratio, m, junk share, name)
for sp in SPECS:
    for fam, name, r in sets(sp):
        g = r[G[r]]; m = len(g)
        if m < NMIN[0]: continue
        mu = Z[g].mean(0); js = 1. - m / len(r)
        for e, (c, s) in enumerate(ENDS):
            ratio = s * mu[c] / EST[e][m - 1]
            for nmin in NMIN:
                k = (sp[0], fam, e, nmin)
                if m >= nmin and (k not in best or ratio > best[k][0]): best[k] = (ratio, m, js, name)
fams = ['leaves', 'units', 'site poles (labels)', 'site poles (post > 0.5)', 'site score tails', 'flex latent tails', 'fflex latent tails']
with open(OUT + '.tsv', 'w') as f:
    f.write('tree\tfamily\tend\tn_min\tratio\tgenuine\tjunk_share\tset\n')
    for (t, fam, e, nmin), (ratio, m, js, name) in sorted(best.items(), key=lambda x: (x[0][0], fams.index(x[0][1]), x[0][2], x[0][3])):
        f.write('%s\t%s\t%s\t%d\t%.3f\t%d\t%.3f\t%s\n' % (t, fam, EN[e], nmin, ratio, m, js, name))
print('oracle e*(n), median over the %d ends: %s' % (len(ENDS), ', '.join('n %d: %.2f' % (n, np.median([E[n - 1] for E in EST])) for n in NT)))
print('| tree | family | ' + ' | '.join('n_min %d: median ratio, ends >= 0.5 / >= 0.7' % n for n in NMIN) + ' |'); print('|---|---|' + '---|' * len(NMIN))
for sp in SPECS:
    for fam in fams:
        cells = []
        for nmin in NMIN:
            rs = [best[(sp[0], fam, e, nmin)][0] if (sp[0], fam, e, nmin) in best else 0. for e in range(len(ENDS))]
            cells.append('%.2f, %d / %d' % (np.median(rs), sum(r >= 0.5 for r in rs), sum(r >= 0.7 for r in rs)) if any((sp[0], fam, e, nmin) in best for e in range(len(ENDS))) else '-')
        print('| %s | %s | %s |' % (sp[0], fam, ' | '.join(cells)))
for sp in SPECS:
    print('== %s: best set per end at n_min 5000 (ratio, genuine, junk share, set) per family' % sp[0])
    for e in range(len(ENDS)):
        print('   %-12s e*(5k) %.2f | %s' % (EN[e], EST[e][4999], ' | '.join('%s %.2f %d %.2f %s' % (fam.split()[0] + ('+' if 'post' in fam else ''), *best[(sp[0], fam, e, 5000)]) for fam in fams if (sp[0], fam, e, 5000) in best)))
