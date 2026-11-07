#!/usr/bin/env python3
"""Score a hetero_tree on a multi-motion toy (make_mm_maps.py + qtoy_merge.py truth): how the tree represents several INDEPENDENT motions.
usage: mm_score.py TREE_DIR [truth=SIM/truth_qtoy.txt] [factors=MAPS/factors.txt] [tag=NAME]
TREE_DIR holds hetero_tree_labels.txt / _leaves.txt / _merged.txt (one label per project row, 0 = no leaf).
Reports, on genuine (non-junk) particles:
 1 leaves: size, junk share, the majority state and purity of every factor, top combinations;
 2 combinations: best-leaf F1 per combination, how many leaves hold >= 10 % of it (fragmentation), combinations that are a leaf (F1 >= 0.7);
 3 splits: per internal node, Cramer's V of the child label with every factor; the factor each split follows (V >= 0.3), and REPEATS
   (the same factor split again in disjoint branches = a motion rediscovered below another motion's split);
 4 factor states (the 'extremes' of a two-state motion): best single leaf, best tree node and best merged unit F1 for each state, and
   the leaf-majority accuracy (each leaf votes its majority state: does the leaf set still carry the motion, however fragmented);
 5 merge-up: units, and for every unit of >= 2 leaves whether its leaves share a combination (dedup) or differ in a factor (mixed).
A continuous factor (factors.txt cont_angle varying) is scored by eta^2 of the angle on leaves and its 10 % tails' best-leaf F1.
One SUMMARY line at the end."""
import sys, os, numpy as np
args = dict(a.split('=', 1) if '=' in a else ('tree', a) for a in sys.argv[1:]); T = args['tree']; tag = args.get('tag', os.path.basename(T.rstrip('/')))
truth = args.get('truth', '/mnt/beegfs/elmlund/afan/routes/toy55/simmm3/truth_qtoy.txt'); facf = args.get('factors', '/mnt/beegfs/elmlund/afan/routes/toy55/maps160/mm3/factors.txt')
COLS = 'row combo junk A_domain B_blob C_closed jtype conformer s1 closed_arm holo bfac snrq'.split()
t = np.loadtxt(truth, ndmin=2); rows = t[:, 0].astype(int); combo = t[:, 1].astype(int); junk = t[:, 2] > 0
F = np.loadtxt(facf, ndmin=2); nf = F.shape[1] - 3; fmap = {int(r[0]): r[1:1 + nf].astype(int) for r in F}; amap = {int(r[0]): r[1 + nf] for r in F}
gen = ~junk; fs = np.zeros((len(rows), nf), int); ang = np.full(len(rows), np.nan)
for i in np.where(gen)[0]: fs[i] = fmap[combo[i]]; ang[i] = amap[combo[i]]
cont = np.nanstd(ang) > 1e-6
leafnames = [l.split()[0] for l in open(T + '/hetero_tree_leaves.txt') if l.strip()]
lab = np.loadtxt(T + '/hetero_tree_labels.txt', dtype=int, ndmin=1)[rows - 1]
unit = np.loadtxt(T + '/hetero_tree_merged.txt', dtype=int, ndmin=1)[rows - 1] if os.path.exists(T + '/hetero_tree_merged.txt') else lab.copy()
def nmi(a, b):
    a = np.unique(a, return_inverse=True)[1]; b = np.unique(b, return_inverse=True)[1]; n = len(a)
    p = np.zeros((a.max() + 1, b.max() + 1)); np.add.at(p, (a, b), 1); p /= n; pa, pb = p.sum(1), p.sum(0); nz = p > 0
    mi = (p[nz] * np.log(p[nz] / (pa[:, None] * pb[None, :])[nz])).sum(); ha = -(pa[pa > 0] * np.log(pa[pa > 0])).sum(); hb = -(pb[pb > 0] * np.log(pb[pb > 0])).sum()
    return mi / max(np.sqrt(ha * hb), 1e-12)
def cramer(a, b):
    a = np.unique(a, return_inverse=True)[1]; b = np.unique(b, return_inverse=True)[1]
    if a.max() == 0 or b.max() == 0: return 0.
    o = np.zeros((a.max() + 1, b.max() + 1)); np.add.at(o, (a, b), 1); e = o.sum(1)[:, None] * o.sum(0)[None] / o.sum()
    return float(np.sqrt(((o - e) ** 2 / np.maximum(e, 1e-12)).sum() / o.sum() / (min(o.shape) - 1)))
def f1(sel, tgt):
    tp = (sel & tgt).sum(); return 2. * tp / max(sel.sum() + tgt.sum(), 1)
combos = sorted(fmap); ul = [k for k in np.unique(lab) if k > 0]; uu = [k for k in np.unique(unit) if k > 0]
gl = gen & (lab > 0)
print('== %s: %d leaves, %d units over %d rows (%d genuine, %d junk); %d factors%s, %d combinations' % (tag, len(ul), len(uu), (lab > 0).sum(), gen.sum(), junk.sum(), nf, ' + continuous angle' if cont else '', len(combos)))
# 1 leaves
print('-- 1 leaves (genuine purity per factor: majority state:purity)')
print('%-14s %6s %5s  %s  %s' % ('leaf', 'rows', 'junk', '  '.join('f%d' % (k + 1) for k in range(nf)), 'top combos'))
lpur = {}
for k in sorted(ul, key=lambda k: -(lab == k).sum()):
    m = lab == k; mg = m & gen; s = []; lpur[k] = []
    for j in range(nf):
        if mg.sum() == 0: s.append('  -   '); lpur[k].append((0, 0.)); continue
        v, c = np.unique(fs[mg, j], return_counts=True); i = np.argmax(c); s.append('%d:%.2f' % (v[i], c[i] / mg.sum())); lpur[k].append((v[i], c[i] / mg.sum()))
    cb, cn = np.unique(combo[m], return_counts=True); o = np.argsort(-cn)[:3]
    print('%-14s %6d %5.2f  %s  %s' % (leafnames[k - 1], m.sum(), junk[m].mean(), ' '.join(s), ' '.join('%d(%.2f)' % (cb[i], cn[i] / m.sum()) for i in o)))
# 2 combinations
print('-- 2 combinations: n, best leaf F1 (leaf), leaves holding >= 10 % of the combination, best unit F1')
cleaf = 0; frag = []
for c in combos:
    tg = gen & (combo == c); n = tg.sum()
    if n == 0: continue
    sc = [(f1(lab == k, tg), k) for k in ul]; b, kb = max(sc); su = max(f1(unit == u, tg) for u in uu)
    nfr = sum(1 for k in ul if (tg & (lab == k)).sum() >= 0.1 * n); frag.append(nfr); cleaf += b >= 0.7
    print('%6d %6d  leaf F1 %.2f (%s)  spread over %d leaves  unit F1 %.2f' % (c, n, b, leafnames[kb - 1], nfr, su))
# 3 splits: node rows from the leaf-name paths (rX contains every leaf whose name starts with rX)
def node_of(name): return name.split('_')[0]
nodes = sorted({node_of(n)[:d] for n in leafnames for d in range(1, len(node_of(n)) + 1)}, key=lambda x: (len(x), x))
def child_key(name, node):   # child of node that leaf `name` falls in: the next path digit, or the stage state suffix
    p = node_of(name)
    if len(p) > len(node): return p[len(node)]
    return name.split('_', 1)[1] if p == node else None
print('-- 3 splits (genuine rows): node, rows, children, Cramer V per factor, junk V; follows = factors with V >= 0.3')
follow = {j: [] for j in range(nf)}; nsplit = 0
for nd in nodes:
    ks = [k for k in ul if node_of(leafnames[k - 1]).startswith(nd)]
    ch = {}
    for k in ks:
        ck = child_key(leafnames[k - 1], nd)
        if ck is not None: ch.setdefault(ck, []).append(k)
    if len(ch) < 2: continue
    nsplit += 1; m = np.isin(lab, ks); cl = np.zeros(len(lab), int)
    for i, (ck, kk) in enumerate(sorted(ch.items())): cl[np.isin(lab, kk)] = i + 1
    mg = m & gen; V = [cramer(cl[mg], fs[mg, j]) for j in range(nf)]; vj = cramer(cl[m], junk[m].astype(int))
    fol = [j for j in range(nf) if V[j] >= 0.3]; [follow[j].append(nd) for j in fol]
    ea = ''
    if cont:
        a = ang[mg]; mu = np.array([a[cl[mg] == i].mean() for i in np.unique(cl[mg])]); ea = ' angle eta2 %.2f' % (sum((cl[mg] == i).sum() * (mu[q] - a.mean()) ** 2 for q, i in enumerate(np.unique(cl[mg]))) / max(((a - a.mean()) ** 2).sum(), 1e-12))
    print('%-8s %6d  K %d  %s  junk %.2f%s  follows %s' % (nd, m.sum(), len(ch), ' '.join('f%d %.2f' % (j + 1, V[j]) for j in range(nf)), vj, ea, ','.join('f%d' % (j + 1) for j in fol) or '-'))
def disjoint(ns): return [n for n in ns if not any(o != n and o.startswith(n) for o in ns)]   # deepest nodes of the list = its largest antichain (branches)
rep = {j: len(disjoint(follow[j])) for j in range(nf)}
print('-- repeats: ' + '; '.join('f%d split at %s (%d disjoint branch%s)' % (j + 1, ','.join(follow[j]) or 'none', rep[j], '' if rep[j] == 1 else 'es') for j in range(nf)))
# 4 factor states
print('-- 4 factor states: n, best leaf F1 (leaf), best node F1 (node), best unit F1, leaf-majority accuracy / NMI of the factor')
nodemask = {nd: np.isin(lab, [k for k in ul if node_of(leafnames[k - 1]).startswith(nd)]) for nd in nodes}
summ = []
for j in range(nf):
    acc = 0
    for k in ul:
        mg = gen & (lab == k)
        if mg.sum(): acc += np.bincount(fs[mg, j]).max()
    acc /= max(gl.sum(), 1); nm = nmi(lab[gl], fs[gl, j])
    for s_ in (1, 2):
        tg = gen & (fs[:, j] == s_); bl = max((f1(lab == k, tg), k) for k in ul); bn = max((f1(nodemask[nd], tg), nd) for nd in nodes if nd != 'r'); bu = max(f1(unit == u, tg) for u in uu)
        print('f%d=%d %6d  leaf %.2f (%s)  node %.2f (%s)  unit %.2f   | factor f%d leaf-majority acc %.3f, NMI %.3f' % (j + 1, s_, tg.sum(), bl[0], leafnames[bl[1] - 1], bn[0], bn[1], bu, j + 1, acc, nm))
        summ.append('f%d=%d %.2f/%.2f/%.2f' % (j + 1, s_, bl[0], bn[0], bu))
    summ.append('f%d acc %.2f' % (j + 1, acc))
if cont:
    a = ang[gl]; mu = {k: ang[gen & (lab == k)].mean() for k in ul if (gen & (lab == k)).any()}
    eta = sum((gen & (lab == k)).sum() * (mu[k] - a.mean()) ** 2 for k in mu) / max(((a - a.mean()) ** 2).sum(), 1e-12)
    lo, hi = np.nanpercentile(ang[gen], [10, 90]); tl = gen & (ang <= lo); th = gen & (ang >= hi)
    print('continuous angle: leaf eta2 %.3f; 10 %% tails best leaf F1 low %.2f high %.2f' % (eta, max(f1(lab == k, tl) for k in ul), max(f1(lab == k, th) for k in ul)))
    summ.append('angle eta2 %.2f' % eta)
# 5 merge-up
print('-- 5 merge-up: %d leaves -> %d units' % (len(ul), len(uu)))
dedup = mixed = 0
for u in uu:
    ks = sorted(set(lab[unit == u]) - {0})
    if len(ks) < 2: continue
    maj = []
    for k in ks:
        mg = gen & (lab == k); maj.append(np.bincount(combo[mg] if mg.any() else np.array([0])).argmax())
    diff = [j for j in range(nf) if len({lpur[k][j][0] for k in ks}) > 1]
    same = len(set(maj)) == 1; dedup += same; mixed += not same
    print('unit %d: %s  majority combos %s  -> %s' % (u, ' + '.join(leafnames[k - 1] for k in ks), ' '.join(str(m) for m in maj), 'same combination' if same else 'differ in ' + ','.join('f%d' % (j + 1) for j in diff)))
# 6 site layer (factorial/table.txt: one kept site per line, its level's rows and per-row pole labels / posteriors): every pole's F1
# against every factor state over the level's genuine rows, the best match, and the top-n posterior tail (n = that state's size)
sl = []
ft = T + '/factorial/table.txt'
def loc(pth): return pth if os.path.exists(pth) else os.path.join(T, pth[pth.index('1_hetero_tree') + len('1_hetero_tree') + 1:]) if '1_hetero_tree' in pth else pth
if os.path.exists(ft):
    print('-- 6 site layer: site level poles; per pole the best factor state by F1 over the level (and its top-posterior tail F1)')
    for line in open(ft):
        if line.startswith('#') or not line.strip(): continue
        nm, lev, lf, pf, rf = line.split()[:5]
        try: lr = np.loadtxt(loc(rf), dtype=int, ndmin=1); ll = np.loadtxt(loc(lf), dtype=int, ndmin=1); pp = np.loadtxt(loc(pf), ndmin=2)
        except OSError: print('%s: files missing' % nm); continue
        if pp.shape[0] != len(ll) and pp.shape[1] == len(ll): pp = pp.T
        if len(ll) != len(lr): ll, pp = ll[lr - 1], pp[lr - 1]   # labels / posteriors per project row (-1 outside the level)
        pos = {r_: i for i, r_ in enumerate(rows)}; idx = np.array([pos.get(r_, -1) for r_ in lr]); ok = idx >= 0; idx, ll, pp = idx[ok], ll[ok], pp[ok]
        g = gen[idx]; out = []
        for j in range(1, ll.max() + 1):
            sel = ll == j; best = (0, '')
            for q in range(nf):
                for s_ in (1, 2):
                    tg = g & (fs[idx, q] == s_); f = 2. * (sel & tg).sum() / max(sel.sum() + tg.sum(), 1)
                    if f > best[0]:
                        o = np.argsort(-pp[:, j - 1])[:tg.sum()]; tail = 2. * tg[o].sum() / max(len(o) + tg.sum(), 1); best = (f, 'f%d=%d' % (q + 1, s_), tail)
            out.append('pole %d (%d rows, junk %.2f): %s F1 %.2f tail %.2f' % (j, sel.sum(), junk[idx][sel].mean() if sel.any() else 0, best[1], best[0], best[2] if len(best) > 2 else 0))
            sl.append((nm, best[1], best[0]))
        print('%-8s level %-5s %6d rows  %s' % (nm, lev, len(idx), '; '.join(out)))
sz = np.array([(gen & (lab == k)).sum() for k in ul]); gsz = sz[sz >= 0.5 * np.array([(lab == k).sum() for k in ul])]
print('SUMMARY %s leaves %d (genuine-majority %d, median %d, min %d) units %d | NMI combo leaves %.3f units %.3f | combination leaves (F1>=0.7) %d/%d, median spread %.1f | splits %d, repeats %s | %s | merge-up dedup %d mixed %d' %
      (tag, len(ul), len(gsz), np.median(gsz) if len(gsz) else 0, gsz.min() if len(gsz) else 0, len(uu), nmi(lab[gl], combo[gl]), nmi(unit[gen & (unit > 0)], combo[gen & (unit > 0)]), cleaf, len(combos), np.median(frag), nsplit,
       ' '.join('f%d:%d' % (j + 1, rep[j]) for j in range(nf)), ' '.join(summ), dedup, mixed) + (' | sites ' + ' '.join('%s:%s %.2f' % x for x in sl) if sl else ''))
