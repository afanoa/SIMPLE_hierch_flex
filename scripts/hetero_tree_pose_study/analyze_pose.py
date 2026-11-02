#!/usr/bin/env python3
"""analyze_pose.py PAIRS.tsv NODES.tsv SAVED.tsv [RNASEP_PAIRS.tsv]: joins pose_split.py's per-pair statistics to split_effect.py's
effect-size table (by tree, node and child sizes), labels every pair (needed / unneeded by analyze_splits.py's rule; micelle fixture
trees = pose artefact; real PfCRT = pose-suspect; OLE unknown), and reports the distributions, the AUCs and what a pose gate would do.
Pure Python, no numpy."""
import sys, math
from collections import defaultdict


def read_tsv(p):
    rows = []
    with open(p) as f:
        hdr = f.readline().rstrip('\n').split('\t')
        for ln in f:
            v = ln.rstrip('\n').split('\t')
            if len(v) == len(hdr):
                rows.append(dict(zip(hdr, v)))
    return rows


def fl(x):
    try:
        return float(x)
    except Exception:
        return float('nan')


def auc(pos, neg):
    pos, neg = [p for p in pos if not math.isnan(p)], [n for n in neg if not math.isnan(n)]
    if not pos or not neg:
        return float('nan')
    s = 0.
    for p in pos:
        for n in neg:
            s += 1. if p > n else 0.5 if p == n else 0.
    return s / (len(pos) * len(neg))


def label(r):
    t, nd = r['tree'], r['node']
    if t.startswith('fx_mic') or 'micnull' in t:
        return 'artefact', 'micelle fixture'
    if t.startswith('pfcrt_'):
        return 'suspect', 'real PfCRT'
    if 'null' in t:
        return 'unneeded', 'null'
    V, z = fl(r.get('cramerV', 'nan')), fl(r.get('permz', 'nan'))
    if not math.isnan(V):
        return ('needed' if (V >= 0.3 and z >= 5) else 'unneeded'), 'truth V %.2f z %.0f' % (V, z)
    if 'r131' in t:
        return ('needed', 'swing / closed corner') if nd in ('r', 'r1') else ('unknown', 'refused r2')
    if 'ole' in t:
        return 'ole', 'no truth'
    if t == 'rnasep692':
        return 'production', '692k'
    return 'unknown', '-'


pairs, nodes, saved = read_tsv(sys.argv[1]), read_tsv(sys.argv[2]), read_tsv(sys.argv[3])
if len(sys.argv) > 4:
    pairs += read_tsv(sys.argv[4])
# join: effect-size rows by (tree, node, {n_a, n_b})
eff = defaultdict(list)
for r in saved:
    eff[(r['tree'], r['node'])].append(r)
for p in pairs:
    key = (p['tree'], p['node']); sz = {int(p['n_a']), int(p['n_b'])}
    hit = [r for r in eff.get(key, []) if {int(r['n_a']), int(r['n_b'])} == sz]
    if not hit:  # the stage-B states refined carry the same rows as the split groups up to junk rows: match within 2 %
        hit = [r for r in eff.get(key, []) if all(any(abs(int(a) - b) <= 0.02 * b + 5 for b in sz) for a in (r['n_a'], r['n_b']))]
    for k in ('decision', 'P_over_S', 'P_over_N', 'L', 'ratio', 'cramerV', 'permz'):
        p[k] = hit[0][k] if hit else 'nan'
    p['label'], p['why'] = label(p)

print('# %d pairs; labels: %s' % (len(pairs), ', '.join('%s %d' % (k, sum(p['label'] == k for p in pairs)) for k in ('needed', 'unneeded', 'artefact', 'suspect', 'ole', 'production', 'unknown'))))
print('\n## every pair (sorted by z_dflip)')
print('tree\tnode\tkind\tpair\tn_a\tn_b\tflip_a\tflip_b\tdflip\tz_dflip\tks_ang\tz_ks\tz_ks_shift\tview_tv\tP/N\tL/Lref\tV\tlabel\twhy')
for p in sorted(pairs, key=lambda p: -fl(p['z_dflip'])):
    print('\t'.join([p['tree'], p['node'], p['kind'], p['pair'], p['n_a'], p['n_b'], p['flip_a'], p['flip_b'], p['dflip'], p['z_dflip'],
                     p['ks_ang'], p['z_ks_ang'], p['z_ks_shift'], p['view_tv'], p['P_over_N'], p['ratio'], p['cramerV'], p['label'], p['why']]))

print('\n## node-level flip fractions (stage B, > 30 deg) by tree')
bytree = defaultdict(list)
for n in nodes:
    bytree[n['tree']].append(fl(n['flip30_stB']))
for t, v in sorted(bytree.items()):
    v = [x for x in v if not math.isnan(x)]
    if v:
        print('   %-22s nodes %3d  flip30 max %.4f median %.4f' % (t, len(v), max(v), sorted(v)[len(v) // 2]))

stats = ('z_dflip', 'dflip', 'z_ks_ang', 'ks_ang', 'z_ks_shift', 'view_tv')
print('\n## AUC: pose-artefact pairs (micelle fixtures) against needed splits, and against every legitimate-tree pair (needed + unneeded + ole)')
art = [p for p in pairs if p['label'] == 'artefact']
for name, neg in (('needed', [p for p in pairs if p['label'] == 'needed']),
                  ('all legitimate (needed, unneeded, ole)', [p for p in pairs if p['label'] in ('needed', 'unneeded', 'ole')])):
    print('   against %s (%d vs %d):' % (name, len(art), len(neg)))
    for s in stats:
        a = [fl(p[s]) for p in art]; b = [fl(p[s]) for p in neg]
        print('      %-12s AUC %.3f   artefact min %.3f median %.3f | %s max %.3f median %.3f' % (
            s, auc(a, b), min(a) if a else float('nan'), sorted(a)[len(a) // 2] if a else float('nan'), name.split()[0],
            max(b) if b else float('nan'), sorted(b)[len(b) // 2] if b else float('nan')))

print('\n## what a pose gate would do: pairs with z_dflip >= Z, by label (kept/leaves decisions only = the tree delivered them)')
for Z in (3, 5, 8, 12):
    sel = [p for p in pairs if fl(p['z_dflip']) >= Z]
    kept = [p for p in sel if p['decision'] in ('kept', 'leaves')]
    print('   Z >= %2d: %3d pairs flagged (%s); delivered: %s' % (
        Z, len(sel), ', '.join('%s %d' % (k, sum(p['label'] == k for p in sel)) for k in ('artefact', 'suspect', 'needed', 'unneeded', 'ole', 'production', 'unknown')),
        ', '.join('%s %d' % (k, sum(p['label'] == k for p in kept)) for k in ('needed', 'unneeded', 'ole'))))
    for p in kept:
        if p['label'] in ('needed', 'unneeded', 'ole'):
            print('        %s %s %s n %s/%s flips %s/%s z %s (%s)' % (p['tree'], p['node'], p['pair'], p['n_a'], p['n_b'], p['flip_a'], p['flip_b'], p['z_dflip'], p['label']))

print('\n## real PfCRT and 692k pairs')
for p in pairs:
    if p['label'] in ('suspect', 'production') and (fl(p['z_dflip']) >= 3 or fl(p['dflip']) >= 0.05):
        print('   %s %s %s %s n %s/%s flips %s/%s dflip %s z %s ks %s z %s view_tv %s' % (
            p['tree'], p['node'], p['kind'], p['pair'], p['n_a'], p['n_b'], p['flip_a'], p['flip_b'], p['dflip'], p['z_dflip'], p['ks_ang'], p['z_ks_ang'], p['view_tv']))
