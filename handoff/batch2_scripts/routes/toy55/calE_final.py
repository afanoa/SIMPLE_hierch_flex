#!/usr/bin/env python3
"""usage: calE_final.py RELF_FILE MAPTESTS_FILE.. [--nbsw 0.030] [--relmin X] [--res 15] [--picks PICKS ..]
Joins, per series, phaseE_maptests.py's SUMMARY line and SF lines with phaseF_rel.py's RELF line, then applies the at-node rule:
1. quality: the worst bin's FSC 0.143 must be at least as fine as the test resolution, else 'insufficient';
2. none: rel_loc (15 A) < RELMIN;
3. swing: nb_SF (15 A, max over bins 2-4) >= NBSW and the SF residual z of some middle bin >= 5;
4. else switch, with the flags graded (dominant step share < 0.40) and dipole (1/3 < +/- < 3).
Truth class from the name: null_* and cdisc_B_* are nulls, m52M_* / m8c17_* swings, CAHRA / cdisc switches, m33_* / ole_* real data.
Prints one row per series, then each class's calls and the per-statistic class ranges."""
import re, sys, numpy as np
a = sys.argv[1:]; opt = {}
for k in ('--nbsw', '--relmin', '--res'):
    if k in a: i = a.index(k); opt[k] = float(a[i + 1]); a = a[:i] + a[i + 2:]
picks = {}
if '--picks' in a:
    i = a.index('--picks'); pf = a[i + 1:]; a = a[:i]
    for f in pf:
        tag = 'm8c17' if 'm8c17' in f else 'm52M'
        for l in open(f):
            w = l.split()
            if w and w[0] == 'PICK': picks['%s_%s_p%s_%s' % (tag, w[2], w[3], w[1].replace('_', ''))] = abs(float(w[4]))
NBSW, RELMIN, TRES = opt.get('--nbsw', 0.030), opt.get('--relmin', 0.0), opt.get('--res', 15.)
RELF, MT = a[0], a[1:]; S = {}
def fl(s): return [float(x) for x in s.split('/')]
for f in MT:
    cur = None
    for l in open(f):
        if l.startswith('######## '): cur = l.split()[1].rstrip(':'); S.setdefault(cur, {})
        elif cur and ' 15 A SF' in l:
            S[cur]['esnrF'] = float(re.search(r'ends SNR ([\d.]+)', l).group(1)); zz = [float(x) for x in re.findall(r'\(z (-?[\d.]+)\)', l)]; S[cur]['zres'] = max(zz[1:]) if len(zz) > 1 else np.nan
        elif l.startswith('SUMMARY'):
            nm = l.split()[1]; v = S.setdefault(nm, {}); v['rows'] = int(re.search(r'min rows/bin (\d+)', l).group(1))
            v['pm'] = float(re.search(r'\+/- ([\d.]+)', l).group(1))
            for res, blk in zip((15, 10), re.split(r'\| (?=1[05] A:)', l)[1:]):
                m = re.search(r'steps \S+ share ([\d.]+) \| S (\S+) nb (\S+) \| SF (\S+) nb (\S+)', blk)
                v[res] = dict(share=float(m.group(1)), R=max(fl(m.group(2))), nb=max(fl(m.group(3))), nbF=max(fl(m.group(5))))
for l in open(RELF):
    if not l.startswith('RELF'): continue
    nm = l.split()[1]; v = S.setdefault(nm, {}); v['fsc'] = float(re.search(r'FSC 0.143 ([\d.]+)', l).group(1))
    for res, blk in zip((15, 10), re.split(r'\| (?=1[05] A:)', l)[1:]):
        m = re.search(r'rel ([-\d.]+) cos ([-+\d.]+) rel_loc ([-\d.]+) z_loc ([-\d.]+) hmin ([-\d.]+)', blk); v['rf%d' % res] = [float(x) for x in m.groups()]
def cls(nm):
    if nm.startswith('null_') or nm.startswith('cdisc_B'): return 'null'
    if nm.startswith('m33_') or nm.startswith('ole_'): return 'real'
    return 'swing' if ('m52M_' in nm or 'm8c17_' in nm) else 'switch'
def call(v):
    if 'fsc' not in v or 15 not in v: return 'n/a'
    if v['fsc'] > TRES: return 'insufficient'
    if v['rf15'][2] < RELMIN: return 'none'
    if v[15]['nbF'] >= NBSW and v.get('zres', 0) >= 5: return 'swing'
    fl_ = [x for x, c in (('graded', v[15]['share'] < 0.40), ('dipole', 1 / 3. < v['pm'] < 3)) if c]
    return 'switch' + ('(%s)' % ','.join(fl_) if fl_ else '')
print('%-27s %-6s %6s %5s | %6s %6s %7s %6s %5s | %6s %6s %6s %5s %5s | %s' % ('series', 'class', 'rows', 'fsc', 'rel', 'cos', 'rel_loc', 'z_loc', 'hmin', 'R15', 'nbF15', 'nbF10', 'zres', 'step', 'call'))
rows = sorted(S, key=lambda n: ({'switch': 0, 'swing': 1, 'null': 2, 'real': 3}[cls(n)], n))
for nm in rows:
    v = S[nm]
    if 15 not in v or 'fsc' not in v: print('%-27s %-6s (incomplete: %s)' % (nm, cls(nm), ', '.join(sorted(map(str, v)))) ); continue
    r = v['rf15']; print('%-27s %-6s %6d %5.1f | %6.4f %+6.3f %7.4f %6.1f %5.3f | %6.2f %6.3f %6.3f %5.1f %5.2f | %s%s' % (nm, cls(nm), v['rows'], v['fsc'], r[0], r[1], r[2], r[3], r[4], v[15]['R'],
          v[15]['nbF'], v[10]['nbF'], v.get('zres', np.nan), v[15]['share'], call(v), '  (pick %.2f)' % picks[nm] if nm in picks else ''))
print('\ncalls per class (rule: fsc <= %.0f A, rel_loc >= %.4f, swing at nb_SF >= %.3f with residual z >= 5):' % (TRES, RELMIN, NBSW))
for c in ('switch', 'swing', 'null', 'real'):
    cs = [call(S[n]).split('(')[0] for n in rows if cls(n) == c and 15 in S[n] and 'fsc' in S[n]]
    print('  %-6s %2d: %s' % (c, len(cs), ', '.join('%s %d' % (k, cs.count(k)) for k in sorted(set(cs)))))
print('\nclass ranges:')
for key, get in (('fsc', lambda v: v['fsc']), ('rel_loc15', lambda v: v['rf15'][2]), ('cos15', lambda v: v['rf15'][1]), ('z_loc15', lambda v: v['rf15'][3]), ('nbF15', lambda v: v[15]['nbF']),
                 ('nbF10', lambda v: v[10]['nbF']), ('zres', lambda v: v.get('zres', np.nan)), ('R15', lambda v: v[15]['R'])):
    print('  %-9s ' % key + ' | '.join('%s %s' % (c, '%.4g-%.4g' % (min(x), max(x)) if x else '-') for c in ('switch', 'swing', 'null', 'real')
          for x in [[get(S[n]) for n in rows if cls(n) == c and 15 in S[n] and 'fsc' in S[n]]]))
