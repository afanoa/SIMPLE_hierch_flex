#!/usr/bin/env python3
"""old vs new site_rank verdicts on the same inputs: for each site name present in both roots, the poles (site_poles.txt),
the kept poles with their margins (site_pole_rule.txt) and the wall time (rerun_manifest.txt).
usage: compare_sites.py OLD_ROOT NEW_ROOT"""
import sys, os, re
old, new = sys.argv[1:3]
def poles(d):
    p = os.path.join(d, 'site_poles.txt')
    if not os.path.exists(p): return '(none)'
    return ' ; '.join(l.strip() for l in open(p) if l.strip()) or '(empty)'
def rule(d):
    p = os.path.join(d, 'site_pole_rule.txt'); out = []
    if not os.path.exists(p): return '-'
    for l in open(p):
        w = l.split()
        if not w or l.startswith('#'): continue
        out.append(f'p{w[0]} share {w[2]} z {w[3]}/{w[4]} {w[-1]}')
    return '; '.join(out) or '-'
def wall(d):
    p = os.path.join(d, 'rerun_manifest.txt')
    if not os.path.exists(p): return ''
    m = re.search(r'wall ([0-9.]+) s', open(p).read()); return m.group(1) + ' s' if m else ''
def nfits(d):
    p = os.path.join(d, 'site_rank.log')
    if not os.path.exists(p): p = os.path.join(d, 'hetero_step.log')
    if not os.path.exists(p): return ''
    t = open(p, errors='replace').read()
    return f'{t.count("site EM: 2 classes")} EM + {t.count("site ranking by a given score")} scoring'
names = sorted(set(os.listdir(old)) & set(os.listdir(new)))
agree = 0
for nm in names:
    a, b = os.path.join(old, nm), os.path.join(new, nm)
    pa, pb = poles(a), poles(b)
    same = (pa == pb)
    agree += same
    print(f'## {nm}  [{"same" if same else "DIFF"}]  old {wall(a)} {nfits(a)} | new {wall(b)} {nfits(b)}')
    print(f'   old poles: {pa}\n   old rule : {rule(a)}')
    print(f'   new poles: {pb}\n   new rule : {rule(b)}')
print(f'\n{agree} of {len(names)} sites give the same site_poles.txt')
