#!/usr/bin/env python3
"""Merge up the tree, after one level's split tests: every tested node's surviving groups (children_A/B of splitbench8,
cross-half EM labels, one numbering for both halves) become the new units; a unit that did not survive is merged into the
group that took most of its particles; untested nodes (one unit, or a failed test) pass through. SRC = the level's input:
a consolidation dir (states.txt + groups.txt, zoom level) or the previous level's OUT (units.txt + units_A/B.txt).
Writes OUT/units.txt (unit fam reg leaf nA nB members closed_arm_max orig_states), OUT/units_{A,B}.txt (per half row);
ARM = arm_states.txt of the source consolidation (closed-arm occupancy of its states).
usage: mu_collect.py LEVEL_DIR SRC MRA MRB OUT [ARM]"""
import sys, os, re, numpy as np
ld, src, mra, mrb, out = sys.argv[1:6]; os.makedirs(out, exist_ok=True); arm = {}
if len(sys.argv) > 6 and os.path.isfile(sys.argv[6]):
    for l in open(sys.argv[6]):
        p = l.split()
        if p and not l.startswith('#') and 'state_' in p[0]: arm[int(re.findall(r'state_(\d+)', p[0])[0])] = float(p[1])
name, orig = {}, {}
if os.path.isfile(src + '/groups.txt'):   # zoom level: units = consolidated states
    st = np.loadtxt(src + '/states.txt', dtype=int); MR = {h: np.loadtxt(f).astype(int) - 1 for h, f in (('A', mra), ('B', mrb))}; prev = {h: st[MR[h]] for h in 'AB'}
    for l in open(src + '/groups.txt'):
        if not l.startswith('#') and l.strip(): w = l.split(); name[int(w[0])] = w[2]; orig[int(w[0])] = [int(w[0])]
else:
    prev = {h: np.loadtxt('%s/units_%s.txt' % (src, h), dtype=int) for h in 'AB'}
    for l in open(src + '/units.txt'):
        if not l.startswith('#') and l.strip(): w = l.split(); name[int(w[0])] = w[6]; orig[int(w[0])] = [int(x) for x in w[8].split(',')]
new = {h: np.zeros_like(prev[h]) for h in 'AB'}; rows = []; nid = 0
def key(d):   # node dir name -> (fam, reg, leaf)
    b = os.path.basename(d); m = re.match(r'^fam(\d)_zoomregion(.*)_s(\d+)$', b) or re.match(r'^fam(\d)_region(.*)$', b) or re.match(r'^fam(\d)_(family)$', b)
    return (int(m.group(1)), m.group(2), int(m.group(3)) if m.lastindex == 3 else 0)
def add(f, r, lf, sel, members):
    global nid; nid += 1
    for h in 'AB': new[h][sel[h]] = nid
    rows.append((nid, f, r, lf, int(sel['A'].sum()), int(sel['B'].sum()), '+'.join(name[u] for u in members) or 'reassigned', sorted(x for u in members for x in orig[u])))
for l in open(ld + '/nodes.txt'):
    if l.startswith('#') or not l.strip(): continue
    d, mask, us = l.split(); us = [int(u) for u in us.split(',')]; f, r, lf = key(d)
    if len(us) < 2 or not all(os.path.isfile('%s/children_%s.txt' % (d, h)) for h in 'AB'):   # pass through
        for u in us: add(f, r, lf, {h: prev[h] == u for h in 'AB'}, [u])
        if len(us) >= 2: print('%s: test missing, %d units passed through' % (os.path.basename(d), len(us)))
        continue
    lab = {int(x.split()[0]): int(x.split()[1]) for x in open(d + '/labels.txt') if not x.startswith('#')}
    ch = {h: np.loadtxt('%s/children_%s.txt' % (d, h), dtype=int) for h in 'AB'}; cand = {h: np.loadtxt('%s/cand_%s.txt' % (d, h), dtype=int) for h in 'AB'}
    K = int(max(ch['A'].max(), ch['B'].max())); mem = {g: [] for g in range(1, K + 1)}
    for c, u in lab.items():   # a unit belongs to the group holding most of its particles (both halves)
        cnt = np.bincount(np.r_[ch['A'][cand['A'] == c], ch['B'][cand['B'] == c]], minlength=K + 1); cnt[0] = 0
        if cnt.sum(): mem[int(cnt.argmax())].append(u)
    for g in range(1, K + 1): add(f, r, lf, {h: ch[h] == g for h in 'AB'}, mem[g])
    print('%s: %d units -> %d' % (os.path.basename(d), len(us), K))
for h in 'AB': np.savetxt('%s/units_%s.txt' % (out, h), new[h], fmt='%d')
with open(out + '/units.txt', 'w') as fo:
    fo.write('# unit fam reg leaf nA nB members closed_arm_max orig_states\n')
    for i, f, r, lf, na, nb, nm, og in rows:
        a = max([arm[x] for x in og if x in arm], default=float('nan')); fo.write('%d %d %s %d %d %d %s %.3f %s\n' % (i, f, r, lf, na, nb, nm, a, ','.join(str(x) for x in og) or '0'))
print('%d units after this level (%d particles in half A, %d in B)' % (len(rows), int((new['A'] > 0).sum()), int((new['B'] > 0).sum())))
