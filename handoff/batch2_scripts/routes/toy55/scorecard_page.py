#!/usr/bin/env python3
"""usage: scorecard_page.py OUT TREE_DIR NAME KIND   (score_tree.sbatch's last step; writes OUT/scorecard.txt and prints it)
The one-page usefulness scorecard of a hetero_tree run from score_tree.sbatch's outputs in OUT:
- tree shape (leaves, junk-named leaves, units, rows in leaves, nodes, depth, splits) and cost (sacct: wall, core-h, jobs);
- A (rnasep): the leaves against the 55 deposited apo conformers: recall by the paired arm / core test, mutual best matches,
  closed-state recovery (share of particles in leaves with closed-arm >= 0.6 vs the deposited 11.6 %), arm span (s1 vs -10.1 ..
  +7.5 A), leaves outside the deposited set, per-leaf FSC 0.143, quality-confounded leaves (corr AUC vs the rest <= 0.35);
- B: every split's effect against half-map noise, flagged
    n.s. (z < 5) | junk (a junk child) | blur (splits.txt: dB >= 1/2 the root's quality reference at the same box, >= 3 SE) |
    arm (rnasep: abs ds1 >= 3 A or abs dclosed >= 0.10) | minor (rnasep: significant, none of these) | effect (other kinds);
- C (toys): the truth recovery SUMMARY of validate_tree3.py / validate_conf.py; D (a seed pair): tree_pair_verify.py's agreement lines."""
import os, re, sys, numpy as np
from datetime import datetime
OUT, T, NAME, KIND = sys.argv[1:5]; T = T.rstrip('/'); P = []
def p(s=''): P.append(s)
def rd(f): return open(os.path.join(OUT, f)).read() if os.path.exists(os.path.join(OUT, f)) else ''
LV = [l.split() for l in open(T + '/hetero_tree_leaves.txt') if l.strip()]; L = np.loadtxt(T + '/hetero_tree_labels.txt', dtype=int)
nodes = [l.split()[0] for l in open(T + '/hetero_tree_nodes.txt') if l.strip()]; ui = T + '/final_maps/units_index.txt'
nunits = sum(1 for l in open(ui) if l.strip() and not l.startswith('#')) if os.path.exists(ui) else 0
SPL = [l.split() for l in rd('splits.txt').splitlines() if l.strip() and not l.startswith('#')]
p('SCORECARD %s (%s), %s' % (NAME, KIND, datetime.now().strftime('%Y-%m-%d %H:%M')))
if rd('note.txt').strip(): p('note    ' + rd('note.txt').strip())
p('tree    %d leaves (%d named junk), %d units, %s rows in leaves; %d nodes, depth %d; %d splits scored' % (len(LV), sum(w[1] == 'junk' for w in LV), nunits,
  '{:,}'.format(int((L > 0).sum())), len(nodes), max(len(x) for x in nodes) - 1, len(set(w[0] for w in SPL))))
GQ = [l[5:] for l in rd('qtoyscore.txt').splitlines() if l.startswith('QTOY CLOSED GATE')]; GR = re.search(r'^closed: (\d+) of (\d+) recovered', rd('recall.txt'), re.M)
if GQ: p('GATE    closed state ' + GQ[0][len('CLOSED GATE: '):])
elif GR: p('GATE    closed state: %s of %s closed deposited conformers recalled (paired arm / core test)' % GR.groups())
J = [l.split('|') for l in rd('jobs.txt').splitlines() if l.count('|') >= 6]
if J:
    def dt(s):
        try: return datetime.strptime(s, '%Y-%m-%dT%H:%M:%S')
        except ValueError: return None
    st = [dt(j[5]) for j in J if dt(j[5])]; en = [dt(j[6]) for j in J if dt(j[6])]; ch = sum(float(j[3]) * float(j[4]) for j in J if j[3].isdigit() and j[4].isdigit()) / 3600.
    p('cost    wall %.1f h, %s core-h over %d jobs (sacct, jobs run in the run dir)' % ((max(en) - min(st)).total_seconds() / 3600. if st and en else float('nan'), '{:,.0f}'.format(ch), len(J)))
LS = {}
for l in rd('leaves.txt').splitlines():
    m = re.search(r'LEAF (\S+) (\S+) \| rows (\d+) \| qAUC (\S+) \| median dfx (\S+) \| top-view share (\S+) \| FSC0.143 (\S+) A', l)
    if m: LS[m.group(1)] = dict(kind=m.group(2), n=int(m.group(3)), q=float(m.group(4)), res=float(m.group(7)))
if KIND == 'rnasep':
    p(); p('A. AGAINST THE 55 DEPOSITED APO CONFORMERS (every leaf docked on the rigid core first)')
    hd, R = None, {}
    for l in rd('fa_leaves/maps_table.txt').splitlines():
        w = l.split()
        if l.startswith('#'): hd = w[1:]; continue
        if hd and w: R[w[0]] = dict(zip(hd, w))
    rec = re.findall(r'^(closed|intermediate|open): (\d+) of (\d+) recovered', rd('recall.txt'), re.M)
    mu = re.search(r'mutual best matches[^:]*: (.*)', rd('fa_leaves.log')); mu = mu.group(1).split() if mu else []
    if rec: p('   recalled (paired arm / core test): %d of %d conformers (%s); mutual best matches %d' % (sum(int(r[1]) for r in rec), sum(int(r[2]) for r in rec),
             ', '.join('%s %s/%s' % r for r in rec), len(mu)))
    lv = {k: v for k, v in R.items() if v.get('kind') in ('leaf', 'junk')}; tot = sum(int(v['particles']) for v in lv.values())
    if lv:
        cl = [k for k, v in lv.items() if float(v['closed_arm']) >= 0.6]; ncl = sum(int(lv[k]['particles']) for k in cl)
        p('   closed state: %.1f %% of the leaves\' particles in %d leaves with closed-arm >= 0.6 (deposited set 11.6 %%)' % (100. * ncl / max(tot, 1), len(cl)))
        s1 = {k: float(v['s1_A']) for k, v in lv.items()}; ks = sorted(s1, key=s1.get)
        p('   arm span: s1 %+.1f (%s) .. %+.1f A (%s) over the leaves (deposited -10.1 .. +7.5 A)' % (s1[ks[0]], ks[0], s1[ks[-1]], ks[-1]))
        out = [k for k, v in lv.items() if v.get('paired_outside') == '1']
        p('   outside the deposited set (paired test): %d of %d leaves, %s particles' % (len(out), len(lv), '{:,}'.format(sum(int(lv[k]['particles']) for k in out))))
    if LS:
        rs = np.array([v['res'] for v in LS.values() if np.isfinite(v['res'])]); lo = sorted([(v['n'], k, v['q']) for k, v in LS.items() if v['q'] <= 0.35], reverse=True)
        p('   per-leaf FSC 0.143 (soft root mask): median %.1f A, worst %.1f A' % (np.median(rs), rs.max()) if len(rs) else '   per-leaf FSC: n/a')
        p('   quality: %d leaves (%s particles) at corr AUC <= 0.35 vs the rest at matched defocus%s' % (len(lo), '{:,}'.format(sum(x[0] for x in lo)),
          ('; largest ' + ', '.join('%s %s (%.2f)' % (k, '{:,}'.format(n), q) for n, k, q in lo[:4])) if lo else ''))
elif LS:
    rs = np.array([v['res'] for v in LS.values() if np.isfinite(v['res'])]); lo = [k for k, v in LS.items() if v['q'] <= 0.35]
    p(); p('leaves  FSC 0.143 median %.1f A, worst %.1f A; %d leaves at corr AUC <= 0.35 vs the rest at matched defocus' % (np.median(rs), rs.max(), len(lo)) if len(rs) else 'leaves  n/a')
if SPL:
    p(); p('B. SPLITS: each group\'s stage-B map vs the largest group\'s, 15 A, the node\'s soft mask')
    RF = [(float(w[24]), w) for w in SPL if len(w) > 26 and np.isfinite(float(w[24]))]
    p('   flags: arm (abs ds1 >= 3 A or abs dclosed >= 0.10) | blur = quality: the lower-corr group blurrier where the map does not change')
    p('   (unfiltered halves) by at least 1/2 of the root\'s own top vs bottom corr quintile split at the same box, and >= 3 SE')
    p('   (validated: quality 0.51..0.64 of the reference, real states <= 0.18) | minor (local change under 0.028, the nearest-neighbour')
    p('   deposited pair median) | other (none of these) | junk | n.s. (z < 5). qAUC (corr at matched defocus) is shown, not used.')
    p('   root quality reference: dB_ref %s' % ('%+.0f A^2 (stage-B box, every node inherits it)' % RF[0][0] if RF else 'n/a: blur not judged'))
    rows = []
    for w in SPL:
        X, g, n, chd, ref, nr, chr_, z, snr, share, loc, cos, pm, tx, ty, tz, tr, q, ds1, dcl, dco = w[:21]; z, share, q, ds1, dcl, loc = (float(x) for x in (z, share, q, ds1, dcl, loc))
        dB = float(w[21]) if len(w) > 21 else np.nan; dBse = float(w[23]) if len(w) > 23 else np.nan; rt = float(w[25]) if len(w) > 26 else np.nan; f = []
        if z < 5: f.append('n.s.')
        else:
            if 'junk' in chd or 'junk' in chr_: f.append('junk')
            if KIND == 'rnasep' and np.isfinite(ds1) and (abs(ds1) >= 3 or abs(dcl) >= 0.10): f.append('arm')
            if len(w) > 26 and w[26] == '1': f.append('blur')
            if not f: f.append('minor' if (KIND == 'rnasep' and loc < 0.028) else 'other')
        rows.append((X, chd, int(n), chr_, int(nr), z, share, loc, q, ds1, dcl, dB, rt, dBse, '+'.join(f)))
    tags = ['n.s.', 'junk', 'arm', 'blur', 'minor', 'other']; c = {k: sum(k in r[-1].split('+') for r in rows) for k in tags}
    p('   %d compared pairs: %s' % (len(rows), ', '.join('%s %d' % (k, v) for k, v in c.items() if v)))
    qs = rd('qtoyscore.txt').splitlines(); TQ = {}; hq = [i for i, l in enumerate(qs) if l.startswith('# split node')]
    for l in (qs[hq[0] + 1:] if hq else []):
        if l.startswith('#'): break
        w = l.split(); TQ[w[0]] = (float(w[4]), w[-1])
    if TQ:
        qn = {x for x in TQ if TQ[x][1] == 'quality'}; bn = {r[0] for r in rows if 'blur' in r[-1].split('+')}
        p('   truth (qtoy) per split, call and eta^2 of log q: %s; quality-following splits %d, blur-flagged %d of them; blur flags on state splits %d'
          % (', '.join('%s %s %.2f' % (x, TQ[x][1], TQ[x][0]) for x in sorted(TQ, key=lambda s: (len(s), s))), len(qn), len(qn & bn), len(bn - qn)))
    bl = [r for r in rows if 'blur' in r[-1].split('+')]
    if bl: p('   blur splits (dB / dB_ref): ' + ', '.join('%s %s %.2f (%+.0f +- %.0f)' % (r[0], r[1].replace('node_', ''), r[12], r[11], r[13]) for r in sorted(bl, key=lambda r: -r[12])))
    p('   %-10s %-22s %-22s %5s %6s %6s %6s %6s %7s %6s %6s  %s' % ('node', 'group', 'vs (largest)', 'z', 'share', 'loc', 'qAUC', 'ds1', 'dclosed', 'dB', 'dB/ref', 'flag'))
    for r in sorted(rows, key=lambda r: (len(r[0]), r[0]))[:24]:
        p('   %-10s %-22s %-22s %5.1f %6.3f %6.3f %6.3f %+6.2f %+7.3f %+6.0f %+6.2f  %s' % (r[0], '%s %s' % (r[1].replace('node_', ''), '{:,}'.format(r[2])), '%s %s' % (r[3].replace('node_', ''), '{:,}'.format(r[4])),
          r[5], r[6], r[7], r[8], r[9], r[10], r[11], r[12], r[14]))
    if len(rows) > 24: p('   (%d more in splits.txt)' % (len(rows) - 24))
TR = [l for l in rd('truth.txt').splitlines() if l.startswith('SUMMARY') or 'closed' in l.lower() and l.startswith('   leaf')]
AR = [l[5:] for l in rd('armscore.txt').splitlines() if l.startswith('ARMS ')] + [l[5:] for l in rd('qtoyscore.txt').splitlines() if l.startswith('QTOY ')]
if TR or AR: p(); p('C. TRUTH RECOVERY'); [p('   ' + l[:220]) for l in TR[:6]]; [p('   ' + l[:700]) for l in AR]
PR = [l for l in rd('pair.txt').splitlines() if re.search(r'ARI|NMI|same|agree', l)]
if PR: p(); p('D. SEED REPRODUCIBILITY (tree_pair_verify.py)'); [p('   ' + l.strip()[:200]) for l in PR[:6]]
open(os.path.join(OUT, 'scorecard.txt'), 'w').write('\n'.join(P) + '\n'); print('\n'.join(P))
