#!/usr/bin/env python3
"""Phase-3 side by side per site from score_tree_sites.py JSON scorecards (tree leaves, site-layer combination classes, any other
labelling): NMI(state; class), best-class purity and F80 of each state, the nodes splitting the site (tree), combination coverage and
the leaf-map occupancy fidelity; optional oracle (oracle_sites.py output lines are parsed for the AUC per site and band).
usage: compare_tree_sites.py NAME:score.json [NAME:score.json ...] [--oracle oracle.out]"""
import sys, json, re
args = [a for a in sys.argv[1:] if not a.startswith('--')]; orc = sys.argv[sys.argv.index('--oracle') + 1] if '--oracle' in sys.argv else ''
if orc in args: args.remove(orc)
S = [(a.split(':', 1)[0], json.load(open(a.split(':', 1)[1]))) for a in args]; sites = list(S[0][1]['sites'].keys())
orac = {}
if orc:
    for l in open(orc):
        m = re.match(r'lp\s+(\S+) site (\S+)\s+state-1 share (\S+): d\' (\S+), AUC (\S+), Bayes accuracy (\S+)', l)
        if m: orac.setdefault(m.group(2), []).append('lp %s d\' %s AUC %s' % (m.group(1), m.group(4), m.group(5)))
for nm, sc in S: print('== %-14s %s: %d classes (%d pending)%s' % (nm, sc.get('state', ''), sc['leaves'], sc.get('pending', 0), '; coverage %d of %d combinations >= 5k' % (len(sc['coverage']['covered']), len(sc['coverage']['big_combos'])) if 'coverage' in sc else ''))
for s in sites:
    print('-- site %s%s' % (s, ('   oracle: ' + '; '.join(orac[s])) if s in orac else ''))
    for nm, sc in S:
        r = sc['sites'][s]; st = ['state %s: purity %.2f F80 %d' % (k[-1], r[k]['best_purity'], r[k]['F80']) for k in ('state1', 'state0') if k in r]
        mp = r.get('map', {}); sp = r.get('split_nodes', [])
        print('   %-14s NMI %.3f eta2 %.3f | %s%s%s' % (nm, r['nmi'], r['eta2'], ' | '.join(st), (' | split at %d node(s)' % len(sp)) if sp else '', (' | map occupancy vs share corr %.2f rmse %.2f' % (mp['occ_vs_share_corr'], mp['occ_rmse'])) if mp else ''))
