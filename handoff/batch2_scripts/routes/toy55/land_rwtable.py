#!/usr/bin/env python3
"""land_rwtable.py: tables of a land_reweight.py synthetic run: per setting (variant x temperature) and penalty, mean TV, every shape's
barrier estimate (true in the header), the shoulder's minor mass, and the penalty minimising the mean |barrier error| over dw2, dw4, flat
and broad (no false barrier allowed on flat: its estimate must stay below +0.3 kT). usage: land_rwtable.py REWEIGHT.json [more.json]"""
import sys, json
for f in sys.argv[1:]:
    R = json.load(open(f)); print('== %s (pixel-split tau %.3f)' % (f, R['tau_split']))
    for key, tab in R['table'].items():
        shapes = list(tab); lams = list(tab[shapes[0]]); tb = {s: R['chosen'][key][s]['truth']['barrier'] for s in shapes}
        print('  %s (tau %.3f); true barriers %s; shoulder minor %.3f' % (key, R['chosen'][key]['tau'], ' '.join('%s %+.2f' % (s, tb[s]) for s in shapes),
              R['chosen'][key]['shoulder']['truth']['minor'] if 'shoulder' in shapes else float('nan')))
        best = None
        for lm in lams:
            tv = sum(tab[s][lm]['tv'] for s in shapes) / len(shapes); be = [abs(tab[s][lm]['barrier'] - tb[s]) for s in shapes if s in ('dw2', 'dw4', 'flat', 'broad')]
            mbe = sum(be) / max(len(be), 1); ok = tab.get('flat', {}).get(lm, {}).get('barrier', 0.) < 0.3
            if ok and (best is None or mbe < best[1]): best = (lm, mbe)
            print('    lam %-6s TV %.3f | barriers %s | shoulder minor %.3f | mean |barrier err| %.2f%s' % (lm, tv, ' '.join('%s %+.2f' % (s, tab[s][lm]['barrier']) for s in shapes),
                  tab['shoulder'][lm]['minor'] if 'shoulder' in shapes else float('nan'), mbe, '' if ok else ' (false barrier on flat)'))
        print('    -> barrier-optimal lam %s (mean |barrier err| %.2f kT); TV-chosen lam %g' % (best[0], best[1], R['chosen'][key]['lam']))
