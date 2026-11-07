#!/usr/bin/env python3
"""ppinfo_table.py: the lever tables from the ppinfo runs. usage: ppinfo_table.py RUNS_DIR [RUN ...]
Run names are <config>_<site> (config = the lever, site = L065s1 L065s2 L08s1 L08s2 L10s1 qtoynull olenull). Prints (1) one row per
endpoint (TSV: run, endpoint, N*, share, r, z, SE, bar, pole, AUC, class >= 0.5 n / P / R, top-N* P / R, iterations) and (2) the
config x site matrix: the best rare-ward endpoint's z (over endpoints with class share <= 0.5) with its r and AUC, the number of
accepted poles, and for the nulls the poles and the largest z."""
import glob, json, os, sys

root = sys.argv[1]
runs = sys.argv[2:] or sorted(os.path.basename(d) for d in glob.glob(os.path.join(root, '*')) if os.path.exists(os.path.join(d, 'summary.json')))
SITES = ['L065s1', 'L065s2', 'L08s1', 'L08s2', 'L10s1', 'qtoynull', 'olenull']
rows, cell = [], {}
print('\t'.join(['run', 'endpoint', 'Nstar', 'share', 'r', 'z', 'SE', 'bar_r', 'neff', 'accept', 'seeds', 'seed_z', 'kept', 'AUC', 'hi_n', 'hi_P', 'hi_R', 'top_P', 'top_R', 'its', 'box', 'kb', 'khi']))
for run in runs:
    fn = os.path.join(root, run, 'summary.json')
    if not os.path.exists(fn):
        continue
    S = json.load(open(fn))
    site = next((s for s in SITES if run.endswith('_' + s)), run.rsplit('_', 1)[-1])
    cfg = run[:len(run) - len(site) - 1]
    best = None
    npole = 0
    zmax = -9.
    for e in S['endpoints']:
        has = 'auc' in e
        print('\t'.join([run, e['tag'], str(e['nstar']), '%.3f' % e['share'], '%.3f' % e['r'], '%.2f' % e['z'], '%.3f' % e['se'], '%.3f' % S['bar_r'],
                         '%.1f' % S['neff'], 'yes' if e['accepted'] else '-', str(e.get('seed_n', '-')), '%.2f' % e['seed_z'] if 'seed_z' in e else '-', 'KEPT' if e.get('kept') else '-',
                         '%.3f' % e['auc'] if has else '-', str(e['hi_n']) if has else '-', '%.3f' % e['hi_P'] if has else '-', '%.3f' % e['hi_R'] if has else '-',
                         '%.3f' % e['top_P'] if has else '-', '%.3f' % e['top_R'] if has else '-', str(e['nit']), str(S['box']), str(S['kb']), str(S['khi'])]))
        npole += int(e.get('kept', e['accepted']))
        if e['share'] <= 0.5:
            zmax = max(zmax, e['z'])
            if has:
                key = (e['auc'] >= 0.5, e['z'])   # the rare-ward endpoints first, then the largest z
                if best is None or key > best[0]:
                    best = (key, e)
            elif best is None or e['z'] > best[1]['z']:
                best = ((True, e['z']), e)
    cell[(cfg, site)] = (best[1] if best else None, npole, zmax, S)
cfgs = sorted({c for c, _ in cell}, key=lambda c: (c != 'base', c))
print('\n# config x site: best rare-ward endpoint z (r, AUC, class>=0.5 purity) | kept poles (accepted endpoint with >= 6 seeds whose deviation reproduces at z >= 3; the acceptance alone where the run predates the seed test); nulls: kept poles / max z')
hdr = ['config'] + SITES
print('\t'.join(hdr))
bar = ['bar_r'] + ['%.3f' % cell[(c, s)][3]['bar_r'] if (c, s) in cell else '-' for c in cfgs[:1] for s in SITES]
print('\t'.join(bar))
for c in cfgs:
    out = [c]
    for s in SITES:
        if (c, s) not in cell:
            out.append('-'); continue
        e, npole, zmax, S = cell[(c, s)]
        if e is None:
            out.append('no endpoint'); continue
        if 'auc' in e:
            out.append('z %.2f (r %.3f, AUC %.3f, P %.2f n %d) | %d pole%s' % (e['z'], e['r'], e['auc'], e['hi_P'], e['hi_n'], npole, '' if npole == 1 else 's'))
        else:
            out.append('%d pole%s / max z %.2f' % (npole, '' if npole == 1 else 's', zmax))
    print('\t'.join(out))
