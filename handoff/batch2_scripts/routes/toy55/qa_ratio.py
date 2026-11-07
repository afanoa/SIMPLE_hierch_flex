#!/usr/bin/env python3
"""usage: qa_ratio.py cases|table TREE_DIR KIND NODES OUTDIR REFPREFIX [EXCL]
Option (a) of the quality screen in each node's own frame (recB0/b.simple: node_motion's poses and corr). For a node_motion candidate
(motion/cK: bins.txt, recvol_state01 / 05 with halves and _unfil) the lower-corr end bin L (1 or 5) and its corr gap come from the
median corr differences inside dfx deciles, averaged over the deciles; the reference is quality_ref.sbatch's top vs bottom corr quintile
of the same project at the same box (REFPREFIX<node>/: dump.npz, states.txt, dBref.txt), its gap measured the same way.
Expected pure-quality blur = dB_ref x gap / reference gap; ratio = dB / expected; flag at ratio >= 1 and dB >= 3 SE.
NODES: comma list (r,r1,..); KIND cdisc | igm | rnasep | other; EXCL: comma list of node:cand to skip (duplicate bins).
cases: OUTDIR/cases.txt + qscreen.txt (quality_screen3.py: L's blur vs H from the unfiltered halves) + cases_info.tsv, and for rnasep
OUTDIR/dock_maplist.txt (fa_maps10.py --dock on both end bins).
table: + OUTDIR/qs3.txt (rnasep: OUTDIR/dock/maps_table.txt) -> OUTDIR/ratios.txt. Labels: cdisc / igm from the truth composition of the
end bins (L minus H: junk share; cdisc A / B / C site shares, igm arm-angle eta^2 over the genuine rows; supported >= 0.2 share or
0.1 eta^2, weak >= 0.1 / 0.03); rnasep from the docked end bins (real:arm at abs ds1 >= 3 A or abs dclosed >= 0.10, else quality:core
when the core CC drops >= 0.03)."""
import os, re, sys, numpy as np
MODE, T, KIND, NODES, OUT, RP = sys.argv[1:7]; T = T.rstrip('/'); NODES = NODES.split(','); EXCL = set(sys.argv[7].split(',')) if len(sys.argv) > 7 else set()
K = '/mnt/beegfs/elmlund/afan/rnasep_klem'; TRUTH = {'cdisc': K + '/factorial/junk_test/simsites_truth_jt.txt', 'igm': K + '/factorial/armsigm1_truth/arms_truth.txt'}
def gap(c, dd, a, b):
    d, v = [], []
    for j in range(10):
        x, y = c[a & (dd == j)], c[b & (dd == j)]; x, y = x[np.isfinite(x)], y[np.isfinite(y)]
        if len(x) < 20 or len(y) < 20: continue
        d.append(np.median(x) - np.median(y)); v.append(1.5708 * (x.var() / len(x) + y.var() / len(y)))
    return (np.mean(d), np.sqrt(np.sum(v)) / len(d)) if d else (np.nan, np.nan)
def node(X):
    N = '%s/node_%s' % (T, X); R = RP + X; D = np.load(R + '/dump.npz'); rows = np.loadtxt(N + '/rows.txt', dtype=int); st = np.loadtxt(R + '/states.txt', dtype=int)
    idx = np.arange(len(D['corr'])) if len(D['corr']) == len(rows) else rows - 1; c, f, s = D['corr'][idx], D['dfx'][idx], st[idx]; ok = np.isfinite(f)
    dd = np.full(len(f), -1); dd[ok] = np.searchsorted(np.quantile(f[ok], np.linspace(0.1, 0.9, 9)), f[ok])
    m = re.search(r'dBx ([+-]?[0-9.]+) A\^2 \(SE ([0-9.]+)\)', open(R + '/dBref.txt').read()); g, gse = gap(c, dd, s == 1, s == 2)
    return dict(N=N, rows=rows, c=c, dd=dd, dBref=float(m.group(1)), dBref_se=float(m.group(2)), gref=g, gref_se=gse)
def cands(X, N):
    for k in sorted([k for k in os.listdir(N + '/motion') if re.match(r'c\d+$', k)], key=lambda k: int(k[1:])):
        M = '%s/motion/%s' % (N, k); f = [M + '/recvol_state%02d%s.mrc' % (b, x) for b in (1, 5) for x in ('', '_even', '_odd')]
        u = [M + '/recvol_state%02d_%s_unfil.mrc' % (b, x) for b in (1, 5) for x in ('even', 'odd')]
        if '%s:%s' % (X, k) not in EXCL and all(os.path.exists(p) for p in f + u + [M + '/bins.txt']): yield k, M, f
TC = {}
def fl(x):
    try: return float(x)
    except (TypeError, ValueError): return np.nan
def label(rows, b, L, H):
    if KIND not in TRUTH: return '-', '-'
    if not TC: TF = TRUTH[KIND]; TC['hd'] = [l.split()[1:] for l in open(TF) if l.startswith('#')][0]; TC['T'] = np.genfromtxt(TF, comments='#')
    hd, Tt = TC['hd'], TC['T']; rr = Tt[:, 0].astype(int) - 1
    def col(cn): v = np.full(rr.max() + 1, np.nan); v[rr] = Tt[:, hd.index(cn)]; return v[rows - 1]
    junk = col('junk') > 0; dj = junk[b == L].mean() - junk[b == H].mean(); J = abs(dj)
    if KIND == 'cdisc':
        jt = col('jtype'); fac = {'blob': junk & (jt == 1), 'broken': junk & (jt == 2), 'A': col('A_domain') > 0, 'B': col('B_blob') > 0, 'C': col('C_closed') > 0}
        d = {k: v[b == L].mean() - v[b == H].mean() for k, v in fac.items()}; sk = max('ABC', key=lambda k: abs(d[k])); S = abs(d[sk])
        lb = 'junk' if J >= 0.2 and J >= S else 'real:' + sk if S >= 0.2 else ('weak:' + ('junk' if J >= S else sk)) if max(J, S) >= 0.1 else 'none'
        return lb, 'djunk %+.2f ' % dj + ' '.join('d%s %+.2f' % (k, d[k]) for k in fac)
    eta = {}
    for a in [h for h in hd if re.match(r'a\d_(th|ph)$', h)]:
        v = col(a); s = ((b == L) | (b == H)) & ~junk & np.isfinite(v); x, gl = v[s], b[s] == L
        eta[a] = (gl.mean() * (x[gl].mean() - x.mean()) ** 2 + (~gl).mean() * (x[~gl].mean() - x.mean()) ** 2) / x.var() if x.var() > 0 else 0.
    em = max(eta, key=eta.get); E = eta[em]
    lb = 'mixed' if J >= 0.2 and E >= 0.1 else 'junk' if J >= 0.2 else 'real:' + em if E >= 0.1 else ('weak:' + ('junk' if J >= 0.1 else em)) if (J >= 0.1 or E >= 0.03) else 'none'
    return lb, 'djunk %+.2f eta2 %s %.3f' % (dj, em, E)
IH = 'name node n_node cand call L nL nH gap gap_se gref gref_se dBref dBref_se label detail'.split()
if MODE == 'cases':
    fc, fq, fi = open(OUT + '/cases.txt', 'w'), open(OUT + '/qscreen.txt', 'w'), open(OUT + '/cases_info.tsv', 'w'); fi.write('\t'.join(IH) + '\n')
    fd = open(OUT + '/dock_maplist.txt', 'w') if KIND == 'rnasep' else None; tg = os.path.basename(RP.rstrip('_'))
    for X in NODES:
        nd = node(X)
        for k, M, f in cands(X, nd['N']):
            b = np.loadtxt(M + '/bins.txt'); b = b[nd['rows'] - 1] if len(b) != len(nd['c']) and len(b) >= nd['rows'].max() else b  # per project row
            assert len(b) == len(nd['c']), '%s: bins %d, node rows %d' % (M, len(b), len(nd['c']))
            g, gse = gap(nd['c'], nd['dd'], b == 1, b == 5); L, H = (5, 1) if g > 0 else (1, 5); nm = '%s_%s_%s' % (tg, X, k); n = {1: (b == 1).sum(), 5: (b == 5).sum()}
            call = open(M + '/motion_call.txt').readline().strip() if os.path.exists(M + '/motion_call.txt') else '-'
            fc.write('%s cand %s %s - -\n' % (nm, ' '.join(f), nd['N'] + '/mask.mrc')); fq.write('QSCREEN %s cand | L = %s (%d) vs H (%d)\n' % (nm, 'B' if L == 5 else 'A', n[L], n[H]))
            lb, det = label(nd['rows'], b, L, H)
            fi.write('\t'.join(str(x) for x in (nm, X, len(b), k, call.replace(' ', '_') or '-', L, n[L], n[H], '%.5f' % abs(g), '%.5f' % gse, '%.5f' % nd['gref'], '%.5f' % nd['gref_se'],
                     nd['dBref'], nd['dBref_se'], lb, det)) + '\n')
            if fd: fd.write('%s_b1 other cand %d %s\n%s_b5 other cand %d %s\n' % (nm, n[1], ' '.join(f[:3]), nm, n[5], ' '.join(f[3:])))
            print('%s: L = bin %d, gap %.4f +- %.4f (reference %.4f), %s' % (nm, L, abs(g), gse, nd['gref'], lb))
else:
    QS = {}
    for l in open(OUT + '/qs3.txt'):
        m = re.match(r'QS3 (\S+)\s.*dBx ([+-]?[0-9.]+) A\^2 \(SE ([0-9.]+)\)', l)
        if m: QS[m.group(1)] = (float(m.group(2)), float(m.group(3)))
    FA, hd = {}, None
    if os.path.exists(OUT + '/dock/maps_table.txt'):
        for l in open(OUT + '/dock/maps_table.txt'):
            w = l.split()
            if l.startswith('#'): hd = w[1:]; continue
            if hd and w: FA[w[0]] = dict(zip(hd, w))
    rows = []
    for l in list(open(OUT + '/cases_info.tsv'))[1:]:
        r = dict(zip(IH, l.rstrip('\n').split('\t'))); nm = r['name']
        if nm not in QS: continue
        dB, se = QS[nm]; g, gse, gr, grse, dr, drse = (float(r[k]) for k in ('gap', 'gap_se', 'gref', 'gref_se', 'dBref', 'dBref_se')); ex = dr * g / gr
        ra = dB / ex if ex != 0 else np.nan; rse = np.sqrt((se / ex) ** 2 + ra ** 2 * ((drse / dr) ** 2 + (gse / g) ** 2 + (grse / gr) ** 2)) if ex != 0 else np.nan
        lb, det = r['label'], r['detail']
        if KIND == 'rnasep':
            L, H = r['L'], ('5' if r['L'] == '1' else '1'); a, h = FA.get('%s_b%s' % (nm, L)), FA.get('%s_b%s' % (nm, H))
            if a and h:
                d = {k: fl(a.get(k)) - fl(h.get(k)) for k in ('s1_A', 'closed_arm', 'core_cc_mean')}
                lb = 'real:arm' if abs(d['s1_A']) >= 3 or abs(d['closed_arm']) >= 0.10 else 'quality:core' if d['core_cc_mean'] <= -0.03 else 'none'
                det = 'ds1 %+.2f dclosed %+.3f dcore %+.3f' % (d['s1_A'], d['closed_arm'], d['core_cc_mean'])
        rows.append((nm, r['node'], int(r['n_node']), r['cand'], r['call'], int(r['L']), dB, se, dr, drse, g, gse, gr, ex, ra, rse, dB / dr, int(ra >= 1 and dB >= 3 * se), int(dB / dr >= 0.5 and dB >= 3 * se), lb, det))
    with open(OUT + '/ratios.txt', 'w') as fo:
        fo.write('# name node n_node cand call L dB dB_se dB_ref dB_ref_se gap gap_se ref_gap expected ratio_a ratio_a_se dB/dB_ref flag_a flag_half label detail\n')
        for x in rows: fo.write(' '.join(('%.4g' % v) if isinstance(v, float) else str(v) for v in x[:-1]) + ' | ' + x[-1] + '\n')
        fo.write('# SUMMARY ratio_a by label (flag_a = ratio_a >= 1 and dB >= 3 SE)\n')
        for lb in sorted(set(x[19].split(':')[0] for x in rows)):
            s = [x for x in rows if x[19].split(':')[0] == lb]
            fo.write('# %-8s %2d candidates: ratio_a %s; flagged %d\n' % (lb, len(s), ', '.join('%.2f+-%.2f' % (x[14], x[15]) for x in sorted(s, key=lambda x: -x[14])), sum(x[17] for x in s)))
    print(open(OUT + '/ratios.txt').read())
