#!/bin/bash
#SBATCH -p norm
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:15:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/scorecard/%x.%j.out
# The m33 root sites under the 2D-corr-ranked quality reference (ptcl2D corr: no 3D poses): frame agreement of dB_ref (stage B vs
# recB0), each site's lower-2D-corr end bin, dB / dB_ref and option (a) with 2D-corr gaps. Site dB from quality_a/m79 (recB0 bin maps).
/mnt/nasapps/production/python/3.12/bin/python3 - <<'PYEOF'
import re, numpy as np
Q = '/mnt/beegfs/elmlund/afan/rnasep_klem/scorecard/'; G = '/mnt/beegfs/elmlund/afan/rnasep_klem/stree_port/gate/m33_root_m79/1_hetero_tree/node_r'
def ref(t):
    m = re.search(r'dBx ([+-]?[0-9.]+) A\^2 \(SE ([0-9.]+)\)', open(Q + 'quality_ref/%s/dBref.txt' % t).read()); return float(m.group(1)), float(m.group(2))
rb, rr = ref('m79_root_2d_bframe'), ref('m79_root_2d_recB0'); sb = np.loadtxt(Q + 'quality_ref/m79_root_2d_bframe/states.txt', dtype=int); sr = np.loadtxt(Q + 'quality_ref/m79_root_2d_recB0/states.txt', dtype=int)
print('2D-ranked reference: stage-B frame %+.1f +- %.1f, recB0 frame %+.1f +- %.1f (difference %.1f SE); quintile states identical in both: %s' % (rb + rr + (abs(rb[0] - rr[0]) / np.hypot(rb[1], rr[1]), bool((sb == sr).all()))))
for t in ('m79_root_2d_bframe', 'm79_root_2d_recB0'): print('  %s: %s' % (t, open(Q + 'quality_ref/%s/groups.log' % t).read().strip()))
D = np.load(Q + 'quality_ref/m79_root_2d_recB0/dump.npz'); c2 = np.load(Q + 'quality_ref/m79_root_2d_recB0/dump2d.npz')['corr']; c3 = D['corr']; f = D['dfx']; ok = np.isfinite(f)
dd = np.full(len(f), -1); dd[ok] = np.searchsorted(np.quantile(f[ok], np.linspace(0.1, 0.9, 9)), f[ok])
print('  2D vs recB0 3D corr: r %.3f' % np.corrcoef(c2[np.isfinite(c2) & np.isfinite(c3)], c3[np.isfinite(c2) & np.isfinite(c3)])[0, 1])
def gap(c, a, b):
    d = [np.median(c[a & (dd == j)]) - np.median(c[b & (dd == j)]) for j in range(10) if (a & (dd == j)).sum() > 20 and (b & (dd == j)).sum() > 20]; return np.mean(d)
gr = gap(c2, sr == 1, sr == 2)
I = {l.split('\t')[0]: l.rstrip('\n').split('\t') for l in list(open(Q + 'quality_a/m79/cases_info.tsv'))[1:]}; S = {}
for l in open(Q + 'quality_a/m79/qs3.txt'):
    m = re.match(r'QS3 (\S+)\s.*dBx ([+-]?[0-9.]+) A\^2 \(SE ([0-9.]+)\)', l)
    if m: S[m.group(1)] = (float(m.group(2)), float(m.group(3)))
LBL = {}
for l in open(Q + 'quality_a/m79/ratios.txt'):
    if not l.startswith('#'): w = l.split(' | ')[0].split(); LBL[w[0]] = w[19] + ' ' + l.split(' | ')[1].strip()
print('site  L(3D) L(2D)  2D gap   dB (L by 2D)     dB/dB_ref   ratio_a(2D gaps)   label')
for s in (1, 2, 3, 5, 6, 7):
    nm = 'qa_m79_r_c%d' % s; b = np.loadtxt(G + '/motion/c%d/bins.txt' % s); g15 = gap(c2, b == 1, b == 5); L2 = 5 if g15 > 0 else 1; L3 = int(I[nm][5])
    dB, se = S[nm]; dB2 = dB if L2 == L3 else -dB; ex = rr[0] * abs(g15) / gr
    print('%d     bin%d  bin%d  %.5f  %+7.1f +- %4.1f   %+.3f       %+7.2f            %s' % (s, L3, L2, abs(g15), dB2, se, dB2 / rr[0], dB2 / ex, LBL.get(nm, '')))
print('reference 2D-corr quintile gap %.5f' % gr)
PYEOF
