#!/bin/bash
#SBATCH -J klem_treewiener_real
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# real RNase P (f692r families, merge-up states): (1) tree-weighted maps without truth (tree_wiener3: covariance fit of the tree vs a
# shuffled tree, borrowing readout; raw / own / tree2 / emp2 leaf maps) -> closed-arm readout of each map set per state (does borrowing
# keep the rare closed states' arm?); (2) per-node focus masks (node_masks2) and their overlap with the closed-arm / open-arm regions.
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem; R=$K_/f692r
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16; CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
for F in ${FAMS:-1 2}; do
  C=$R/coordtree_fam$F; RR=$(ls -td $C/rec/*_reconstruct3D | head -1); O=$K_/treewiener_real/fam$F; echo "== f692r family $F ($(awk '$1 > 0' $C/units.txt | sort -n | uniq | wc -l) states)"
  [ -s $O/summary.txt ] || $PY $T/tree_wiener3.py $O --rdir $RR --labels $C/units.txt --tree $C/tree/merge_tree.txt --mskdiam 232 --maps raw,own,tree2,emp2 2>&1 | grep -v Warning
  for k in raw own tree2 emp2; do [ -s $O/arm_$k.txt ] || $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/arm_$k.txt $O/maps/$k > /dev/null 2>&1; done
  $PY - <<PYEOF
import numpy as np, re
nm = {int(l.split()[0]): int(l.split()[1]) for l in open('$C/names.txt') if l.strip()}; u = np.loadtxt('$C/units.txt', dtype=int)
A = {}
for k in ('raw', 'own', 'tree2', 'emp2'):
    for l in open('$O/arm_%s.txt' % k):
        if l.startswith('#') or not l.strip(): continue
        w = l.split(); s = int(re.findall(r'state_(\d+)', w[0])[0]); A.setdefault(s, {})[k] = (float(w[1]), float(w[2]))
rows = sorted(A, key=lambda s: -A[s]['raw'][0])
print('closed-arm occupancy per state (raw / own filter / tree2 / emp2), open-arm in brackets; top 8 by raw closed-arm and the 4 smallest:')
sel = rows[:8] + [s for s in sorted(A, key=lambda s: (u == s).sum())[:4] if s not in rows[:8]]
for s in sel: print('   state %2d (orig %3d, %6d ptcls): %s' % (s, nm.get(s, s), (u == s).sum(), '  '.join('%.2f(%.2f)' % A[s][k] for k in ('raw', 'own', 'tree2', 'emp2'))))
c = {k: np.array([A[s][k][0] for s in rows]) for k in ('raw', 'own', 'tree2', 'emp2')}
print('corr of closed-arm occupancy across states with raw: own %.3f tree2 %.3f emp2 %.3f; spread (sd) raw %.3f own %.3f tree2 %.3f emp2 %.3f' % tuple([np.corrcoef(c['raw'], c[k])[0, 1] for k in ('own', 'tree2', 'emp2')] + [c[k].std() for k in ('raw', 'own', 'tree2', 'emp2')]))
PYEOF
  NM=$K_/nodemasks_real/fam$F; [ -s $NM/node_masks.txt ] || $PY $T/node_masks2.py $NM --rdir $RR --labels $C/units.txt --tree $C/tree/merge_tree.txt --mskdiam 232 > /dev/null 2>&1
  $PY - <<PYEOF
import numpy as np, glob, os, sys
sys.path.insert(0, '$T'); from model_class import read_mrc
smpd = 1.464; LP = 10.; cons = read_mrc('$T/maps160/toy_consensus.mrc'); n = cons.shape[0]; f = np.fft.fftfreq(n, d=smpd); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); G = np.exp(-0.5 * (kx**2 + ky**2 + kz**2) * (LP / 1.2)**2)
lp = lambda v: np.real(np.fft.ifftn(np.fft.fftn(v.astype(np.float32)) * G)); cons = lp(cons); supp = cons > 0.2 * cons.max(); outside = cons < 0.02 * cons.max()
def std(v): v = v - v[outside].mean(); return v / max(np.sqrt(np.mean(v[supp] ** 2)), 1e-12)
cs = std(cons); cl = std(lp(read_mrc('$CL'))); d = cl - cs; thr = 0.5 * np.percentile(np.abs(d[supp | (cl > 0.2 * cl.max())]), 99); mC = d > thr; mO = d < -thr
rows = {int(l.split()[0]): l for l in open('$NM/node_masks.txt') if l.strip() and not l.startswith('#')}
print('node masks (sig80) vs the arm regions: node leaves particles | mask %% of the box support | closed-arm voxels covered | open-arm covered | share of the mask on the arm regions')
for u in sorted(rows, key=lambda r: -int(rows[r].split()[1])):
    fm = '$NM/masks/node%03d_sig80.mrc' % u
    if not os.path.isfile(fm): continue
    m = read_mrc(fm) > 0.5; w = rows[u].split()
    print('   node %3d %2s %7s | %5.1f %% | %.2f | %.2f | %.2f' % (u, w[1], w[2], 100 * (m & supp).sum() / supp.sum(), (m & mC).sum() / max(mC.sum(), 1), (m & mO).sum() / max(mO.sum(), 1), (m & (mC | mO)).sum() / max(m.sum(), 1)))
PYEOF
done
echo DONE
