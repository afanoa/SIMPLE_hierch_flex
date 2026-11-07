#!/bin/bash
#SBATCH -J klem_nodemasks3
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# split-direction focus masks (node_masks3): panel sets against truth, then f692r families against the closed / open arm regions
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16; CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
for t in swing30hs prodihs; do P=$K_/panel_$t; echo "== $t"
  $PY $T/node_masks3.py $K_/nodemasks3/$t --rdir $P/1_refine3D_states --labels $P/refined_states.txt --tree $P/smerge/merge_tree.txt --gt $T/sim$t/gt_particles.txt --truth $T/sim$t/maps --mskdiam 220 2>&1 | grep -v Warning | grep -o "^[0-9]* [0-9]* [0-9]* | [0-9.]*\|true support [0-9.]* |\|sig80 [0-9.]*% rec [0-9.]*\|split [^|]*" | paste - - - - | sed 's/^/   /'
done
for F in 1 2; do C=$K_/f692r/coordtree_fam$F; RR=$(ls -td $C/rec/*_reconstruct3D | head -1); NM=$K_/nodemasks3/f692r_fam$F; echo "== f692r family $F"
  [ -s $NM/node_masks.txt ] || $PY $T/node_masks3.py $NM --rdir $RR --labels $C/units.txt --tree $C/tree/merge_tree.txt --mskdiam 232 > /dev/null 2>&1
  $PY - <<PYEOF
import numpy as np, os, sys
sys.path.insert(0, '$T'); from model_class import read_mrc
smpd = 1.464; LP = 10.; cons = read_mrc('$T/maps160/toy_consensus.mrc'); n = cons.shape[0]; f = np.fft.fftfreq(n, d=smpd); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); G = np.exp(-0.5 * (kx**2 + ky**2 + kz**2) * (LP / 1.2)**2)
lp = lambda v: np.real(np.fft.ifftn(np.fft.fftn(v.astype(np.float32)) * G)); cons = lp(cons); supp = cons > 0.2 * cons.max(); outside = cons < 0.02 * cons.max()
def std(v): v = v - v[outside].mean(); return v / max(np.sqrt(np.mean(v[supp] ** 2)), 1e-12)
cs = std(cons); cl = std(lp(read_mrc('$CL'))); d = cl - cs; thr = 0.5 * np.percentile(np.abs(d[supp | (cl > 0.2 * cl.max())]), 99); mC = d > thr; mO = d < -thr; arm = mC | mO
print('arm regions: closed %d + open %d voxels = %.1f %% of the support' % (mC.sum(), mO.sum(), 100 * (arm & supp).sum() / supp.sum()))
rows = {int(l.split()[0]): l.split() for l in open('$NM/node_masks.txt') if l.strip() and not l.startswith('#')}
print('   node leaves particles | members mask: %% support, closed / open arm covered, share on arm | split mask: same')
for u in sorted(rows, key=lambda r: -int(rows[r][1])):
    out = []
    for sfx in ('_sig80', '_split'):
        fm = '$NM/masks/node%03d%s.mrc' % (u, sfx)
        if not os.path.isfile(fm): out.append('-'); continue
        m = read_mrc(fm) > 0.5; out.append('%5.1f %.2f/%.2f %.2f' % (100 * (m & supp).sum() / supp.sum(), (m & mC).sum() / max(mC.sum(), 1), (m & mO).sum() / max(mO.sum(), 1), (m & arm).sum() / max(m.sum(), 1)))
    if int(rows[u][1]) >= 3: print('   node %3d %2s %7s | %s | %s' % (u, rows[u][1], rows[u][2], out[0], out[1]))
PYEOF
done
echo DONE
