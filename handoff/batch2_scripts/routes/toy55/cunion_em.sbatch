#!/bin/bash
#SBATCH -J klem_cunion
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# port fulltoy_cdisc: does merging the port's partial site-C candidates recover the closed state? Root site EM (30 it) from SEED's pole with the
# union of the C-like candidate masks vs the true C mask (ceiling), same start; closed precision / recall / AUC / top-N vs truth. env O (port_sites_cdisc) PT SEED
set -uo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; O=${O:-$K_/factorial/port_sites_cdisc}; PT=${PT:-$K_/stree_port/fulltoy_cdisc/1_hetero_tree}
D=$O/cunion; mkdir -p $D; H=$T/simsites; SEED=${SEED:-$O/fits/r1_s4@r/hi/top_star.txt}; MEAN=$O/view/node_r/mean/1_reconstruct3D/recvol_state01.mrc
$PY - <<PYEOF || exit 1
import numpy as np
def rd(p): h = np.fromfile(p, np.int32, 256); n = int(h[0]); return h, np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n)
ms = ['$PT/node_r/sites/site_02.mrc', '$PT/node_r1/sites/site_03.mrc', '$PT/node_r1/sites/site_04.mrc', '$PT/node_r1/sites/site_05.mrc']
h, U = rd(ms[0])
for p in ms[1:]: U = np.maximum(U, rd(p)[1])
hh = h.copy(); hh[23] = 0; hf = hh.view(np.float32); hf[19:22] = U.min(), U.max(), U.mean()
with open('$D/union_C.mrc', 'wb') as f: f.write(hh.tobytes()); f.write(U.astype(np.float32).tobytes())
C = rd('$T/maps160/simsites/mask_C.mrc')[1] > 0.5; u = U > 0.5; print('union of the 4 C-like candidates: %d voxels > 0.5, covers %.2f of true C, Dice %.2f' % (u.sum(), (u & C).sum() / C.sum(), 2. * (u & C).sum() / (u.sum() + C.sum())))
s = np.loadtxt('$SEED', dtype=int) > 0
for h_ in 'AB': r = np.loadtxt('${H}_%s/meta_rows.txt' % h_).astype(int) - 1; np.savetxt('$D/start_%s.txt' % h_, 1 + s[r].astype(int), fmt='%d')
print('seed: %d particles' % s.sum())
PYEOF
for m in union:$D/union_C.mrc true:$T/maps160/simsites/mask_C.mrc; do nm=${m%%:*}; mk=${m#*:}
  $PY -u $T/minor_em3.py --sim ${H}_A,${H}_B --poses $O/node_r/poses_A.txt,$O/node_r/poses_B.txt --labels $D/start_A.txt,$D/start_B.txt --cache $T/benchcache_sites --ctfabs 0 --cons_ref $MEAN \
    --mask $mk --soft 1 --iters ${ITERS:-30} --em_lim 15 --post 1 --nomaps 1 --out $D/em_$nm 2>&1 | grep -v Warning > $D/em_$nm.log & done; wait
for nm in union true; do grep -a "^iteration" $D/em_$nm.log | tail -2 | cut -c1-160; done
$PY - <<PYEOF
import numpy as np
T = np.loadtxt('$K_/factorial/junk_test/simsites_truth_jt.txt'); hd = [l.split()[1:] for l in open('$K_/factorial/junk_test/simsites_truth_jt.txt') if l.startswith('#')][0]
c = lambda n: T[:, hd.index(n)]; junk = c('junk') > 0; y = (c('C_closed') > 0) & ~junk; g = ~junk; NM = len(T)
def auc(s, yy): r_ = np.argsort(np.argsort(s)) + 1.; n1 = yy.sum(); n0 = len(yy) - n1; return (r_[yy].sum() - n1 * (n1 + 1) / 2) / max(n1 * n0, 1)
for nm in ('union', 'true'):
    P = np.full(NM, np.nan)
    for h_ in 'AB': r = np.loadtxt('${H}_%s/meta_rows.txt' % h_).astype(int) - 1; P[r] = np.loadtxt('$D/em_%s/post_soft_%s.txt' % (nm, h_))[:, 1]
    cl = P >= 0.5; o = np.argsort(-np.where(g, P, -np.inf))
    print('%s mask: class 2 = %d (junk %.3f); closed precision %.3f recall %.3f (base %.3f); closed AUC %.3f; top-2.5k purity %.3f, top-10k %.3f'
          % (nm, cl.sum(), (cl & junk).sum() / max(cl.sum(), 1), (cl & y).sum() / max((cl & g).sum(), 1), (cl & y).sum() / y.sum(), y.sum() / g.sum(), auc(P[g], y[g]), y[o[:2500]].mean(), y[o[:10000]].mean()))
PYEOF
