#!/bin/bash
#SBATCH -J klem_cunion2
#SBATCH -p csbdevel
#SBATCH -c 24
#SBATCH --mem=96G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# port fulltoy_cdisc: root site EMs (30 it) on data-driven unions of the port's site candidates; SETS 'name:maskA,maskB,..:seed_top_star ..'
# (masks relative to PT) -> closed precision / recall / AUC / top-N vs truth. env O (port_sites_cdisc) PT SETS ITERS
set -uo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; O=${O:-$K_/factorial/port_sites_cdisc}; PT=${PT:-$K_/stree_port/fulltoy_cdisc/1_hetero_tree}
D=$O/cunion2; mkdir -p $D; H=$T/simsites; MEAN=$O/view/node_r/mean/1_reconstruct3D/recvol_state01.mrc
for e in ${SETS:?}; do IFS=: read nm ms sd <<< "$e"; mkdir -p $D/$nm
  $PY - <<PYEOF || continue
import numpy as np
def rd(p): h = np.fromfile(p, np.int32, 256); n = int(h[0]); return h, np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n)
ms = ['$PT/' + m for m in '$ms'.split(',')]; h, U = rd(ms[0])
for p in ms[1:]: U = np.maximum(U, rd(p)[1])
hh = h.copy(); hh[23] = 0; hf = hh.view(np.float32); hf[19:22] = U.min(), U.max(), U.mean()
with open('$D/$nm/mask.mrc', 'wb') as f: f.write(hh.tobytes()); f.write(U.astype(np.float32).tobytes())
C = rd('$T/maps160/simsites/mask_C.mrc')[1] > 0.5; u = U > 0.5; s = np.loadtxt('$sd', dtype=int) > 0
print('$nm: union of %s: covers %.2f of true C (Dice %.2f); seed %d particles' % ('$ms', (u & C).sum() / C.sum(), 2. * (u & C).sum() / (u.sum() + C.sum()), s.sum()))
for h_ in 'AB': r = np.loadtxt('${H}_%s/meta_rows.txt' % h_).astype(int) - 1; np.savetxt('$D/$nm/start_%s.txt' % h_, 1 + s[r].astype(int), fmt='%d')
PYEOF
  $PY -u $T/minor_em3.py --sim ${H}_A,${H}_B --poses $O/node_r/poses_A.txt,$O/node_r/poses_B.txt --labels $D/$nm/start_A.txt,$D/$nm/start_B.txt --cache $T/benchcache_sites --ctfabs 0 --cons_ref $MEAN \
    --mask $D/$nm/mask.mrc --soft 1 --iters ${ITERS:-30} --em_lim 15 --post 1 --nomaps 1 --out $D/$nm/em 2>&1 | grep -v Warning > $D/$nm/em.log & done; wait
$PY - <<PYEOF
import numpy as np, os
T = np.loadtxt('$K_/factorial/junk_test/simsites_truth_jt.txt'); hd = [l.split()[1:] for l in open('$K_/factorial/junk_test/simsites_truth_jt.txt') if l.startswith('#')][0]
c = lambda n: T[:, hd.index(n)]; junk = c('junk') > 0; y = (c('C_closed') > 0) & ~junk; g = ~junk; NM = len(T)
def auc(s, yy): r_ = np.argsort(np.argsort(s)) + 1.; n1 = yy.sum(); n0 = len(yy) - n1; return (r_[yy].sum() - n1 * (n1 + 1) / 2) / max(n1 * n0, 1)
for e in '$SETS'.split():
    nm = e.split(':')[0]; P = np.full(NM, np.nan)
    if not os.path.exists('$D/%s/em/post_soft_B.txt' % nm): print(nm, 'EM FAILED'); continue
    for h_ in 'AB': r = np.loadtxt('${H}_%s/meta_rows.txt' % h_).astype(int) - 1; P[r] = np.loadtxt('$D/%s/em/post_soft_%s.txt' % (nm, h_))[:, 1]
    cl = P >= 0.5; o = np.argsort(-np.where(g, P, -np.inf))
    print('%s: class 2 = %d (junk %.3f); closed precision %.3f recall %.3f (base %.3f); closed AUC %.3f; top-2.5k purity %.3f, top-10k %.3f; %s'
          % (nm, cl.sum(), (cl & junk).sum() / max(cl.sum(), 1), (cl & y).sum() / max((cl & g).sum(), 1), (cl & y).sum() / y.sum(), y.sum() / g.sum(), auc(P[g], y[g]), y[o[:2500]].mean(), y[o[:10000]].mean(),
             open('$D/%s/em.log' % nm).read().strip().split('\n')[-2][:70]))
PYEOF
