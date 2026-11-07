#!/bin/bash
#SBATCH -J klem_joint_em_seeds
#SBATCH -p csbdevel
#SBATCH -c 24
#SBATCH --mem=200G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# seeding test for pole-seeded refinement on node r1: are seeds from ONE joint P-class soft EM (minor_em3, shared input poses, site_02 =
# 1f1 mask, started from each run's argmax pole seeds) more reproducible across runs than the argmax seeds themselves? P = 4 (z1-closed,
# open A, open B, z3-closed; argmax seeds disagree, port1 has 2.1k z1-closed rows) and P = 3 (closed merged) for script / port1 / port2.
# Writes rnasep_klem/repro_r1/jointem/<tag>/labels.txt (argmax of the joint responsibilities per master row) and compare.txt. env PARS
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=4
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K=$A/rnasep_klem; BASE=$K/factorial/f692st_r1_inposes
H=$T/simrnasep692; CP=$T/benchcache_rnasep692; MK=$BASE/picked/site_02.mrc; MEAN=$(ls -td $K/rflex_f692st/node_r1/mean/*_reconstruct3D | head -1)/recvol_state01.mrc
O=$K/repro_r1/jointem; mkdir -p $O; log(){ echo "[$(date '+%F %T') joint em seeds] $*"; }
run(){ local tag=$1 lab=$2 D=$O/$1; mkdir -p $D
  [ -s $D/start_B.txt ] || $PY - <<PYEOF || return 1
import numpy as np
use = np.loadtxt('$BASE/rows_use.txt', dtype=int) > 0; lab = np.loadtxt('$lab', dtype=int)
for h, mr in (('A', '${H}_A/meta_rows.txt'), ('B', '${H}_B/meta_rows.txt')):
    r = np.loadtxt(mr).astype(int) - 1; np.savetxt('$D/start_%s.txt' % h, np.where(use[r], lab[r], 0), fmt='%d')
PYEOF
  [ -s $D/em/post_soft_B.txt ] || $PY -u $T/minor_em3.py --sim ${H}_A,${H}_B --poses $BASE/poses_A.txt,$BASE/poses_B.txt --labels $D/start_A.txt,$D/start_B.txt \
    --cache $CP --ctfabs 1 --cons_ref $MEAN --mask $MK --soft 1 --iters ${ITERS:-12} --em_lim 15 --post 1 --nomaps 1 --out $D/em > $D/em.log 2>&1 || { log "$tag EM FAILED: $(tail -2 $D/em.log)"; return 1; }
  $PY - <<PYEOF
import numpy as np
use = np.loadtxt('$BASE/rows_use.txt', dtype=int) > 0; out = np.zeros(len(use), int); pm = np.zeros(len(use))
for h in ('A', 'B'):
    r = np.loadtxt('${H}_%s/meta_rows.txt' % h).astype(int) - 1; P = np.loadtxt('$D/em/post_soft_%s.txt' % h, ndmin=2)
    k = P.sum(1) > 0; out[r[k]] = 1 + P[k].argmax(1); pm[r[k]] = P[k].max(1)
out[~use] = 0; np.savetxt('$D/labels.txt', out, fmt='%d'); print('$tag: joint-EM states %s, mean max responsibility %.2f' % (np.bincount(out[use])[1:], pm[use].mean()))
PYEOF
  log "$tag done"; }
for p in 4 3; do for r in script port1 port2; do ( run p${p}_$r $K/repro_r1/poleseeds/poleseed${p}_$r.txt ) & while [ $(jobs -rp | wc -l) -ge ${PARS:-3} ]; do sleep 20; done; done; done; wait
$PY - <<PYEOF | tee $O/compare.txt
import numpy as np, itertools, os
use = np.loadtxt('$BASE/rows_use.txt', dtype=int) > 0; old = np.loadtxt('$K/f692r/coordtree_fam1/units.txt', dtype=int); old = ((old == 15) | (old == 16))
def nmi(a, b):
    ka, kb = np.unique(a), np.unique(b); n = len(a); M = np.array([[np.sum((a == i) & (b == j)) for j in kb] for i in ka], float) / n
    pa, pb = M.sum(1), M.sum(0); nz = M > 0; I = (M[nz] * np.log(M[nz] / np.outer(pa, pb)[nz])).sum()
    ha, hb = -(pa[pa > 0] * np.log(pa[pa > 0])).sum(), -(pb[pb > 0] * np.log(pb[pb > 0])).sum(); return 2 * I / (ha + hb)
for p in (4, 3):
    for kind in ('argmax seeds', 'joint EM'):
        L = {}
        for r in ('script', 'port1', 'port2'):
            f = '$K/repro_r1/poleseeds/poleseed%d_%s.txt' % (p, r) if kind == 'argmax seeds' else '$O/p%d_%s/labels.txt' % (p, r)
            if os.path.exists(f): L[r] = np.loadtxt(f, dtype=int)[use]
        print('==== P = %d, %s: sizes %s' % (p, kind, ' | '.join('%s %s' % (r, np.bincount(l, minlength=p + 1)[1:]) for r, l in L.items())))
        for a, b in itertools.combinations(L, 2):
            J = [len(np.flatnonzero((L[a] == k) & (L[b] == k))) / max(1, len(np.flatnonzero((L[a] == k) | (L[b] == k)))) for k in range(1, p + 1)]
            print('  %s vs %s: NMI %.2f | per-state Jaccard %s' % (a, b, nmi(L[a], L[b]), ' '.join('%.2f' % j for j in J)))
        for r, l in L.items():
            base = old[use].mean(); print('  %s old-closed enrichment per state: %s' % (r, ' '.join('%.1fx' % (old[use][l == k].mean() / base) for k in range(1, p + 1))))
PYEOF
log DONE
