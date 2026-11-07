#!/bin/bash
#SBATCH -J klem_site_reseed
#SBATCH -p csbdevel
#SBATCH -c 24
#SBATCH --mem=160G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# site EM of BASE's site J re-seeded per SEEDS 'name=labelfile@id+id | name=random@N | name=flex@N' -> ranking agreement with BASE and top-N maps; env RT NODE SRC FLEXDIR FLEXARM FLEXNAME TOPN PARS MSKDIAM
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=4
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; RT=${RT:?}; R=$K_/rflex_$RT; SRC=${SRC:?}; ND=${NODE:-r}
H=${H:-$T/simrnasep692}; MRA=${MRA:-${H}_A/meta_rows.txt}; MRB=${MRB:-${H}_B/meta_rows.txt}; CP=${CP:-$T/benchcache_rnasep692}; CTFABS=${CTFABS:-1}; MD=${MSKDIAM:-220}
BASE=${BASE:?}; J=${J:-2}; O=${OUT:-$BASE/reseed}; V=$K_/view_rflex/reseed_$(basename $BASE); mkdir -p $O $V; cd $O; log(){ echo "[$(date '+%F %T') site reseed $RT] $*"; }
SD=$BASE/picked; [ -d $SD ] || SD=$BASE/discovery; MK=$SD/site_$(printf %02d $J).mrc; NM=$(grep -c . $BASE/rows_use.txt); MEAN=$(ls -td $R/node_$ND/mean/*_reconstruct3D | head -1)/recvol_state01.mrc
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1); CONS=$T/maps160/toy_consensus.mrc; OLD=${OLD:-$K_/f692r/coordtree_fam1/units.txt}
flexrank(){ [ -s $O/flexrank.txt ] && return 0; until [ -s ${FLEXARM:?} ]; do sleep 30; done; $PY - <<PYEOF
import numpy as np, glob, os
fd = sorted(glob.glob('${FLEXDIR:?}/*_flex_pca'), key=os.path.getmtime)[-1]; C = np.loadtxt(fd + '/flex_pca_coordinates.txt', comments='#'); r, lab, Z = C[:, 0].astype(int) - 1, C[:, 2].astype(int), C[:, 5:]
arm = {os.path.basename(w[0])[:-4]: float(w[1]) for w in (l.split() for l in open('$FLEXARM')) if w and not w[0].startswith('#')}
s = max((k for k in arm if k.startswith('${FLEXNAME:-closedsite}_state')), key=lambda k: arm[k]); sid = int(s.split('_state')[1]); u = Z[lab == sid].mean(0) - Z.mean(0)
sc = np.full($NM, -1e30); sc[r] = (Z - Z.mean(0)) @ u; np.savetxt('$O/flexrank.txt', sc, fmt='%.6g'); print('flex ranking: closed-most state %d of %s (closed-arm %.2f)' % (sid, fd, arm[s]))
PYEOF
}
seed(){ local nm=$1 spec=$2 D=$O/$1; mkdir -p $D; [ "${spec%@*}" = flex ] && flexrank
  [ -s $D/start_B.txt ] || $PY - <<PYEOF || return 1
import numpy as np
use = np.loadtxt('$BASE/rows_use.txt', dtype=int) > 0; kind, arg = '$spec'.rsplit('@', 1); seed = np.zeros(len(use), bool)
if kind == 'random': seed[np.random.default_rng(7).choice(np.flatnonzero(use), int(arg), replace=False)] = True
elif kind == 'flex': sc = np.loadtxt('$O/flexrank.txt'); sc[~use] = -1e30; seed[np.argsort(-sc)[:int(arg)]] = True
else: seed = np.isin(np.loadtxt(kind, dtype=int), [int(x) for x in arg.split('+')])
seed &= use; np.savetxt('$D/seed.txt', seed.astype(int), fmt='%d')
for h, mr in (('A', '$MRA'), ('B', '$MRB')):
    r = np.loadtxt(mr).astype(int) - 1; np.savetxt('$D/start_%s.txt' % h, np.where(use[r], 1 + seed[r].astype(int), 0), fmt='%d')
print('seed $nm ($spec): %d rows' % seed.sum())
PYEOF
  [ -s $D/em/post_soft_B.txt ] || $PY -u $T/minor_em3.py --sim ${H}_A,${H}_B --poses $BASE/poses_A.txt,$BASE/poses_B.txt --labels $D/start_A.txt,$D/start_B.txt --cache $CP --ctfabs $CTFABS \
    --cons_ref $MEAN --mask $MK --soft 1 --iters ${ITERS:-10} --em_lim ${EMLIM:-15} --post 1 --nomaps 1 --out $D/em > $D/em.log 2>&1 || { log "seed $nm EM FAILED: $(tail -2 $D/em.log)"; return 1; }
  [ -s $D/post_master.txt ] || $PY -u $T/site_size.py --em $D/em --start $D/start_A.txt,$D/start_B.txt --sim ${H}_A,${H}_B --poses $BASE/poses_A.txt,$BASE/poses_B.txt --mr $MRA,$MRB --nmaster $NM \
    --mask $MK --out $D --cache $CP --ctfabs $CTFABS --em_lim ${EMLIM:-15} > $D/size.log 2>&1 || { log "seed $nm sizing FAILED: $(tail -2 $D/size.log)"; return 1; }
  log "seed $nm: $(grep -a '^N\*' $D/size.log | cut -c1-120)"; }
log "site $J of $BASE ($MK) re-seeded: ${SEEDS:?}"
for e in $SEEDS; do ( seed ${e%%=*} ${e#*=} ) & while [ $(jobs -rp | wc -l) -ge ${PARS:-2} ]; do sleep 20; done; done; wait
RK="base=$BASE/site_$J/post_master.txt"; for e in $SEEDS; do RK="$RK ${e%%=*}=$O/${e%%=*}/post_master.txt"; done; [ -s $O/flexrank.txt ] && RK="$RK flexrank=$O/flexrank.txt"
$PY - <<PYEOF | tee $O/agreement.txt
import numpy as np, os
use = np.loadtxt('$BASE/rows_use.txt', dtype=int) > 0; old = np.loadtxt('$OLD', dtype=int); old = ((old == 15) | (old == 16)) & use; NS = [int(x) for x in '${TOPN:-1000,2500,5000,10000,15000}'.split(',')]
rk = {}
for e in '$RK'.split():
    nm, f = e.split('=', 1)
    if os.path.exists(f): p = np.loadtxt(f); p[~use] = -1e30; rk[nm] = p
top = {nm: {N: set(np.argsort(-p)[:N]) for N in NS} for nm, p in rk.items()}; rank = {nm: np.argsort(np.argsort(-p[use])) for nm, p in rk.items()}
print('rows %d, old f692r closed (15+16) %d (%.1f %%); top-N sets: Jaccard per pair at N %s; Spearman over all rows' % (use.sum(), old.sum(), 100. * old.sum() / use.sum(), NS))
nms = list(rk)
for i in range(len(nms)):
    for j in range(i + 1, len(nms)):
        a, b = nms[i], nms[j]; print('%-10s vs %-10s Spearman %+.3f | Jaccard %s' % (a, b, np.corrcoef(rank[a], rank[b])[0, 1], ' '.join('%.2f' % (len(top[a][N] & top[b][N]) / len(top[a][N] | top[b][N])) for N in NS)))
print('old-closed enrichment of the top-N (x base rate) and share of all old closed captured:')
for nm in nms: print('%-10s %s' % (nm, ' '.join('%d: %.2fx %.0f%%' % (N, old[list(top[nm][N])].mean() / old[use].mean(), 100. * old[list(top[nm][N])].sum() / old.sum()) for N in NS)))
for nm in nms:
    f = '$O/%s/seed.txt' % nm
    if os.path.exists(f): s = np.loadtxt(f, dtype=int) > 0; n = s.sum(); t = set(np.argsort(-rk[nm])[:n]); print('%-10s seed %d rows: %.2f of them in its own final top-%d' % (nm, n, len(t & set(np.flatnonzero(s))) / n, n))
PYEOF
rec(){ local d=$1 lab=$2; mkdir -p $d; [ -n "$(ls $d/*_reconstruct3D/recvol_state01_even_unfil.mrc 2>/dev/null)" ] && return 0; cd $d && cp -f $SRC r.simple
  simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$lab mkdir=no > sel.log 2>&1 && simple_exec prg=update_project projfile=r.simple qsys_name=local nparts=3 nthr=8 mkdir=no > up.log 2>&1 \
    && simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=$MD nstates=1 nparts=3 nthr=8 > rec.log 2>&1 || log "reconstruct $d FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-150)"; cd $O; }
: > $O/spec_top.txt
for e in $(echo $RK | tr ' ' '\n' | grep -v '^base='); do nm=${e%%=*}; for N in ${MAPN:-2500 5000 10000}; do
  [ -s $O/top_${nm}_$N.txt ] || $PY -c "import numpy as np; p = np.loadtxt('${e#*=}'); u = np.loadtxt('$BASE/rows_use.txt') > 0; p[~u] = -1e30; l = np.zeros(len(p), int); l[np.argsort(-p)[:$N]] = 1; np.savetxt('$O/top_${nm}_$N.txt', l, fmt='%d')"
  rec $O/rec_${nm}_top$N $O/top_${nm}_$N.txt; RD=$(ls -td $O/rec_${nm}_top$N/*_reconstruct3D 2>/dev/null | head -1); [ -n "$RD" ] && echo "${nm}_top$N $RD/recvol_state01_even_unfil.mrc $RD/recvol_state01_odd_unfil.mrc" >> $O/spec_top.txt; done; done
$PY $T/make_view.py $V $O/spec_top.txt --box 160 > /dev/null 2>&1; $PY $T/arm_readout.py $CONS $CL $O/arm_top.txt $(awk -v v=$V '{print v "/" $1 ".mrc"}' $O/spec_top.txt) 2>&1 | grep -v Warning | tail -1
grep -v '^#' $O/arm_top.txt | awk '{n = $1; sub(/.*\//, "", n); printf "%s closed-arm %s open-arm %s %s\n", n, $2, $3, $4}'
log "maps (lp 10, unsharpened, full maps only): $V"; log DONE
