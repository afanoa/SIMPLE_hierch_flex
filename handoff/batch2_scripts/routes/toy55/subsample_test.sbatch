#!/bin/bash
#SBATCH -J klem_subsample
#SBATCH -p csbdevel
#SBATCH -c 32
#SBATCH --mem=128G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# site EM fitted on a random FRACS share of the node (same seed, same mask / iterations; minor_em3s.py) and every particle scored once,
# against the existing full-node fit: Spearman (all / held-out particles), top-N overlap, class-2 share, truth AUC (toys). env FRACS PAR
set -uo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K=$A/rnasep_klem; F=$K/factorial; O=$F/subsample_test; mkdir -p $O
P=$F/port_sites_cdisc; R1=$F/f692st_treesites_r1b; HS=$T/simsites; HR=$T/simrnasep692; MEANR1=$(ls -td $K/rflex_f692st/node_r1/mean/*_reconstruct3D | head -1)/recvol_state01.mrc
# name | sim | cache | ctfabs | poses prefix | mask | cons_ref | start dir | iters | full em dir | truth column
cat > $O/cases.txt <<CASES
toyA $HS $T/benchcache_sites 0 $P/node_r/poses $P/node_r/disc/site_01.mrc $P/view/node_r/mean/1_reconstruct3D/recvol_state01.mrc $P/fits/r_s1@r/lo 10 $P/fits/r_s1@r/lo/em A_domain
toyC_refined $HS $T/benchcache_sites 0 $P/node_r/poses $P/maskref/round1/mask.mrc $P/view/node_r/mean/1_reconstruct3D/recvol_state01.mrc $P/maskref 30 $P/maskref/round1/em C_closed
toyC_true $HS $T/benchcache_sites 0 $P/node_r/poses $T/maps160/simsites/mask_C.mrc $P/view/node_r/mean/1_reconstruct3D/recvol_state01.mrc $P/cunion 30 $P/cunion/em_true C_closed
r1_1f1_lo $HR $T/benchcache_rnasep692 1 $R1/node_r1/poses $R1/node_r1/disc/site_01.mrc $MEANR1 $R1/fits/r1_s1@r1/lo 10 $R1/fits/r1_s1@r1/lo/em -
r1_1f1_hi $HR $T/benchcache_rnasep692 1 $R1/node_r1/poses $R1/node_r1/disc/site_01.mrc $MEANR1 $R1/fits/r1_s1@r1/hi 10 $R1/fits/r1_s1@r1/hi/em -
CASES
run(){ nm=$1; fr=$2; read _ H CP CA PO MK CR SD IT FULL TC <<< "$(grep "^$nm " $O/cases.txt)"; D=$O/${nm}_f$fr; mkdir -p $D
  [ -s $D/em/post_soft_B.txt ] && return 0
  $PY -c "
import numpy as np
rng = np.random.default_rng(1)
for h in 'AB': s = np.loadtxt('$SD/start_%s.txt' % h, dtype=int); np.savetxt('$D/fit_%s.txt' % h, ((s > 0) & (rng.random(len(s)) < $fr)).astype(int), fmt='%d')"
  $PY -u $T/minor_em3s.py --sim ${H}_A,${H}_B --poses ${PO}_A.txt,${PO}_B.txt --labels $SD/start_A.txt,$SD/start_B.txt --cache $CP --ctfabs $CA --cons_ref $CR --mask $MK --soft 1 --iters $IT \
    --em_lim 15 --post 1 --nomaps 1 --fit $D/fit_A.txt,$D/fit_B.txt --out $D/em > $D/em.log 2>&1 || echo "$nm $fr FAILED: $(tail -2 $D/em.log)"; }
export -f run; export O PY T
JOBS="toyC_refined:1.0"; for nm in toyA toyC_refined toyC_true r1_1f1_lo r1_1f1_hi; do for fr in ${FRACS:-0.25 0.5}; do JOBS="$JOBS $nm:$fr"; done; done
echo $JOBS | tr ' ' '\n' | xargs -P ${PAR:-4} -I{} bash -c 'j={}; run ${j%%:*} ${j#*:}'
$PY - <<PYEOF
import numpy as np, os
O = '$O'; C = {l.split()[0]: l.split() for l in open(O + '/cases.txt')}; TJ = '$K/factorial/junk_test/simsites_truth_jt.txt'
hd = [l.split()[1:] for l in open(TJ) if l.startswith('#')][0]; TT = np.loadtxt(TJ)
rk = lambda x: np.argsort(np.argsort(x)).astype(float)
def sp(x, y): return np.corrcoef(rk(x), rk(y))[0, 1]
def auc(s, y): r_ = rk(s) + 1.; n1 = y.sum(); n0 = len(y) - n1; return (r_[y].sum() - n1 * (n1 + 1) / 2) / max(n1 * n0, 1)
def master(d, H, what):
    NM = 199996 if 'simsites' in H else 932496; v = np.full(NM, np.nan)
    for h in 'AB':
        r = np.loadtxt('%s_%s/meta_rows.txt' % (H, h)).astype(int) - 1; x = np.loadtxt('%s/%s_%s.txt' % (d, what, h)); v[r] = x[:, 1] if x.ndim == 2 else x
    return v
for d in sorted(x for x in os.listdir(O) if os.path.isdir(O + '/' + x)):
    nm, fr = d.rsplit('_f', 1); c = C[nm]; H, SD, FULL, TC = c[1], c[7], c[9], c[10]
    if not os.path.exists('%s/%s/em/post_soft_B.txt' % (O, d)): print(d, 'missing'); continue
    full, sub, st, fit = master(FULL, H, 'post_soft'), master(O + '/' + d + '/em', H, 'post_soft'), master(SD, H, 'start'), master(O + '/' + d, H, 'fit')
    m = st > 0; ho = m & (fit == 0)
    s = '%-14s fit %4.0f %% (%6d of %6d): Spearman all %.3f, held-out %.3f' % (nm, 100 * float(fr), (fit[m] > 0).sum(), m.sum(), sp(sub[m], full[m]), sp(sub[ho], full[ho]) if ho.sum() > 10 else float('nan'))
    for N in (2500, 10000):
        a_ = set(np.flatnonzero(m)[np.argsort(-sub[m])[:N]]); b_ = set(np.flatnonzero(m)[np.argsort(-full[m])[:N]]); s += ', top-%d overlap %.3f' % (N, len(a_ & b_) / N)
    s += ', class-2 share %.3f vs full %.3f' % ((sub[m] >= 0.5).mean(), (full[m] >= 0.5).mean())
    if TC != '-':
        y = np.zeros(len(TT), bool); y[TT[:, 0].astype(int) - 1] = (TT[:, hd.index(TC)] > 0) & (TT[:, hd.index('junk')] == 0); g = m.copy(); g[TT[:, 0].astype(int) - 1] &= TT[:, hd.index('junk')] == 0
        s += ' | %s AUC sub %.3f full %.3f' % (TC, auc(sub[g], y[g]), auc(full[g], y[g]))
    print(s)
PYEOF
echo "[$(date '+%F %T') subsample test] DONE"
