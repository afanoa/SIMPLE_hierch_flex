#!/bin/bash
#SBATCH -J klem_port_sites
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 10:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# factorial site post-pass (tree_sites.py) on a PORT hetero_tree, read-only: a view in the rflex layout (node.done, n.simple, mean map,
# rows_<u>.txt 0/1, final_labels = hetero_tree_labels) with the port's own per-node site candidates as the discovery -> per-site flex,
# two-tailed poles, pooling -> scoring vs truth (SCORE: sites = score_port_sites.py, flex = score_site_flex.py). env PT OUT H MRA MRB CP CTFABS MSKDIAM NSN GATE TRUTH TRUTHJT
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=${QPART:-csbdevel} OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=4
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; O=${OUT:?}; V=$O/view; MD=${MSKDIAM:-232}; mkdir -p $V; cd $O
$PY - <<PYEOF || exit 1
import numpy as np, os, glob
P, V, O = '${PT:?}', '$V', '$O'; NM = sum(1 for _ in open(P + '/hetero_tree_labels.txt'))
for D in sorted(glob.glob(P + '/node_*')):
    u = D.split('/node_')[-1]; need = [D + f for f in ('/n.simple', '/rows.txt', '/sites/sites.txt', '/mean/recvol_state01.mrc')]
    if not all(os.path.exists(f) for f in need): print('node %s skipped (missing %s)' % (u, ' '.join(os.path.basename(f) for f in need if not os.path.exists(f)))); continue
    os.makedirs(V + '/node_%s/mean/1_reconstruct3D' % u, exist_ok=True); open(V + '/node_%s/node.done' % u, 'w').close()
    for src, dst in ((D + '/n.simple', V + '/node_%s/n.simple' % u), (D + '/mean/recvol_state01.mrc', V + '/node_%s/mean/1_reconstruct3D/recvol_state01.mrc' % u)):
        if not os.path.lexists(dst): os.symlink(src, dst)
    r = np.loadtxt(D + '/rows.txt', dtype=int); L = np.zeros(NM, int); L[r - 1] = 1; np.savetxt(V + '/rows_%s.txt' % u, L, fmt='%d')
    dd = O + '/node_%s/disc' % u; os.makedirs(dd, exist_ok=True)
    for f in glob.glob(D + '/sites/*'):
        t = dd + '/' + os.path.basename(f)
        if not os.path.lexists(t): os.symlink(f, t)
    print('node %s: %d rows, %d port site candidates' % (u, L.sum(), len(glob.glob(D + '/sites/site_*.mrc'))))
if not os.path.lexists(V + '/final_labels.txt'): os.symlink(P + '/hetero_tree_labels.txt', V + '/final_labels.txt')
PYEOF
$PY -u $T/tree_sites.py --tree $V --out $O --h ${H:?} --cache ${CP:?} ${MRA:+--mra $MRA} ${MRB:+--mrb $MRB} --ctfabs ${CTFABS:-0} --mskdiam $MD --nsn ${NSN:-8} --gate ${GATE:-0.5} \
  --pars ${PARS:-4} --parf ${PARF:-2} --np ${NP:-4} || { echo "tree_sites FAILED"; exit 1; }
[ -s $O/sites_final.txt ] || exit 0
[ -n "${TRUTHJT:-}" ] && $PY $T/score_port_sites.py --out $O --truth $TRUTHJT --tree_labels $V/final_labels.txt 2>&1 | grep -v Warning | tee $O/score_sites.txt
if [ -n "${TRUTH:-}" ]; then FL=$(awk '!/^#/ {printf "%s=%s/fits/%s@%s/flex ", $1, o, $1, $2}' o=$O $O/sites_final.txt); EM=$(awk '!/^#/ {printf "%s_lo=%s/fits/%s@%s/lo %s_hi=%s/fits/%s@%s/hi ", $1, o, $1, $3, $1, o, $1, $3}' o=$O $O/sites_final.txt)
  LB="tree=$V/final_labels.txt"; [ -s $O/hybrid_labels.txt ] && LB="$LB treesites=$O/hybrid_labels.txt"
  $PY $T/score_site_flex.py --truth $TRUTH --rows $V/rows_r.txt --flex $FL --em $EM --labels $LB --json $O/score_flex.json 2>&1 | grep -v Warning | tee $O/score_flex.txt; fi
echo "[$(date '+%F %T') port tree sites] DONE"
