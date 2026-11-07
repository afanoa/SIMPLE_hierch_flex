#!/bin/bash
#SBATCH -J land_leafpath
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 03:00:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_leafpath.%j.out
# side readout for the extremes test: a path between one run's two extreme leaves on one pooled stable axis (land_leafext.sbatch views,
# lo -> hi), OT morph inside the node variance region, r1 scored at the shared input poses of node_r1/n.simple, reweighting + basins
# (lam 0.001 / 0.01, 0.6 kT, half check). args: RUN AXIS
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; L=$K_/landscape; H=$T/simrnasep692
R=${1:?}; X=${2:?}; O=$L/leafpath_${R}_ax$X; V=$K_/view_rflex/leafext_r1; mkdir -p $O && cd $O; log(){ echo "[$(date '+%F %T') leafpath $R ax$X] $*"; }
grep "state_${R}_ax${X}_lo.mrc" $L/leafext/arm.txt > $O/arm.txt; grep "state_${R}_ax${X}_hi.mrc" $L/leafext/arm.txt >> $O/arm.txt; log "ends: $(cut -c1-200 $O/arm.txt | tr '\n' '|')"
[ -s $O/morph/path_30.mrc ] || $PY -u $T/land_morph2.py $V/state_${R}_ax${X}_lo.mrc $V/state_${R}_ax${X}_hi.mrc $O/morph --mask $K_/rflex_f692st/node_r1/rsmask_variance.mrc --mask_thr 0.2 2>&1 | grep --line-buffered -v Warning
ls $O/morph/path_*.mrc > $O/maps.txt
[ -s $O/ps.npz ] || $PY -u $T/land_pscore.py --sim ${H}_A,${H}_B --mr ${H}_A/meta_rows.txt,${H}_B/meta_rows.txt --cache $K_/lsamp/caches/rnasep692 \
   --poses $L/leafext/poses_n.txt --rowflags $K_/rflex_f692st/rows_r1.txt --ctfabs 1 --mskdiam 232 --maps $O/maps.txt --out $O/ps.npz 2>&1 | grep --line-buffered -v "Warning\|warn(\|block"
for lam in 0.001 0.01; do log "basins lam $lam"; $PY -u $T/land_basins.py --scores $O/ps.npz --lam $lam --thr 0.5,0.6,0.75 2>&1 | grep --line-buffered -v Warning; done
$PY -u $T/land_reweight.py --scores $O/ps.npz --out $O/rw_lmarg_split_0.001 --real --variant lmarg --tau split --lam 0.001 --nm 932496 --nboot 20 2>&1 | grep --line-buffered -v "Warning\|warn("
log "DONE: $O"
