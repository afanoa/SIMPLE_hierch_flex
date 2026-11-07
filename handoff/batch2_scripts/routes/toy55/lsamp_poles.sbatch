#!/bin/bash
#SBATCH -J klem_lsamp_poles
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 02:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# landscape stage with pole maps in each run's score space, then the cross-run comparison. Per run: the node's leaves (leaf-name prefix) +
# the run's pole maps (lsamp_polemaps.py: top-N of each pole's posterior, at the run's poses, amplitude-matched) -> scores (lsamp_score2.sbatch)
# -> lsamp2h.py with the defaults -> arm readouts -> equal-N bins (lsamp_repro_bins2.py, oriented by the first run); at the end the bins
# comparison (bins/bins_summary.txt) and the per-sample comparison (repro_summary.txt).
# env RUNS (space-separated NAME:TREE:POSES:PREFIX:POLE=ENDPOINT+POLE=ENDPOINT..), PS (pole-set dir: ENDPOINT/post_master.txt + nstar.txt),
# NTOP (star = each endpoint's N*, or a number), OUT
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; O=${OUT:?}; mkdir -p $O/bins
H=$T/simrnasep692; MR=${H}_A/meta_rows.txt,${H}_B/meta_rows.txt; CP=$K_/lsamp/caches/rnasep692; CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
log(){ echo "[$(date '+%F %T') lsamp_poles] $*"; }; NAMES=""; REF=""; CMP=""
for spec in ${RUNS:?}; do IFS=: read -r nm tree ps pre poles <<< "$spec"; D=$O/${nm}_tree; SC=$O/${nm}_ps; OD=$O/$nm; mkdir -p $D $OD
  [ -n "$REF" ] || REF=$OD; NAMES="$NAMES$nm,"; CMP="$CMP $nm:$OD"
  awk -v p=$pre -v d=$D 'FNR == NR { if (index($1, p) == 1) { m[FNR] = ++k; print > (d "/base_leaves.txt") }; next } { print (($1 in m) ? m[$1] : 0) > (d "/base_labels.txt") }' $tree/leaves_list.txt $tree/final_labels.txt
  PL=""; for e in $(echo $poles | tr '+' ' '); do pn=${e%%=*}; ep=${e#*=}; nt=${NTOP:-star}; [ $nt = star ] && nt=$(cat ${PS:?}/$ep/nstar.txt | tr -d ' '); PL="$PL,$pn:$PS/$ep/post_master.txt:$nt"; done
  log "$nm: $(grep -c . $D/base_leaves.txt) leaves ($pre), poles ${PL#,}"
  [ -s $D/leaves_list.txt ] || $PY -u $T/lsamp_polemaps.py --leaves $D/base_leaves.txt --labels $D/base_labels.txt --poses $ps --poles ${PL#,} --prefix $pre --out $D \
      --sim ${H}_A,${H}_B --mr $MR --cache $CP --ctfabs 1 --mskdiam 232 2>&1 | grep -v "Warning\|warn("
  [ -s $SC/meta.npz ] || { TREE=$D OUT=$SC SIM=${H}_A,${H}_B MR=$MR CP=$CP CTFABS=1 EMLIM=12 MSKDIAM=232 POSES=$ps bash $T/lsamp_score2.sbatch 2>&1 | grep "own leaf ranks\|score band\|DONE\|Error\|Traceback"; }
  $PY -u $T/lsamp2h.py --scores $SC --sim ${H}_A,${H}_B --mr $MR --cache $CP --poses $ps --ctfabs 1 --mskdiam 232 --out $OD --use $pre --naxes 0 --ends centre --tag $nm 2>&1 \
    | grep --line-buffered -v "Warning\|warn(" > $OD/lsamp.log
  grep "axes (lam\|^  *[1-4]  \|sampled axes\|dropping\|resolvable samples\|CHOSEN\|steps CV\|maps (cells\|Error\|Traceback" $OD/lsamp.log | cut -c1-200
  for d in heldout_A heldout_B final_all; do $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $OD/arm_$d.txt $OD/$d 2>&1 | grep -v Warning > /dev/null
    echo "  $d: $(awk '!/^#/ {k = $1; sub(/.*state_/, "", k); sub(/\.mrc/, "", k); printf "%s %.2f/%.2f  ", k + 0, $2, $3}' $OD/arm_$d.txt)"; done
  $PY -u $T/lsamp_repro_bins2.py --name $nm --rundir $OD --scores $SC --poses $ps --ref_rundir $REF --sim ${H}_A,${H}_B --mr $MR --cache $CP --ctfabs 1 --mskdiam 232 --out $O/bins 2>&1 | grep -v "Warning\|warn("
  for f in A B; do $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/bins/$nm/arm_bins_$f.txt $O/bins/$nm/bins_$f 2>&1 | grep -v Warning > /dev/null; done
done
$PY -u $T/lsamp_repro_bins2.py --cmp 1 --names $NAMES --out $O/bins 2>&1 | grep -v "Warning\|warn(\|def lp" > $O/bins/bins_summary.txt; cat $O/bins/bins_summary.txt
$PY -u $T/lsamp_repro_cmp2.py $O/repro_summary.txt $CMP 2>&1 | grep -v "Warning\|warn(\|def prep" > /dev/null; cat $O/repro_summary.txt
log DONE
