#!/bin/bash
#SBATCH -J klem_lsamp_bins2
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 01:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# lsamp_repro_bins.sbatch with every run oriented by the first one (lsamp_repro_bins2.py --ref_rundir; lsamp_repro_cmp2.py): equal-N bins of
# each run's fold coordinates -> maps from the held-out pixels, arm readouts, the bins comparison and the per-sample comparison.
# env RUNS (space-separated NAME:RUNDIR:SCOREDIR:POSES; the first run sets the orientation), OUT
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; O=${OUT:?}; mkdir -p $O
H=$T/simrnasep692; MR=${H}_A/meta_rows.txt,${H}_B/meta_rows.txt; CP=$K_/lsamp/caches/rnasep692; CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1); NAMES=""; REF=""; CMP=""
for spec in ${RUNS:?}; do IFS=: read -r nm rd sc ps <<< "$spec"; NAMES="$NAMES$nm,"; CMP="$CMP $nm:$rd"; [ -n "$REF" ] || REF=$rd
  echo "[$(date '+%F %T')] $nm"; $PY -u $T/lsamp_repro_bins2.py --name $nm --rundir $rd --scores $sc --poses $ps --ref_rundir $REF --sim ${H}_A,${H}_B --mr $MR --cache $CP --ctfabs 1 --mskdiam 232 --out $O 2>&1 | grep -v "Warning\|warn("
  for f in A B; do $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/$nm/arm_bins_$f.txt $O/$nm/bins_$f 2>&1 | grep -v Warning > /dev/null; done
done
$PY -u $T/lsamp_repro_bins2.py --cmp 1 --names $NAMES --out $O 2>&1 | grep -v "Warning\|warn(\|def lp" > $O/bins_summary.txt; cat $O/bins_summary.txt
$PY -u $T/lsamp_repro_cmp2.py $O/repro_summary.txt $CMP 2>&1 | grep -v "Warning\|warn(\|def prep" > /dev/null; cat $O/repro_summary.txt
echo "[$(date '+%F %T')] DONE"
