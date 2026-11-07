#!/bin/bash
#SBATCH -J klem_lsamp_bins
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 01:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# equal-N comparison of independent tree runs' landscapes of one node (lsamp_repro.sbatch outputs): per run the same quantile bins of its
# fold coordinates -> maps from the held-out pixels (lsamp_repro_bins.py), arm readouts, then the cross-run comparison (--cmp 1).
# env RUNS (space-separated NAME:RUNDIR:SCOREDIR:POSES), OUT
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; O=${OUT:?}; mkdir -p $O
H=$T/simrnasep692; MR=${H}_A/meta_rows.txt,${H}_B/meta_rows.txt; CP=$K_/lsamp/caches/rnasep692; CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1); NAMES=""
for spec in ${RUNS:?}; do IFS=: read -r nm rd sc ps <<< "$spec"; NAMES="$NAMES$nm,"
  echo "[$(date '+%F %T')] $nm"; $PY -u $T/lsamp_repro_bins.py --name $nm --rundir $rd --scores $sc --poses $ps --sim ${H}_A,${H}_B --mr $MR --cache $CP --ctfabs 1 --mskdiam 232 --out $O 2>&1 | grep -v "Warning\|warn("
  for f in A B; do $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/$nm/arm_bins_$f.txt $O/$nm/bins_$f 2>&1 | grep -v Warning > /dev/null; done
done
$PY -u $T/lsamp_repro_bins.py --cmp 1 --names $NAMES --out $O 2>&1 | grep -v "Warning\|warn(" > $O/bins_summary.txt; cat $O/bins_summary.txt
echo "[$(date '+%F %T')] DONE"
