#!/bin/bash
#SBATCH -J lsamp_psub
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# the Python stage (lsamp2p.py --coords, closed - open A pole axis) on the particle subset the Fortran landscape_sample pair uses (state 3 =
# open B left out), for a like-for-like check of the stage logic against repro_r1/fortran/<run>/landscape_1_2.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; W=$K_/lsamp/repro_r1
H=$T/simrnasep692; MR=${H}_A/meta_rows.txt,${H}_B/meta_rows.txt; CP=$K_/lsamp/caches/rnasep692; CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
for p in p1 p2; do O=$W/pq/stage_c_oA_sub/$p; mkdir -p $O
  awk '{print ($1 == 3) ? 1 : 0}' $W/pq/${p}_tree/final_labels.txt > $O/exclude.txt
  $PY -u $T/lsamp2p.py --scores $W/pq/${p}_ps --sim ${H}_A,${H}_B --mr $MR --cache $CP --poses $W/pq/${p}_tree/dump.txt --ctfabs 1 --mskdiam 232 \
      --out $O --use r_ --axes 1 --coords $W/pq/c_oA/$p/coords.npz --exclude_rows $O/exclude.txt --ends centre --tag ${p}_sub 2>&1 | grep --line-buffered -v "Warning\|warn(" > $O/lsamp.log
  echo "== $p"; grep "excluded rows\|resolvable samples\|CHOSEN\|steps CV\|check 3\|fold . at n\|neighbour steps" $O/lsamp.log | cut -c1-230
  for d in heldout_A heldout_B; do $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/arm_$d.txt $O/$d 2>&1 | grep -v Warning > /dev/null
    echo "  $d: $(awk '!/^#/ {k = $1; sub(/.*state_/, "", k); sub(/\.mrc/, "", k); printf "%s %.2f/%.2f  ", k + 0, $2, $3}' $O/arm_$d.txt)"; done
  grep "final_all maps" $O/lsamp.log | grep -o "sizes.*"
done
echo DONE
