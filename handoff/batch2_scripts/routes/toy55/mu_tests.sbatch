#!/bin/bash
#SBATCH -J klem_mutest
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=110G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# merge up the tree: the gold-standard split tests of one level's nodes (LEVEL dir of mu_prep.py) with index % NCHUNK == CHUNK,
# PAR at a time: splitbench8 on the half sets SIMP (real RNase P halves; CACHE, CTFABS, GTX for a synthetic set), family-refined poses (POSES prefix), the node's focus mask, cross-half EM (XEM)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
one(){ d=$1; m=$2; [ -s $d/children_B.txt ] && return 0
  $PY -u $W/splitbench8.py --sim ${SIMP:-$W/simrnasep692}_A,${SIMP:-$W/simrnasep692}_B --poses ${POSES:?}_A.txt,${POSES}_B.txt --cand $d/cand_A.txt,$d/cand_B.txt --out $d \
    --cache ${CACHE:-$W/benchcache_rnasep692} --ctfabs ${CTFABS:-1} ${GTX:-} --pair identity --xem ${XEM:-3} --softem 0 $( [ "$m" != none ] && [ -s "$m" ] && echo "--focus $m" ) > $d/bench.log 2>&1 \
    && echo "$(basename $d): $(grep -c . $d/labels.txt | awk '{print $1 - 1}') units -> K* $(grep -a 'K\* =' $d/bench.log | tail -1 | grep -o '[0-9]*' | head -1)" || echo "$(basename $d) FAILED: $(grep -a -m1 'Error' $d/bench.log)"; }
i=0; while read d m us; do case "$d" in \#*) continue;; esac; [ -s $d/cand_A.txt ] || continue
  if [ $((i % ${NCHUNK:-1})) = ${CHUNK:-0} ]; then while [ $(jobs -rp | wc -l) -ge ${PAR:-2} ]; do sleep 5; done; one $d $m & fi; i=$((i + 1)); done < ${LEVEL:?}/nodes.txt
wait; echo "MUTEST CHUNK ${CHUNK:-0} DONE"
