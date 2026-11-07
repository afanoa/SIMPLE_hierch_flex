#!/bin/bash
#SBATCH -J klem_dscore
#SBATCH -p csbdevel
#SBATCH -c 24
#SBATCH --mem=120G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# particle-level width rule test for node f692st r1: per run (RUNS: script port1 port2 [port3]) every r1 particle scored cross-half against
# that run's own 8 stage-B leaf maps at that run's stage-B poses (score_embed.py, band to 15 A), then decisive_merge.py merges leaves until
# every remaining pair is decided per particle and compares the merged partitions across runs. Outputs rnasep_klem/repro_r1/dscore_<run>/,
# report rnasep_klem/repro_r1/decisive_merge.txt. env RUNS
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$PATH OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 NUMBA_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; H=$T/simrnasep692
log(){ echo "[$(date '+%F %T') dscore] $*"; }
score_run(){ local run=$1 RD O=$K_/repro_r1/dscore_$1
  case $run in script) RD=$K_/rflex_f692st/node_r1/fstageB/2_refine3D;; port1) RD=$K_/stree_port/tree_f692r1/1_hetero_tree/node_r/stageB;;
    port2) RD=$K_/stree_port/tree_f692r1_seed2/1_hetero_tree/node_r/stageB;; port3) RD=$K_/stree_port/tree_f692r1_rep3/1_hetero_tree/node_r/stageB;; esac
  mkdir -p $O; [ -s $O/scores_B.npy ] && { log "$run: scores exist"; return 0; }
  simple_exec prg=print_project_field projfile=$RD/b.simple oritype=ptcl3D mkdir=no > $O/dump.txt 2> /dev/null || { log "$run: dump FAILED"; return 1; }
  awk '{s = 0; for (i = 1; i <= NF; i++) if ($i ~ /^state=/) {split($i, a, "="); s = a[2] + 0}; print s}' $O/dump.txt > $O/labels.txt
  for h in A B; do awk 'NR == FNR {l[FNR] = $0; next} {print l[$1 + 0]}' $O/dump.txt ${H}_$h/meta_rows.txt > $O/poses_$h.txt; done
  : > $O/leaves.txt; for j in $(seq 1 8); do jj=$(printf %02d $j); echo "s$jj $RD/recvol_state${jj}_even_unfil.mrc $RD/recvol_state${jj}_odd_unfil.mrc $(awk -v j=$j '$1 == j' $O/labels.txt | wc -l)" >> $O/leaves.txt; done
  $PY -u $T/score_embed.py --sim ${H}_A,${H}_B --cache $T/benchcache_rnasep692 --poses $O/poses_A.txt,$O/poses_B.txt --leaves $O/leaves.txt \
    --labels $O/labels.txt --out $O --em_lim 15 --nsub 5000 > $O/score.log 2>&1 || { log "$run: scoring FAILED: $(grep -a -m2 'Error\|Traceback' $O/score.log)"; return 1; }
  log "$run: $(grep -a 'own leaf ranks first' $O/score.log | tr '\n' ' ')"; }
for r in ${RUNS:?}; do score_run $r & done; wait
$PY $T/decisive_merge.py $RUNS --thr 2,3,4 > $K_/repro_r1/decisive_merge.txt 2>&1; cat $K_/repro_r1/decisive_merge.txt
log DONE
