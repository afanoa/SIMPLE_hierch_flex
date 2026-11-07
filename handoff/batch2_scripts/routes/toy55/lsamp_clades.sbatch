#!/bin/bash
#SBATCH -J klem_lsamp_clades
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# landscape sampling of a FLAT leaf set (one K-state refinement, no tree): complete-linkage tree over the states' half-maps (state_merge.py),
# clades of the cut CUT (4) with >= 2 states play the terminal nodes; lsamp2.py on each clade (lsamp_score2 scores of all states).
# env RD (refine3D_states dir), LAB (state labels, master rows), SC, OUT, SIM, MR, CP, POSES, CTFABS, MSKDIAM, GT, ARGS (extra lsamp2.py flags), CUT
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; O=${OUT:?}; CUT=${CUT:-4}; mkdir -p $O
[ -s $O/smerge/merge_tree.txt ] || $PY $T/state_merge.py $O/smerge --rdir ${RD:?} --labels ${LAB:?} ${GT:+--gt $GT} --unfil --link complete --ncut 8,$CUT,2 2>&1 | grep -v Warning | grep "^== " | cut -c1-200
awk '!/^#/ && NF >= 3' $O/smerge/cut$CUT/groups.txt | while read g n mem; do
  k=$(echo $mem | tr '+' '\n' | wc -l); [ $k -ge 2 ] || continue
  use=$(echo $mem | tr '+' '\n' | awk '{printf "%ss%02d", (NR > 1 ? "," : ""), $1}'); D=$O/clade$g; mkdir -p $D
  echo "== clade $g: $n particles, states $mem"
  $PY -u $T/lsamp2.py --scores ${SC:?} --sim ${SIM:?} --mr ${MR:?} --cache ${CP:?} --poses ${POSES:?} --ctfabs ${CTFABS:-0} --mskdiam ${MSKDIAM:-220} --out $D --use $use \
      --naxes 0 --tag clade$g ${GT:+--gt $GT} ${ARGS:-} 2>&1 | grep --line-buffered -v "Warning\|warn(" > $D/lsamp.log
  grep "sampled axes\|junk screen\|resolvable\|CHOSEN\|steps CV\|TRUTH .*span\|TRUTH |corr\|TRUTH best\|TRUTH closed\|no landscape\|  sample .*coded" $D/lsamp.log | cut -c1-300
done
echo "DONE"
