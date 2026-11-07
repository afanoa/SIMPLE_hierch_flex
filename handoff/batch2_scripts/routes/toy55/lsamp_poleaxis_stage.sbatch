#!/bin/bash
#SBATCH -J klem_lsamp_pastage
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# the landscape stage on the pole axis (lsamp2p.py --coords repro_r1/poleaxis/<run>/coords.npz --axes 1, the poles_nstar score spaces):
# sample count, spacing, checks and held-out maps per run, arm readouts, then the per-sample comparison across runs (lsamp_repro_cmp2.py).
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; W=$K_/lsamp/repro_r1; O=$W/poleaxis_stage; mkdir -p $O
H=$T/simrnasep692; MR=${H}_A/meta_rows.txt,${H}_B/meta_rows.txt; CP=$K_/lsamp/caches/rnasep692; CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1); CMP=""
for r in script port1 port2; do pre=r_; ps=$W/${r}_ps/poses_combined.txt; [ $r = script ] && { pre=r1_; ps=$K_/lsamp/f692st_ps/poses_combined.txt; }; mkdir -p $O/$r; CMP="$CMP $r:$O/$r"
  $PY -u $T/lsamp2p.py --scores $W/poles_nstar/${r}_ps --sim ${H}_A,${H}_B --mr $MR --cache $CP --poses $ps --ctfabs 1 --mskdiam 232 --out $O/$r --use $pre --axes 1 \
      --coords $W/poleaxis/$r/coords.npz --ends centre --tag ${r}_poleaxis 2>&1 | grep --line-buffered -v "Warning\|warn(" > $O/$r/lsamp.log
  echo "== $r"; grep "coordinate from\|resolvable samples\|CHOSEN\|steps CV\|check 3\|maps (cells\|Error\|Traceback" $O/$r/lsamp.log | cut -c1-220
  for d in heldout_A heldout_B final_all; do $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/$r/arm_$d.txt $O/$r/$d 2>&1 | grep -v Warning > /dev/null
    echo "  $d: $(awk '!/^#/ {k = $1; sub(/.*state_/, "", k); sub(/\.mrc/, "", k); printf "%s %.2f/%.2f  ", k + 0, $2, $3}' $O/$r/arm_$d.txt)"; done; done
$PY -u $T/lsamp_repro_cmp2.py $O/repro_summary.txt $CMP 2>&1 | grep -v "Warning\|warn(\|def prep" > /dev/null; cat $O/repro_summary.txt
echo "[$(date '+%F %T')] DONE"
