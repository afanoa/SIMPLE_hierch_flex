#!/bin/bash
#SBATCH -J klem_lsamp_pqf
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# three-pole coordinates of the factorial-pole runs (script / port1 / port2 poles_nstar score spaces: c = pA closed, oA = pB z1-open,
# oB = pD z3-open; lsamp_pq_coords.py) for one coordinate TYPE, then their equal-N held-out bins (lsamp_repro_bins2.py, oriented by script).
# env TYPE (c_oA, c_oB, oA_oB or c_rest)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; W=$K_/lsamp/repro_r1; O=$W/pq; TY=${TYPE:?}
H=$T/simrnasep692; MR=${H}_A/meta_rows.txt,${H}_B/meta_rows.txt; CP=$K_/lsamp/caches/rnasep692; CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
for r in script port1 port2; do pre=r_; ps=$W/${r}_ps/poses_combined.txt; [ $r = script ] && { pre=r1_; ps=$K_/lsamp/f692st_ps/poses_combined.txt; }; SC=$W/poles_nstar/${r}_ps
  $PY -u $T/lsamp_pq_coords.py --scores $SC --name $r --out $O --roles ${pre}pA,${pre}pB,${pre}pD --types $TY 2>&1 | grep -v "Warning\|warn("
  $PY -u $T/lsamp_repro_bins2.py --name $r --rundir $O/$TY/$r --scores $SC --poses $ps --ref_rundir $O/$TY/script --sim ${H}_A,${H}_B --mr $MR --cache $CP --ctfabs 1 --mskdiam 232 --out $O/$TY/bins 2>&1 | grep -v "Warning\|warn("
  for f in A B; do $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/$TY/bins/$r/arm_bins_$f.txt $O/$TY/bins/$r/bins_$f 2>&1 | grep -v Warning > /dev/null; done; done
echo "[$(date '+%F %T')] DONE $TY"
