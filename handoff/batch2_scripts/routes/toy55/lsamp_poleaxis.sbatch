#!/bin/bash
#SBATCH -J klem_lsamp_poleaxis
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# (1) pole-axis coordinates (lsamp_poleaxis.py: closed-pole minus open-pole score, from the poles_nstar score spaces) -> equal-N bins across
# runs (lsamp_repro_bins2.py, held-out maps, arm readouts, comparison); (2) view / defocus dependence of every landscape coordinate tested
# tonight (lsamp_viewcheck.py). env W (repro_r1), K_ (rnasep_klem)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; W=${W:-$K_/lsamp/repro_r1}; O=$W/poleaxis; mkdir -p $O/bins
H=$T/simrnasep692; MR=${H}_A/meta_rows.txt,${H}_B/meta_rows.txt; CP=$K_/lsamp/caches/rnasep692; CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
PS_script=$K_/lsamp/f692st_ps/poses_combined.txt; PS_port1=$W/port1_ps/poses_combined.txt; PS_port2=$W/port2_ps/poses_combined.txt; PS_port3=$W/port3_ps/poses_combined.txt
for r in script port1 port2; do pre=r_; [ $r = script ] && pre=r1_; ps=PS_$r; ps=${!ps}; SC=$W/poles_nstar/${r}_ps
  $PY -u $T/lsamp_poleaxis.py --scores $SC --closed ${pre}pA --open ${pre}pB --out $O/$r 2>&1 | grep -v "Warning\|warn("
  $PY -u $T/lsamp_repro_bins2.py --name $r --rundir $O/$r --scores $SC --poses $ps --ref_rundir $O/script --sim ${H}_A,${H}_B --mr $MR --cache $CP --ctfabs 1 --mskdiam 232 --out $O/bins 2>&1 | grep -v "Warning\|warn("
  for f in A B; do $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/bins/$r/arm_bins_$f.txt $O/bins/$r/bins_$f 2>&1 | grep -v Warning > /dev/null; done; done
$PY -u $T/lsamp_repro_bins2.py --cmp 1 --names script,port1,port2 --out $O/bins 2>&1 | grep -v "Warning\|warn(\|def lp" > $O/bins/bins_summary.txt; cat $O/bins/bins_summary.txt
V=""; for r in script port1 port2 port3; do ps=PS_$r; ps=${!ps}; V="$V axis1_$r:$W/$r:$ps axis2_$r:$W/${r}_ax2:$ps"; done
for r in script port1 port2; do ps=PS_$r; ps=${!ps}; V="$V polesNstar_$r:$W/poles_nstar/$r:$ps poles2500_$r:$W/poles_top2500/$r:$ps poleaxis_$r:$O/$r:$ps"; done
$PY -u $T/lsamp_viewcheck.py $V 2>&1 | grep -v "Warning\|warn(" | tee $O/viewcheck.txt
echo "[$(date '+%F %T')] DONE"
