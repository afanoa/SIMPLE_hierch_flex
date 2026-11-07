#!/bin/bash
#SBATCH -J klem_arm_childpoles
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:15:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# closed-arm / open-arm readouts (10 A) of the closed child's top-N* pole maps from the per-child discovery check (stree_port/childdisc)
export OMP_NUM_THREADS=2; PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; C=/mnt/beegfs/elmlund/afan/rnasep_klem/stree_port/childdisc
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
for r in p1 p2; do for e in z1lo z1hi z3lo z3hi; do d=$C/${r}_r1/top_$e; [ -d $d ] || { echo "$r $e: no directory"; continue; }
  m=$(ls $d/recvol_state01_even_unfil.mrc 2>/dev/null || ls $d/recvol_state01.mrc 2>/dev/null); [ -n "$m" ] || m=$(ls -t $d/*/recvol_state01_even_unfil.mrc 2>/dev/null | head -1)
  [ -n "$m" ] || { echo "$r $e: no map"; continue; }
  echo "$r $e $(cat $d/nstar.txt 2>/dev/null) $($PY $T/arm_readout_any.py $T/maps160/toy_consensus.mrc $CL $m 2>&1 | grep -v Warning | tail -1)"; done; done
