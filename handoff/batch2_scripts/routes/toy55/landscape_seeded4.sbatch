#!/bin/bash
#SBATCH -J klem_landscape_seeded4
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# map-level reproducibility of the pole-seeded K=4 stage-A states (z1-closed, open A, open B, z3-closed) across script / port1 / port2
# poles, and per run K=4 against that run's K=3 (is the z3-closed state a split of the K=3 closed state?) (landscape_cmp.py)
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8; K=/mnt/beegfs/elmlund/afan/rnasep_klem; PY=/mnt/nasapps/production/python/3.12/bin/python3
L=/mnt/beegfs/elmlund/afan/routes/toy55/landscape_cmp.py; a(){ d=$K/repro_r1/seeded_$1_$2; echo "$3=$(ls -td $d/stageA/*_refine3D | head -1):$d/stageA_states.txt"; }
{ echo "######## pole-seeded p4 across runs"; $PY $L $(a p4 script script) $(a p4 port1 port1) $(a p4 port2 port2) --lp 15 2>&1
  for r in script port1 port2; do echo "######## $r: p4 vs p3"; $PY $L $(a p4 $r p4) $(a p3 $r p3) --lp 15 2>&1; done
  T=/mnt/beegfs/elmlund/afan/routes/toy55; CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
  for p in p3 p4; do for r in script port1 port2; do d=$K/repro_r1/seeded_${p}_$r; R=$(ls -td $d/stageA/*_refine3D | head -1)
    echo "######## arm readouts (10 A, stage-A box 64 padded to 160) $p $r, state sizes $(awk '$1 > 0' $d/stageA_states.txt | sort -n | uniq -c | awk '{printf "%d ", $1}')"
    $PY $T/arm_readout_any.py $T/maps160/toy_consensus.mrc $CL $R/recvol_state0?_even_unfil.mrc 2>&1 | grep -v Warning; done; done; } > $K/repro_r1/landscape_seeded4.txt; cat $K/repro_r1/landscape_seeded4.txt
