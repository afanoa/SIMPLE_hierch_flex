#!/bin/bash
#SBATCH -J klem_arm_d1
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# closed-arm / open-arm readouts (10 A) of the depth-1 recursion test (poles_r1d1_p1/p2): each r1 child's stage-B states and the delivered maps
export OMP_NUM_THREADS=2; PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; P=/mnt/beegfs/elmlund/afan/rnasep_klem/stree_port
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
for r in p1 p2; do H=$P/poles_r1d1_$r/1_hetero_tree
  for n in $(ls -d $H/node_r? 2>/dev/null); do B=$n/stageB; [ -n "$(ls $B/recvol_state0?_even_unfil.mrc 2>/dev/null)" ] || continue
    echo "######## $r $(basename $n) stage B (box padded to 160): $(cat $n/rows.txt 2>/dev/null | wc -l) rows"
    $PY $T/arm_readout_any.py $T/maps160/toy_consensus.mrc $CL $B/recvol_state0?_even_unfil.mrc 2>&1 | grep -v Warning; done
  F=$H/final_maps; C=$(awk '{print $2}' $F/choice.txt 2>/dev/null); [ -n "$C" ] || continue
  echo "######## $r delivered ($C): $(grep -v '^#' $F/units_index.txt | awk '{printf "%s %s %d | ", $1, $2, $3}')"
  $PY $T/arm_readout_any.py $T/maps160/toy_consensus.mrc $CL $C/unit_*.mrc 2>&1 | grep -v Warning; done
