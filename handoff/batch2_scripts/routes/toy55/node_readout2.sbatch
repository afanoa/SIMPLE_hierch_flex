#!/bin/bash
#SBATCH -J klem_node_readout2
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 00:40:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# node_readout.sbatch with the inputs given instead of taken from the node layout (also run with bash inside other jobs): closed-arm / open-arm
# occupancy of every state map (make_view.py lp 10, box 160) and the old f692r closed states 15+16 (coordtree_fam1/units.txt) traced into each state.
# env RD (dir with recvol_stateNN_{even,odd}_unfil.mrc), LAB (per master-row state labels), OUT (readout dir), VIEW (folder of the full viewing maps), NAME
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; RD=${RD:?}; LAB=${LAB:?}; O=${OUT:?}; V=${VIEW:?}
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1); U=$K_/f692r/coordtree_fam1/units.txt; mkdir -p $O
[ $(wc -l < $LAB) = $(wc -l < $U) ] || { echo "${NAME:-readout}: $LAB has $(wc -l < $LAB) rows, $U $(wc -l < $U): not the master rows"; exit 1; }
NS=$(awk '$1 > 0' $LAB | sort -n | tail -1); : > $O/spec.txt
for s in $(seq 1 $NS); do b=$RD/recvol_state$(printf %02d $s); [ -s ${b}_even_unfil.mrc ] && echo "state_$(printf %02d $s) ${b}_even_unfil.mrc ${b}_odd_unfil.mrc" >> $O/spec.txt; done
$PY $T/make_view.py $V $O/spec.txt --box 160 > /dev/null 2>&1
$PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/arm.txt $V 2>&1 | grep -v Warning | tail -1
echo "${NAME:-readout}: state particles closed-arm open-arm | old closed 15+16 in the state (share of the state, share of the node's closed)"
awk 'FILENAME == ARGV[1] {c[FNR] = ($1 == 15 || $1 == 16); next} FILENAME == ARGV[2] {if ($1 > 0) {n[$1]++; k[$1] += c[FNR]; tk += c[FNR]}; next}
  !/^#/ {s = $1; sub(/.*state_/, "", s); sub(/\.mrc/, "", s); s += 0; printf "   %2d %6d  %.2f %.2f | %5d (%4.1f %%, %4.1f %%)\n", s, n[s], $2, $3, k[s], 100 * k[s] / (n[s] > 0 ? n[s] : 1), 100 * k[s] / (tk > 0 ? tk : 1)}' $U $LAB $O/arm.txt | sort -k3,3gr | tee $O/states.txt
