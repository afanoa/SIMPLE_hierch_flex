#!/bin/bash
#SBATCH -J klem_node_readout
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 00:40:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# closed-state readout of one rflex node's refined states: closed-arm / open-arm occupancy of every state map (lp 10, box 160) and the
# old f692r closed states 15+16 (closed-arm 0.87 / 0.83) traced into each state (count, share of the state, share of the node's closed).
# env RT (rflex_<RT>), NODE (r), STAGE (fstageB | fstageA | softrec | focus | oldfocus)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; RT=${RT:?}; NODE=${NODE:-r}; ST=${STAGE:-fstageB}
D=$K_/rflex_$RT/node_$NODE; RD=$(ls -td $D/$ST/*_refine3D $D/$ST/*_reconstruct3D 2>/dev/null | head -1); LAB=$D/refined_f.txt; [ $ST = softrec ] && LAB=$D/refined_s.txt
[ $ST = focus ] && { RD=$(ls -td $D/focus/*_flex_pca | head -1); LAB=$D/fflex_states_raw.txt; }   # focused flex states (raw labels = the map numbers)
[ $ST = oldfocus ] && { RD=$(ls -td $D/old_flexmask/focus/*_flex_pca | head -1); LAB=$D/old_flexmask/fflex_states_raw.txt; }   # the set-aside focused flex
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1); U=$K_/f692r/coordtree_fam1/units.txt; V=$K_/view_rflex/${RT}_${NODE}_$ST; O=$D/readout_$ST; mkdir -p $O
NS=$(awk '$1 > 0' $LAB | sort -n | tail -1); : > $O/spec.txt
for s in $(seq 1 $NS); do b=$RD/recvol_state$(printf %02d $s); [ -s ${b}_even_unfil.mrc ] && echo "state_$(printf %02d $s) ${b}_even_unfil.mrc ${b}_odd_unfil.mrc" >> $O/spec.txt; done
$PY $T/make_view.py $V $O/spec.txt --box 160 > /dev/null 2>&1
$PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/arm.txt $V 2>&1 | grep -v Warning | tail -1
echo "node $NODE ($ST, $RT): state particles closed-arm open-arm | old closed 15+16 in the state (share of the state, share of the node's closed)"
awk 'FILENAME == ARGV[1] {c[FNR] = ($1 == 15 || $1 == 16); next} FILENAME == ARGV[2] {if ($1 > 0) {n[$1]++; k[$1] += c[FNR]; tk += c[FNR]}; next}
  !/^#/ {s = $1; sub(/.*state_/, "", s); sub(/\.mrc/, "", s); s += 0; printf "   %2d %6d  %.2f %.2f | %5d (%4.1f %%, %4.1f %%)\n", s, n[s], $2, $3, k[s], 100 * k[s] / (n[s] > 0 ? n[s] : 1), 100 * k[s] / (tk > 0 ? tk : 1)}' $U $LAB $O/arm.txt | sort -k3,3gr
