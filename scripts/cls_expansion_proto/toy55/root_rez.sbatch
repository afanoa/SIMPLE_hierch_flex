#!/bin/bash
#SBATCH -J klem_rez
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# re-decide one data-tree node from its existing maps with the z statistic (env MODE, TAG (dtree1), ND node id (r), NSHUF (3))
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; K=$A/rnasep_klem; MODE=${MODE:?}; TAG=${TAG:-dtree1}; ND=${ND:-r}
NA=$K/synth_${MODE}_A/$TAG/$ND; NB=$K/synth_${MODE}_B/$TAG/$ND; CONS=$W/maps160/toy_consensus.mrc; MV="FAMS=0 MASK=sphere NPROC=16 MINREC=100 FSCMAX=30 NUIS=${NUIS:-2} NULLCAL=0.99 STAT=${STAT:-z} MATCHBAND=${MATCHBAND:-low}"; X=${STAT:-z}${NUIS:-2}${MATCHBAND:-low}
for s in $(seq 1 ${NSHUF:-3}); do env $MV LEAFGLOB=regionshuf$s $PY $W/merge_cross_tree.py $NA $NB $CONS $NA/merge${X}_shuf$s > $NA/merge${X}_shuf$s.log 2>&1; echo "shuf$s: $(grep -a "pair statistics" $NA/merge${X}_shuf$s.log | sed 's/.*sorted): //' | cut -c1-140)"; done
thr=$(cat $NA/merge${X}_shuf*.log | grep -a "pair statistics at the start" | sed 's/.*sorted): //; s/; threshold.*//' | tr ' ' '\n' | grep . | sort -g | tail -1)
env $MV LEAFGLOB=regioncand ARCH=1 THRMIN=$thr $PY $W/merge_cross_tree.py $NA $NB $CONS $NA/merge${X}_cand > $NA/merge${X}_cand.log 2>&1
echo "== $MODE node $ND, statistic $X, shuffled-control threshold $thr"; grep -a "pair statistics\|archetype test\|archetypes:\|merges ->" $NA/merge${X}_cand.log | cut -c1-220
for H in A B; do N=$K/synth_${MODE}_$H/$TAG/$ND; $PY $W/group_labels.py $N/fam0/regioncand/level1/states_full.txt $NA/merge${X}_cand_$H.npz $N/children_$X.txt > /dev/null 2>&1
  echo " half $H:"; $PY $W/dtree_compose.py $N/children_$X.txt $W/sim${MODE}_$H/gt_particles.txt 2>&1 | grep -v Warning
  [ -s $W/truthtree_${MODE%_*}_200k.npz ] && [ "${MODE#55}" != "$MODE" ] && $PY $W/clade_compose.py $N/children_$X.txt $W/sim${MODE}_$H/gt_particles.txt $W/truthtree_55_200k.npz 0.45,0.40 2>&1 | grep -v Warning; done
