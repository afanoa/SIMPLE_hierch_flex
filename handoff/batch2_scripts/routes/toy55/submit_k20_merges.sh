#!/bin/bash
# bottom-up trees of the K=20 composition runs (complete and average linkage), truth readout per cut
cd /mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem
for t in ${TAGS:-swing30hs swing30 prodihs}; do R=$K/ff2_$t/comp_k20; D=$(ls -td $R/*_refine3D_states 2>/dev/null | head -1)
  [ -s $R/refined_states.txt ] && [ -n "$D" ] || { echo "$t: not finished"; continue; }
  for l in complete average; do NCUT=16,12,8,6,4,2 sbatch --parsable --export=ALL,RD=$D,LABELS=$R/refined_states.txt,OUT=$R/smerge_$l,GT=/mnt/beegfs/elmlund/afan/routes/toy55/sim$t/gt_particles.txt,LINK=$l state_merge.sbatch; done
done
