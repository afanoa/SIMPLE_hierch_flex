#!/bin/bash
#SBATCH -J klem_p3masks2
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=110G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# focus masks in matched units (problem 3, p3_masks2.py: % of sphere and % of molecule, member + split masks, true supports): synthetic panels and
# the realistic-SNR recursive run's flex states (rflex_s30a, with truth) against the real f692r families and rflex_f692rf flex states (read only).
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem; O=$K_/p3masks; mkdir -p $O
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2; ARM=$T/maps160/toy_consensus.mrc:$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
run(){ $PY $T/p3_masks2.py $O "$@" --reg 0 > $O/log2_$2.txt 2>&1 & }
for nd in r r1 r2 r21; do run --name f692rf_$nd --list $K_/rflex_f692rf/node_$nd/members.txt --mskdiam 232 --arm $ARM; done
for nd in r r1 r2 r11 r12 r21 r22; do run --name s30a_$nd --list $K_/rflex_s30a/node_$nd/members.txt --labels $K_/rflex_s30a/node_$nd/flex_states.txt --gt $T/simswing30/gt_particles.txt --truth $T/simswing30/maps --mskdiam 220; done
for t in swing30hs prodihs; do P=$K_/panel_$t; run --name panel_$t --rdir $P/1_refine3D_states --labels $P/refined_states.txt --tree $P/smerge/merge_tree.txt --gt $T/sim$t/gt_particles.txt --truth $T/sim$t/maps --mskdiam 220 --minleaves 2; done
for F in 1 2; do C=$K_/f692r/coordtree_fam$F; run --name f692r_fam$F --rdir $(ls -td $C/rec/*_reconstruct3D | head -1) --labels $C/units.txt --tree $C/tree/merge_tree.txt --mskdiam 232 --arm $ARM --minleaves $([ $F = 1 ] && echo 3 || echo 8); done
wait; for f in $O/*_masks2.txt; do cat $f; echo; done; grep -l "Traceback\|Error" $O/log2_*.txt && echo "FAILURES above"; echo DONE
