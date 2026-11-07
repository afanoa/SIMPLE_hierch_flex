#!/bin/bash
#SBATCH -J klem_calE_null2
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 00:25:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# More nulls for the switch / swing call, calE_null.sbatch's recipe: the IgM M tree's narrow genuine splits r1111112 and r111212
# (passed the re-test by < 1.5x the null max) and weak ones r12121 and r211112; the cdisc M tree's other null-type nodes r1112,
# r11111 and r111111. 5 at-node bins along the seeding score, the node's own project, recB0 line and mask
set -uo pipefail
NP=2; NT=4; TAG=calE_null2; source /mnt/beegfs/elmlund/afan/routes/toy55/seriesE_lib.sh
OUT=$K/factorial/phaseE/null2; mkdir -p $OUT; CD=$S/fulltoy_m52M_cd1/1_hetero_tree; IG=$S/fulltoy_m52M_igm/1_hetero_tree
ST=$K/factorial/junk_test/simsites_truth_jt.txt; IT=$K/factorial/armsigm1_truth/arms_truth.txt
for e in igm_r1111112:$IG:r1111112:arms:$IT igm_r111212:$IG:r111212:arms:$IT igm_r12121:$IG:r12121:arms:$IT igm_r211112:$IG:r211112:arms:$IT \
         cdisc_r1112:$CD:r1112:cdisc:$ST cdisc_r11111:$CD:r11111:cdisc:$ST cdisc_r111111:$CD:r111111:cdisc:$ST; do IFS=: read nm D nd kind tf <<< "$e"
  N=$D/node_$nd; newdir null_$nm
  { $PY $T/null_score.py $N score.txt > ../score_rule.txt 2>&1 && $PY $T/motion_bins4.py score.txt states.txt $kind $tf tree $D > ../bins.txt 2>&1 \
    && recon null_$nm $N/n.simple $N/recB0/hetero_step.sh $N/mask.mrc; } || log "null_$nm FAILED: $(tail -qn1 ../score_rule.txt ../bins.txt 2>&1)"
done
finish
