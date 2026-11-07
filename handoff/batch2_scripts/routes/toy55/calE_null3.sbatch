#!/bin/bash
#SBATCH -J klem_calE_null3
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 00:15:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# calE_null2.sbatch's nulls whose last-ranked site wrote no site_pole_rule.txt (the fflex fallback seeded them): IgM M r1111112,
# r111212, r12121, r211112 and cdisc M r1112, r111111, scored with null_score2.py; the same at-node recipe otherwise
set -uo pipefail
NP=2; NT=4; TAG=calE_null3; source /mnt/beegfs/elmlund/afan/routes/toy55/seriesE_lib.sh
OUT=$K/factorial/phaseE/null3; mkdir -p $OUT; CD=$S/fulltoy_m52M_cd1/1_hetero_tree; IG=$S/fulltoy_m52M_igm/1_hetero_tree
ST=$K/factorial/junk_test/simsites_truth_jt.txt; IT=$K/factorial/armsigm1_truth/arms_truth.txt
for e in igm_r1111112:$IG:r1111112:arms:$IT igm_r111212:$IG:r111212:arms:$IT igm_r12121:$IG:r12121:arms:$IT igm_r211112:$IG:r211112:arms:$IT \
         cdisc_r1112:$CD:r1112:cdisc:$ST cdisc_r111111:$CD:r111111:cdisc:$ST; do IFS=: read nm D nd kind tf <<< "$e"
  N=$D/node_$nd; newdir null_$nm
  { $PY $T/null_score2.py $N score.txt > ../score_rule.txt 2>&1 && $PY $T/motion_bins4.py score.txt states.txt $kind $tf tree $D > ../bins.txt 2>&1 \
    && recon null_$nm $N/n.simple $N/recB0/hetero_step.sh $N/mask.mrc; } || log "null_$nm FAILED: $(tail -qn1 ../score_rule.txt ../bins.txt 2>&1)"
done
finish
