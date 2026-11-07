#!/bin/bash
#SBATCH -J klem_calE_null
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 00:12:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# The null case of the switch / swing call: at four of the toys' false splits (cdisc M r111 and r1111, fflex fallback; IgM M r212
# and r2122, no truth support) 5 at-node bins along the score that seeded the split (null_score.py; motion_bins4.py, tree junk
# out), the node's own project and recB0 reconstruct3D with nstates=5, postprocess=no; phaseE_maptests.py in the node's mask
set -uo pipefail
NP=2; NT=4; TAG=calE_null; source /mnt/beegfs/elmlund/afan/routes/toy55/seriesE_lib.sh
OUT=$K/factorial/phaseE/null; mkdir -p $OUT; CD=$S/fulltoy_m52M_cd1/1_hetero_tree; IG=$S/fulltoy_m52M_igm/1_hetero_tree
ST=$K/factorial/junk_test/simsites_truth_jt.txt; IT=$K/factorial/armsigm1_truth/arms_truth.txt
for e in cdisc_r111:$CD:r111:cdisc:$ST cdisc_r1111:$CD:r1111:cdisc:$ST igm_r212:$IG:r212:arms:$IT igm_r2122:$IG:r2122:arms:$IT; do IFS=: read nm D nd kind tf <<< "$e"
  N=$D/node_$nd; newdir null_$nm
  { $PY $T/null_score.py $N score.txt > ../score_rule.txt 2>&1 && $PY $T/motion_bins4.py score.txt states.txt $kind $tf tree $D > ../bins.txt 2>&1 \
    && recon null_$nm $N/n.simple $N/recB0/hetero_step.sh $N/mask.mrc; } || log "null_$nm FAILED: $(tail -qn1 ../score_rule.txt ../bins.txt 2>&1)"
done
finish
