#!/bin/bash
#SBATCH -J nodetype_setup
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Node-type dev set: junk-free ("clean") pseudo-nodes rc of the null_50k and two_50k roots (members = truth non-junk; poses, mean and
# flex latent linked from the root node r), then cheap children (cheapkids2.sbatch) for r and rc of both sets.
set -uo pipefail
A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem
for m in null_50k two_50k; do R=$K_/rflex_root_$m
  awk 'NR == FNR {if ($0 !~ /^#/) c[$1] = $2; next} {r++; print ($1 > 0 && c[r] > 0) ? 1 : 0}' $T/sim$m/gt_particles.txt $R/rows_r.txt > $R/rows_rc.txt
  PP=$(awk '$1 == "r" {print $3}' $R/queue.txt); grep -q "^rc " $R/queue.txt || echo "rc 1 $PP $R/rows_rc.txt" >> $R/queue.txt
  mkdir -p $R/node_rc; ln -sfn $R/node_r/mean $R/node_rc/mean; ln -sfn $(ls -td $R/node_r/*_flex_pca | head -1) $R/node_rc/1_flex_pca
  echo "$m: clean node rc = $(awk '$1 > 0' $R/rows_rc.txt | wc -l) of $(awk '$1 > 0' $R/rows_r.txt | wc -l) members"
  for u in r rc; do sbatch --export=ALL,RUN=rflex_root_$m,NODE=$u,TAG=$m $T/cheapkids2.sbatch; done
done
