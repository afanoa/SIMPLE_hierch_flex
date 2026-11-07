#!/bin/bash
#SBATCH -J klem_pole_check
#SBATCH -p csbdevel
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 00:10:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Which pole does a node's kept rule index mean? Per seed dir of the null and RNase P series: site_poles.txt, site_pole_rule.txt, the
# state counts of site_seeds.txt and of the node's seeds/seeds.txt, and the mean of every endpoint posterior over each seed state's rows
set -uo pipefail
S=/mnt/beegfs/elmlund/afan/rnasep_klem/stree_port; IG=$S/fulltoy_m52M_igm/1_hetero_tree; CD=$S/fulltoy_m52M_cd1/1_hetero_tree; M=$S/full692_m33_s1/1_hetero_tree
for e in igm_r212:$IG/node_r212 igm_r2122:$IG/node_r2122 m33_r:$M/node_r m33_r1:$M/node_r1 m33_r2:$M/node_r2 cdisc_r111:$CD/node_r111 cdisc_r1111:$CD/node_r1111 \
         cdisc_r11111:$CD/node_r11111 cdisc_r111111:$CD/node_r111111 cdisc_r1112:$CD/node_r1112 igm_r1111112:$IG/node_r1111112 igm_r111212:$IG/node_r111212 \
         igm_r12121:$IG/node_r12121 igm_r211112:$IG/node_r211112; do nm=${e%%:*}; N=${e#*:}; id=${N##*/node_}
  k=$(head -1 $N/sites/seed_site.txt 2>/dev/null | tr -d ' '); [ -z "$k" ] && k=1
  sd=$N/site; [ $k -gt 1 ] && { sd=$N/nsite_${id}_s$k; [ -d $sd ] || sd=$N/fsite_${id}_s$k; }
  echo "######## $nm: seed site $k, seed dir $(basename $sd)"
  echo "  poles: $(tr '\n' ';' < $sd/site_poles.txt)"
  echo "  rule:  $([ -f $sd/site_pole_rule.txt ] && grep -v '^#' $sd/site_pole_rule.txt | tr '\n' ';' || echo 'no site_pole_rule.txt')"
  echo "  site_seeds states: $(sort $sd/site_seeds.txt | uniq -c | awk '{printf "%s:%s ", $2, $1}')| node seeds states: $(sort $N/seeds/seeds.txt | uniq -c | awk '{printf "%s:%s ", $2, $1}')"
  for ep in z1lo z1hi z3lo z3hi; do f=$sd/site_post_$ep.txt; [ -f $f ] || continue
    echo "  $ep mean posterior by site_seeds state: $(paste $sd/site_seeds.txt $f | awk '$1 > 0 {s[$1] += $2; n[$1]++} END {for (q in s) printf "%s:%.3f ", q, s[q] / n[q]}')"
  done
done
