#!/bin/bash
#SBATCH -J klem_phaseC
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# phase B series re-tested: amplitude matching, focus mask and the per-shell mixture fit (phaseC_maptests.py) on the 6 series
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; S=$K/stree_port
P=$K/factorial/phaseB; I=$K/factorial/motion_series_igm; CA=$S/arm_M_c11_s1/1_hetero_tree; CD=$S/fulltoy_m52M_cd1/1_hetero_tree
MI=$S/fulltoy_m52M_igm/1_hetero_tree/node_r/mask.mrc; MB=$S/fulltoy_igm/1_hetero_tree/node_r/mask.mrc
$PY $T/phaseC_maptests.py $P/cahra_open:$CA/node_r/mask.mrc:160:open $P/cdisc_A:$CD/node_r/mask.mrc:225.456:A $P/cdisc_closed:$CD/node_r/mask.mrc:225.456:closed \
  $I/m52M_r_s1_a4th_q5:$MI:220:mean $I/m52M_r2_s4_a2th_q5:$MI:220:mean $I/m8c17_r2_s1_a2th_q5:$MB:220:mean 2>&1 | grep -v Warning | tee $P/maptests_C.txt
