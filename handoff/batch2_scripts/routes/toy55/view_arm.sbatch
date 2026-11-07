#!/bin/bash
#SBATCH -J klem_view_arm
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 00:40:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# closed-arm readout of the view_rflex folders: the rflex_f692rf box-88 maps redone at box 160 (Fourier pad), explicit file lists
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem; V=$K_/view_rflex; I=$V/info
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
for nd in r r1; do $PY $T/make_view.py $V/f692rf_node_$nd $I/spec_rf_$nd.txt --box 160; done
for d in f692gen_node39_track f692rf_node_r f692rf_node_r1; do $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $I/arm_$d.txt $(ls $V/$d/*.mrc) 2>&1 | grep -v Warning | tail -1 | sed "s/^/$d: /"
  awk '!/^#/ {n = split($1, p, "/"); printf "   %s closed-arm %.2f open-arm %.2f %s\n", p[n], $2, $3, $4}' $I/arm_$d.txt; done
echo DONE
