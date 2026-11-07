#!/bin/bash
#SBATCH -J klem_arm_census
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:15:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# closed-arm / open-arm readouts (10 A) of the r1 pole census endpoints' top-2.5k maps (factorial/f692st_r1_census/rec_*_top2500)
export OMP_NUM_THREADS=2; PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; O=/mnt/beegfs/elmlund/afan/rnasep_klem/factorial/f692st_r1_census
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
for d in $O/rec_*_top2500; do m=$(ls $d/*_reconstruct3D/recvol_state01_even_unfil.mrc 2>/dev/null | head -1); n=$(basename $d); n=${n#rec_}
  [ -n "$m" ] || { echo "$n: no map"; continue; }; echo "$n $($PY $T/arm_readout_any.py $T/maps160/toy_consensus.mrc $CL $m 2>&1 | grep -v Warning | tail -1 | awk '{print $2, $3, $4, $5}')"; done
