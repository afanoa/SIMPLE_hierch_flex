#!/bin/bash
#SBATCH -J pose_hist
#SBATCH -p csbdevel
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; O=$K/polish_rnasep692_gs1
for h in A B; do echo "== real polish half $h: consensus -> polished"; $PY $W/pose_change_hist.py $W/lat_rnasep692_$h/ptcl3D_r0.txt $O/$h/refined/poses_full.txt $O/refstates_$h.txt 2>&1 | grep -v Warning; done
