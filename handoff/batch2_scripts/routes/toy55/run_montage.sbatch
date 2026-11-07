#!/bin/bash
#SBATCH -J toy10_montage
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; export TOY_SMPD=1.464
STATES=1,10,18,27,28,37,46,48,49,51; TRUTH=$(for c in $(echo $STATES | tr , " "); do ls $W/maps160/docked/conf$(printf %02d $c)_*.mrc | head -1; done | tr "\n" , | sed 's/,$//')
for ARM in lat_sup3 lat_sup2 lat; do R=$(ls -d $W/$ARM/maps_r? | sort | tail -1); r=${R##*_r}
  $PY $W/montage_windows.py $R $W/$ARM/windows_r$r.txt $W/maps160/toy_consensus.mrc --truth $TRUTH --states $STATES --out $W/$ARM/fig_windows_r$r.png --title "$ARM round $r: window maps vs dominant truth conformer (60-10 A, z projections)"; done
echo DONE
