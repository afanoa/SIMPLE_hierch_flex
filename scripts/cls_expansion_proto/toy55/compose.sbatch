#!/bin/bash
#SBATCH -J klem_compose
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# composed state maps: parent core (692k consensus halves) outside the state's region, the state's own map inside (env RTAG, CDIR)
set -uo pipefail; PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/${RTAG:?}; C=$R/${CDIR:?}; O=$C/composed
echo "== $RTAG/$CDIR composed with the 692k consensus core"; $PY $W/compose_state_maps.py $R $C $A/rnasep_states/25_refine3D_states/rec_final_state01 $O 2>&1 | grep -v Warning; echo "COMPOSE DONE"
