#!/bin/bash
#SBATCH -J klem_diag_dd
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# diagnostics of the data-derived merge statistic on the zoomed sub-states of RTAG: FSC curves of known pairs (the closed pair = same state)
# under env WHITEN (0/1) and MASK (region/sphere); DIAG_ONLY=1 stops before the merge loop
set -uo pipefail; PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; A=/mnt/beegfs/elmlund/afan/rnasep_klem
export LEAFGLOB=${LEAFGLOB:-zoomregion*} WHITEN=${WHITEN:-0} MASK=${MASK:-region} DIAG_ONLY=1 FLOOR=${FLOOR:-12}
echo "== WHITEN $WHITEN MASK $MASK LEAFGLOB $LEAFGLOB"; $PY $W/merge_data_derived.py $A/${RTAG:-f692c} $W/maps160/toy_consensus.mrc /tmp/diag_$$ 2>&1 | grep -v Warning; echo "DIAG DONE"
