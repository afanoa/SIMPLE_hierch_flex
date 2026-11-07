#!/bin/bash
#SBATCH -J klem_bsample
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# even sampling of one branch of a consolidation: branch_axes.py (axes of the state maps in the family periphery, G x G grid over axes 1-2)
# -> mu_rec.sbatch (each sample reconstructed from all its particles at the family-refined poses) -> halves into maps/halves -> closed-arm
# readout of the sample maps. env RTAG (f692r), C (consolidation dir under RTAG), STATES (comma list), FAM (1), OUT (dir under RTAG), G (3)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; R=/mnt/beegfs/elmlund/afan/rnasep_klem/${RTAG:?}
C=$R/${C:?}; O=$R/${OUT:?}; G=${G:-3}; log(){ echo "[$(date '+%F %T') branch-sample $RTAG $OUT] $*"; }
$PY $W/branch_axes.py $O --rdir $C/fam0/regionall/level1/maps/halves --states ${STATES:?} --groups $C/groups.txt --labels $C/states.txt \
  --mask $R/mergeup/periphery_fam${FAM:-1}.mrc --arm $C/arm_states.txt --k 5 --axis 1 --grid $G 2>&1 | grep -v "Warning\|stddev"
log "reconstructing the grid samples"; RTAG=$RTAG CDIR=$OUT/grid12_g$G CUT=${OUT##*/}_g$G bash $W/mu_rec.sbatch 2>&1 | grep -v "^  *[0-9]* / \|Traceback\|^  File\|^ *\^" | tail -4
bash $W/tidy_halves.sh $O/grid12_g$G
CL=$(ls $W/maps160/docked/conf18_*_m160.mrc | head -1); $PY $W/arm_readout.py $W/maps160/toy_consensus.mrc $CL $O/grid12_g$G/arm_states.txt $O/grid12_g$G/fam0/regionall/level1/states_dir 2>&1 | grep -v Warning | tail -1
awk '!/^#/ {n = split($1, p, "/"); print "   " p[n], "closed-arm", $2, "open-arm", $3, $4}' $O/grid12_g$G/arm_states.txt; log DONE
