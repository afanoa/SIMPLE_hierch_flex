#!/bin/bash
#SBATCH -J igm_mstart
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%A_%a.out
# one start of one region board (array task t -> region (t-1)/NSTART + 1, start (t-1)%NSTART + 1): balanced NLAT x NLAT
# lattice + 1 free, focus = RUN/regions/region_<r>.mrc, random start (seed 100 r + s), cleaned subset SUB, project PS.
# env RUN, SUB, PS, NSTART (4), NLAT (5)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; export FB=$A/rnasep_klem/f692e/klem2_bin SIMPLE_QSYS=local
t=${SLURM_ARRAY_TASK_ID:?}; NS=${NSTART:-4}; r=$(( (t - 1) / NS + 1 )); s=$(( (t - 1) % NS + 1 )); NL=${NLAT:-5}
D=$RUN/arm_$r; mkdir -p $D; cd $D
[ -s $D/board_s$s/klem_latent.txt ] || SERVER_LIST= NPARTS= NTHR=16 N1=$NL N2=$NL NFREE=1 CTFABS=1 PSCALE=1 FOCUSMASK=$RUN/regions/region_$r.mrc INIT=random SEED=$((100 * r + s)) bash $W/fboard.sh $D/board_s$s $PS $SUB "20 14 12 10 10 8 8 8" "1.5 1.2 1.0 1.0 0.8 0.8 0.7 0.7" > $D/board_s$s.out 2>&1 || { echo "region $r start $s FAILED"; exit 1; }
echo "region $r start $s done"
