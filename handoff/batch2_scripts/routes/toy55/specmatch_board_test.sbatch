#!/bin/bash
#SBATCH -J klem_specmatch_test
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# A/B of the 'specmatch 1' board option on an existing region board (env RUN = tree dir, FAM, REG = region tag e.g.
# region_1 / regionmotion_1, M = the synthetic set tag for the truth): the same board.txt with specmatch added, run with
# the current ~/SIMPLE_klem2/build on 10 workers, then the nuisance readout (cells_nuisance.py) of both boards.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; RUN=${RUN:?}; FAM=${FAM:?}; REG=${REG:?}; M=${M:?}
FB=${FB:-$HOME/SIMPLE_klem2/build}; export SIMPLE_PATH=$FB PATH=$FB/bin:$FB/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8
B0=$RUN/fam$FAM/$REG/board; O=$RUN/fam$FAM/${REG}_specmatch/board; mkdir -p $O; log(){ echo "[$(date '+%F %T') specmatch $REG] $*"; }
PS=$(ls -td $RUN/fam$FAM/refine/*_refine3D_states | head -1)/p.simple; [ -f $PS ] || PS=$RUN/src.simple
grep -v "^specmatch" $B0/board.txt > $O/board.txt; echo "specmatch ${SPECMATCH:-1}" >> $O/board.txt
log "baseline $B0 -> $O (project $PS, build $FB)"; diff $B0/board.txt $O/board.txt
cd $O; t0=$(date +%s); rm -f $O/BOARD_DONE
[ -s $O/klem_latent.txt ] || simple_exec prg=klem_board projfile=$PS infile=$O/board.txt dir=$O nthr=8 nparts=${NPARTS:-10} mkdir=no > $O/board.out 2>&1 || { log "board FAILED: $(grep -a -m3 'ERROR\|STOP' $O/board.out | cut -c1-200)"; exit 1; }
log "board done in $(( $(date +%s) - t0 )) s"; grep -a "timing\|argmax changed\|reproducib" $O/board.out | tail -6 | cut -c1-160
TR=$W/sim$M/all.txt; RF=$W/lat_$M/ptcl3D_r0.txt; GT=$W/sim$M/gt_particles.txt
echo "== nuisance readout, BASELINE board"; $PY $W/cells_nuisance.py $B0/klem_latent.txt $TR $RF $GT 2>&1 | grep -v Warning
echo "== nuisance readout, SPECMATCH board"; $PY $W/cells_nuisance.py $O/klem_latent.txt $TR $RF $GT 2>&1 | grep -v Warning
log DONE
