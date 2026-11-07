#!/bin/bash
#SBATCH -J klem_boardvar
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# rerun a chain region board with overrides, locally on 16 threads (env RUN = tree dir, FAM, REG e.g. region_1, M = the
# synthetic set tag, TAG = variant name, LPSO / SIGO = new lps / sigmas lines, PS = project override (default: the
# family's refine3D_states project), SPECMATCH, FB (default the run's frozen klem2_bin)); then the arm-latent readout
# (arms sets) of the chain's board and of the variant
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; RUN=${RUN:?}; FAM=${FAM:?}; REG=${REG:?}; M=${M:?}; TAG=${TAG:?}
FB=${FB:-$RUN/klem2_bin}; export SIMPLE_PATH=$FB PATH=$FB/bin:$FB/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16
B0=$RUN/fam$FAM/$REG/board; O=$RUN/fam$FAM/${REG}_$TAG/board; mkdir -p $O; log(){ echo "[$(date '+%F %T') boardvar $REG $TAG] $*"; }
PS=${PS:-$(ls -td $RUN/fam$FAM/refine/*_refine3D_states | head -1)/p.simple}
grep -v "^lps\|^sigmas\|^specmatch" $B0/board.txt > $O/board.txt
echo "lps ${LPSO:-$(grep "^lps" $B0/board.txt | cut -d' ' -f2-)}" >> $O/board.txt; echo "sigmas ${SIGO:-$(grep "^sigmas" $B0/board.txt | cut -d' ' -f2-)}" >> $O/board.txt
[ -n "${SPECMATCH:-}" ] && echo "specmatch $SPECMATCH" >> $O/board.txt
log "project $PS, build $FB"; diff $B0/board.txt $O/board.txt
cd $O; t0=$(date +%s)
[ -s $O/klem_latent.txt ] || simple_exec prg=klem_board projfile=$PS infile=$O/board.txt dir=$O nthr=16 nparts=1 mkdir=no > $O/board.out 2>&1 || { log "board FAILED: $(grep -a -m3 'ERROR\|STOP' $O/board.out | cut -c1-200)"; exit 1; }
log "board done in $(( $(date +%s) - t0 )) s"
if [ -s $W/truth_$M.txt ]; then
  $PY $W/arm_latent_readout.py $B0/klem_latent.txt $W/sim$M/gt_particles.txt $W/truth_$M.txt chain 2>&1 | grep -v Warning
  $PY $W/arm_latent_readout.py $O/klem_latent.txt $W/sim$M/gt_particles.txt $W/truth_$M.txt $TAG 2>&1 | grep -v Warning
fi
log DONE
