#!/bin/bash
#SBATCH -J rnasep_f692_family
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 20:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# One family of the f692 chain as its own job (so the families run concurrently): child board, k-means,
# refine3D_states, feedback board (f692_family.sh), then the parallel regions stage (f692_regions_par2.sh). Writes
# $R/fam$FAM/DONE. env: R, FAM, SRC, METAALL, FB, NTHR, NP, MAXITS, LPSC, SIGC, MAXREG, NL2, PARMAX
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; log(){ echo "[$(date '+%F %T') famjob$FAM] $*"; }
export R SRC FB NTHR NP MAXITS LPSC SIGC NL2 PARMAX NPARTS SERVER=${SERVER:-0} LOCALDEV=${LOCALDEV:-0}
[ "$LOCALDEV" = 1 ] && { export SIMPLE_QSYS=local; log "LOCALDEV: everything on this node (${SLURM_CPUS_ON_NODE:-?} cores allocated): ${NPARTS:-0} local workers per board, server $SERVER, no queue"; }
stopsrv(){ SL=$R/fam$FAM/server_boards.txt; [ -s $SL.pid ] && kill -0 $(cat $SL.pid) 2>/dev/null && { echo STOP >> $SL; log "board server stopped"; }; }
log "family $FAM stage"; FAM=$FAM META=$METAALL bash $W/f692_family.sh || { log "family stage FAILED"; stopsrv; exit 1; }
[ -s $R/fam$FAM/klem_fb/klem_latent.txt ] || { log "no feedback board"; exit 1; }
SL=$R/fam$FAM/server_boards.txt; if [ "$SERVER" = 1 ] && [ -s $SL.pid ] && kill -0 $(cat $SL.pid) 2>/dev/null; then export SERVER_LIST=$SL; log "board server alive (pid $(cat $SL.pid)): region boards are served"; else unset SERVER_LIST; fi
RS=$W/f692_regions_par3.sh; [ "$LOCALDEV" = 1 ] && RS=$W/f692_regions_local.sh
log "regions stage, family $FAM"; FAM=$FAM META=$R/fam$FAM/meta_fb.npz MAXREG=${MAXREG:-6} NL2=${NL2:-3} bash $RS || log "regions stage FAILED"
[ -n "${SERVER_LIST:-}" ] && { echo STOP >> $SERVER_LIST; while kill -0 $(cat $SERVER_LIST.pid) 2>/dev/null; do sleep 10; done; log "board server finished"; }
D=$R/fam$FAM; if [ -L $D/refine ]; then RW=$(readlink $D/refine); mkdir -p $D/refine_shared; for f in $(cd $RW && find . -maxdepth 2 \( -name 'p.simple' -o -name 'recvol_state*.mrc' -o -name '*.log' -o -name 'fsc_*.bin' -o -name 'FAILED' -o -name 'up.log' \)); do mkdir -p $D/refine_shared/$(dirname $f); cp -f $RW/$f $D/refine_shared/$f; done; rm -f $D/refine; mv $D/refine_shared $D/refine; log "refinement results copied back from $RW"; fi
touch $R/fam$FAM/DONE; log DONE
