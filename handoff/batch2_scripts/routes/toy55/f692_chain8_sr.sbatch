#!/bin/bash
#SBATCH -J rnasep_f692
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# The unsupervised klem tree on the full 692k apo set with the Fortran board (prg=klem_board), one job, stages resume
# from existing outputs: (1) parent board at the consensus poses (7x7 + 2 free, 20 -> 10 A, latent init), (2) families
# by spectral bisection, parent regions, (3) per family: child board, k-means, refine3D_states, feedback board,
# (4) per family: truth-free region boards, states, level 2, (5) landscape count across the leaves.
# env: RTAG (run name), SMOKE (1 = 50k subset, short schedules), NFAM (2), NP (refine3D_states parts, 10; 20 at night), NPARTS (board workers, 10), NTHR (threads per worker, 8), PARMAX (concurrent region jobs, 3)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/${RTAG:-f692}; mkdir -p $R; cd $R
SRC=${SRC:-$A/rnasep_coupled/full692k_2d/rnasep.simple}; METAALL=${METAALL:-$A/rnasep_coupled/basis_span_all/meta_all.npz}; NFAM=${NFAM:-2}; log(){ echo "[$(date '+%F %T') f692] $*"; }
export R SRC
# frozen copy of the Fortran build so the worktree can be rebuilt while this runs
if [ ! -x $R/klem2_bin/bin/simple_exec ]; then mkdir -p $R/klem2_bin && cp -r $HOME/SIMPLE_klem2/build/bin $R/klem2_bin/ && { [ -d $HOME/SIMPLE_klem2/build/scripts ] && cp -r $HOME/SIMPLE_klem2/build/scripts $R/klem2_bin/ || true; }; fi
export FB=$R/klem2_bin NTHR=${NTHR:-8} NP=${NP:-10} PARMAX=${PARMAX:-3} NPARTS=${NPARTS:-10}   # boards distributed over NPARTS workers x NTHR threads
if [ "${SMOKE:-0}" = 1 ]; then LPSP="20 14 10"; SIGP="1.5 1.0 0.7"; LPSC="14 10 8"; SIGC="1.2 0.8 0.7"; MAXITS=1; MAXREG=1; NL2=1; NSUB=50000; else LPSP="20 20 20 14 14 14 10 10 10"; SIGP="1.5 1.5 1.2 1.0 1.0 0.8 0.7 0.7 0.7"; LPSC="14 14 12 12 10 10 8 8 8"; SIGC="1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7"; MAXITS=6; MAXREG=6; NL2=3; NSUB=0; fi
export LPSC SIGC MAXITS NL2
# ---- (1) parent board
$PY - <<PYEOF
import numpy as np
M = np.load('$METAALL'); f = M['flags'].astype(int); rows = np.where(f == 1)[0] + 1; nsub = $NSUB
if nsub > 0: rows = rows[np.random.default_rng(0).choice(len(rows), nsub, replace=False)]; rows.sort()
np.savetxt('$R/subset_all.txt', rows, fmt='%d'); print('parent subset: %d particles (%d flagged)' % (len(rows), (f == 1).sum()))
PYEOF
if [ ! -s $R/parent/klem_latent.txt ]; then NFREE=2 META=$METAALL bash $W/fboard.sh $R/parent $SRC $R/subset_all.txt "$LPSP" "$SIGP" || { log "parent board FAILED"; exit 1; }; fi
echo "== parent board vs the deposited conformers (readout)"; $PY $W/deposited_match.py $R/parent/cells $W/maps160/docked $W/maps160/toy_consensus.mrc 2>&1 | grep -v Warning | tail -12
# ---- (2) families and the parent's own moving region
[ -s $R/expand/family1_subset.txt ] || $PY $W/expand_families.py split $R/parent $R/expand --nfam $NFAM 2>&1 | grep -v Warning | tail -4
wc -l $R/expand/family?_subset.txt
[ -s $R/regions_parent/region_1.mrc ] || $PY $W/derive_regions.py $R/parent $R/regions_parent --maxreg 2 2>&1 | grep -v Warning | tail -4
# ---- (3+4) the families as concurrent jobs: child, refinement, feedback, then the parallel regions stage
for j in $(seq 1 $NFAM); do
  [ -f $R/fam$j/DONE ] && { log "family $j already done"; continue; }
  jid=$(sbatch --parsable -c ${NTHR:-8} --export=ALL,R=$R,FAM=$j,SRC=$SRC,METAALL=$METAALL,FB=$FB,NTHR=${NTHR:-8},NPARTS=${NPARTS:-10},NP=${NP:-10},MAXITS=$MAXITS,LPSC="$LPSC",SIGC="$SIGC",MAXREG=$MAXREG,NL2=$NL2,PARMAX=${PARMAX:-3} -J rnasep_f692_fam$j $W/f692_family_job_sr2.sbatch); log "family $j -> job $jid"
done
n=0; while [ $n -lt 1200 ]; do left=0; for j in $(seq 1 $NFAM); do [ -f $R/fam$j/DONE ] || left=$((left+1)); done; [ $left -eq 0 ] && break; sleep 60; n=$((n+1)); done; log "families finished"
# ---- (5) the landscape: distinct reproducible states across the leaves
DIRS=""; for d in $R/fam?/klem_fb $R/fam?/region*/board; do [ -s $d/klem_halves.npz ] && DIRS="$DIRS $d"; done
for d in $DIRS; do [ -s $d/cell_matches.txt ] || $PY $W/deposited_match.py $d/cells $W/maps160/docked $W/maps160/toy_consensus.mrc --save $d/cell_matches.txt > /dev/null 2>&1; done
echo "== LANDSCAPE over the leaves: $DIRS"; $PY $W/count_states.py $DIRS 2>&1 | grep -v Warning | tail -20
log "ALL DONE ($R)"
