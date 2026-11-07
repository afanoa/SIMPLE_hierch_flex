#!/bin/bash
#SBATCH -J sst_gate4
#SBATCH -p csbdevel
#SBATCH --constraint=x6342
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
# one hetero_tree gate run: a node fixture (FIX=<dir> with input.simple, fixture_rows.txt, fixture_selection.txt; made once from
# SRC=<node dir> when FIX has no manifest yet, the r1211 recipe) or a whole small project (PROJ=<project>, no selection). The run's
# own copy of the project gets the queue settings (update_project), then hetero_tree runs with the dataset keys (DSKEYS) and ARMKEYS.
# The recipe9-lineage split test needs SIMPLE_SPLIT_CONTRAST (the rho floor SIMPLE_SPLIT_RHOMIN is unset in production since 2026-10-29): without it the run refuses to start unless NO_SPLIT_ENV
# names why (a build that reads neither); split_env.txt records them either way. BASEKEYS (default oneflex=yes sta0_npeaks=2 junk_screen=no) can be narrowed.
# usage: sbatch --export=ALL,RUN=..,BUILD=..,DSKEYS=..,ARMKEYS=..,{FIX=..[,SRC=..]|PROJ=..}[,MAXDEPTH=4,NODEPAR=1,SEED=1,MEM=32000] sst_gate.sbatch
set -euo pipefail
: ${RUN:?} ${BUILD:?} ${DSKEYS:?}
if [ -z "${NO_SPLIT_ENV:-}" ]; then : ${SIMPLE_SPLIT_CONTRAST:?the split test needs SIMPLE_SPLIT_CONTRAST}; fi
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
PY=/mnt/nasapps/production/python/3.12/bin/python3; MEM=${MEM:-32000}; PART=${QPART:-csbdevel}
test -x $BUILD/bin/simple_exec || exit 3
if [ -n "${FIX:-}" ] && [ ! -e $FIX/manifest.txt ]; then
  : ${SRC:?}; test ! -e $FIX || { echo "$FIX exists without a manifest"; exit 4; }
  mkdir -p $FIX; cp $SRC/n.simple $FIX/input.simple; cp $SRC/rows.txt $FIX/fixture_rows.txt
  if [ -e $SRC/node_local.txt ]; then awk '{print 1}' $FIX/fixture_rows.txt > $FIX/fixture_selection.txt; mode=row-local
  else
    nrows=$(grep -c . $SRC/../hetero_tree_labels.txt)
    $PY - $SRC $FIX $nrows <<'EOF'
import sys, pathlib
src, fix, n = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2]), int(sys.argv[3]); sel = {int(x) for x in (src / 'rows.txt').read_text().split()}
with open(fix / 'fixture_selection.txt', 'x') as f: f.writelines('1\n' if r in sel else '0\n' for r in range(1, n + 1))
EOF
    mode=root-sized
  fi
  printf '%s\n' "source_node=$SRC" "mode=$mode" "rows=$(grep -c . $FIX/fixture_rows.txt)" "selected=$(awk '$1 > 0' $FIX/fixture_selection.txt | wc -l)" \
    "source_project_sha256=$(sha256sum $SRC/n.simple | cut -c1-64)" "made=$(date -Is)" > $FIX/manifest.txt
fi
test ! -e $RUN || { echo "$RUN exists"; exit 4; }
mkdir -p $RUN; cd $RUN; cp "$0" launcher.sbatch
printf '%s\n' "SIMPLE_SPLIT_CONTRAST=${SIMPLE_SPLIT_CONTRAST:-unset}" "SIMPLE_SPLIT_RHOMIN=${SIMPLE_SPLIT_RHOMIN:-unset}" "SIMPLE_SPLIT_ROUTE=${SIMPLE_SPLIT_ROUTE:-unset}" "SIMPLE_SPLIT_LOCALITY=${SIMPLE_SPLIT_LOCALITY:-unset}" "NO_SPLIT_ENV=${NO_SPLIT_ENV:-}" > split_env.txt
if [ -n "${FIX:-}" ]; then cp $FIX/input.simple input.simple; cp $FIX/fixture_rows.txt fixture_rows.txt; INF="infile=$FIX/fixture_selection.txt"
else cp ${PROJ:?} input.simple; INF=""; fi
export SIMPLE_PATH=$BUILD PATH=$BUILD/bin:$BUILD/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=$PART SIMPLE_EMAIL=nil OMP_NUM_THREADS=8
export SIMPLE_SEED=${SEED:-1} SBATCH_CONSTRAINT=x6342 SBATCH_PARTITION=$PART
X="nthr=8 nodepar=${NODEPAR:-1} maxdepth=${MAXDEPTH:-4} ${BASEKEYS-oneflex=yes sta0_npeaks=2 junk_screen=no} $DSKEYS ${ARMKEYS:-}"
printf '%s\n' "input=${FIX:-$PROJ}" "build=$BUILD" "simple_exec_sha256=$(sha256sum $BUILD/bin/simple_exec | cut -c1-64)" "seed=$SIMPLE_SEED" \
  "keys=$X $INF" "split_env=$(env | grep '^SIMPLE_SPLIT' | tr '\n' ' ' || true)" "rho_floor=${SIMPLE_SPLIT_RHOMIN:-unset (off)}" "job=$SLURM_JOB_ID" "started=$(date -Is)" > run_manifest.txt
simple_exec prg=update_project projfile=input.simple qsys_name=slurm qsys_partition=$PART job_memory_per_task=$MEM walltime=14400 mkdir=no \
  > update_project.log 2>&1 || { echo "update_project FAILED" >> run_manifest.txt; exit 5; }
simple_exec prg=hetero_tree projfile=input.simple $INF $X > tree.log 2>&1 || echo "hetero_tree FAILED" >> run_manifest.txt
echo "finished=$(date -Is)" >> run_manifest.txt
