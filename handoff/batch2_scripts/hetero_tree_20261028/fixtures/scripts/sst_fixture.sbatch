#!/bin/bash
#SBATCH -J sst_fixture
#SBATCH -p norm
#SBATCH --constraint=x6342
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
# closed-subtree fixture from a finished full692 tree node (the r1211 recipe: the node's own project and rows, no new project):
# FIX=<fixture dir> is made once from SRC=<node dir> (input.simple, fixture_rows.txt, fixture_selection.txt, manifest.txt) and never
# written again; each run copies input.simple into RUN=<run dir>. A row-local node project (node_local.txt) gets an all-ones
# selection, and fixture_rows.txt maps its rows to the 932,496 project rows; a root-sized one gets the 932,496-row 0/1 selection.
# usage: sbatch --export=ALL,SRC=..,FIX=..,RUN=..,BUILD=..,MAXDEPTH=..,ARMKEYS='site_layer=no' sst_fixture.sbatch
set -euo pipefail
: ${SRC:?} ${FIX:?} ${RUN:?} ${BUILD:?} ${MAXDEPTH:?} ${ARMKEYS:?}
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
PY=/mnt/nasapps/production/python/3.12/bin/python3
test -x $BUILD/bin/simple_exec || exit 3
if [ ! -e $FIX/manifest.txt ]; then
  test ! -e $FIX || { echo "$FIX exists without a manifest"; exit 4; }
  mkdir -p $FIX
  cp $SRC/n.simple $FIX/input.simple; cp $SRC/rows.txt $FIX/fixture_rows.txt
  if [ -e $SRC/node_local.txt ]; then
    awk '{print 1}' $FIX/fixture_rows.txt > $FIX/fixture_selection.txt; mode=row-local
  else
    $PY - $SRC $FIX <<'EOF'
import sys, pathlib
src, fix = map(pathlib.Path, sys.argv[1:3]); sel = {int(x) for x in (src / 'rows.txt').read_text().split()}
with open(fix / 'fixture_selection.txt', 'x') as f: f.writelines('1\n' if r in sel else '0\n' for r in range(1, 932497))
EOF
    mode=root-sized
  fi
  printf '%s\n' "source_node=$SRC" "mode=$mode" "rows=$(grep -c . $FIX/fixture_rows.txt)" "selected=$(awk '$1 > 0' $FIX/fixture_selection.txt | wc -l)" \
    "source_project_sha256=$(sha256sum $SRC/n.simple | cut -c1-64)" "made=$(date -Is)" > $FIX/manifest.txt
fi
test ! -e $RUN || { echo "$RUN exists"; exit 4; }
mkdir -p $RUN; cd $RUN
cp $FIX/input.simple input.simple; cp $FIX/fixture_rows.txt fixture_rows.txt; cp "$0" launcher.sbatch
export SIMPLE_PATH=$BUILD PATH=$BUILD/bin:$BUILD/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=${QPART:-norm,csbdevel} SIMPLE_EMAIL=nil OMP_NUM_THREADS=8
export SIMPLE_SEED=${SEED:-1} SBATCH_CONSTRAINT=x6342 SBATCH_PARTITION=${QPART:-norm,csbdevel}
X="mskdiam=232 nparts=4 nthr=8 nodepar=${NODEPAR:-3} maxdepth=$MAXDEPTH minpop=3000 oneflex=yes sta0_npeaks=2 junk_screen=no $ARMKEYS"
printf '%s\n' "fixture=$FIX" "build=$BUILD" "simple_exec_sha256=$(sha256sum $BUILD/bin/simple_exec | cut -c1-64)" "seed=$SIMPLE_SEED" \
  "keys=$X" "split_env=$(env | grep '^SIMPLE_SPLIT' | tr '\n' ' ' || true)" "job=$SLURM_JOB_ID" "started=$(date -Is)" > run_manifest.txt
simple_exec prg=hetero_tree projfile=input.simple infile=$FIX/fixture_selection.txt $X > tree.log 2>&1 || echo "hetero_tree FAILED" >> run_manifest.txt
echo "finished=$(date -Is)" >> run_manifest.txt
