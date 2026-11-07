#!/bin/bash
# launch_gates.sh BUILD ARM [DS ...]: the gate fixtures (qtoy2 md1, r131 md1, qtoynull / olenull / micnull md0) on one build, maps
# kept for site_rank replays; runs land in math_reworks/gates/<ARM>_<ds>. Keys = the recipe18d gate runs' (the latest accepted arms).
set -euo pipefail
BUILD=$1; ARM=$2; shift 2; DS=${@:-qtoy2 r131 qtoynull olenull micnull ole}
M=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks; F=/mnt/beegfs/elmlund/afan/hetero_tree_20261028
test -x $BUILD/bin/simple_exec || { echo "no simple_exec in $BUILD"; exit 3; }
mkdir -p $M/gates $M/logs
export SIMPLE_SPLIT_CONTRAST=maps SIMPLE_TREE_KEEPMAPS=yes
for ds in $DS; do
  RUN=$M/gates/${ARM}_${ds}${RUNSFX:-}; X=""
  case $ds in
    qtoy2)    X="PROJ=/mnt/beegfs/elmlund/afan/routes/toy55/simqtoy2/src.simple,DSKEYS=nparts=4 mskdiam=232 minpop=900 site_layer=yes site_states=yes kway_maxdepth=2,MAXDEPTH=1,NODEPAR=4,SEED=1,MEM=32000" ;;
    qtoyflex) X="PROJ=/mnt/beegfs/elmlund/afan/routes/toy55/simqtoyflex/src.simple,DSKEYS=nparts=4 mskdiam=232 minpop=900 site_layer=yes site_states=yes kway_maxdepth=2,MAXDEPTH=2,NODEPAR=4,SEED=1,MEM=32000" ;;
    r131)     X="FIX=$F/ancestral/fixtures/r131,DSKEYS=nparts=12 mskdiam=232 minpop=3000 site_layer=yes site_states=yes kway_maxdepth=-1 split_ntests=460,MAXDEPTH=1,NODEPAR=1,SEED=2,MEM=32000" ;;
    qtoynull) X="PROJ=/mnt/beegfs/elmlund/afan/routes/toy55/simqtoynull/src.simple,DSKEYS=nparts=4 mskdiam=232 minpop=900 site_layer=yes site_states=yes,MAXDEPTH=0,NODEPAR=1,SEED=1,MEM=32000" ;;
    olenull)  X="PROJ=/mnt/beegfs/elmlund/afan/hetero_tree_regression_20261026/null_controls/ole_null_v1/src.simple,DSKEYS=nparts=5 pgrp=c2 mskdiam=240 minpop=3000 site_layer=yes site_states=yes,MAXDEPTH=0,NODEPAR=1,SEED=1,MEM=48000" ;;
    micnull)  X="PROJ=$F/fixtures/pfcrt_micelle/micnull/src_in_fp.simple,DSKEYS=nparts=5 pgrp=c1 mskdiam=200 minpop=3000 site_layer=yes site_states=yes kway_maxdepth=2,MAXDEPTH=0,NODEPAR=4,SEED=1,MEM=32000" ;;
    ole)      X="PROJ=/mnt/beegfs/elmlund/afan/ole_12707/5_refine3D_auto/ole_12707.simple,DSKEYS=nparts=5 pgrp=c2 mskdiam=240 minpop=3000 site_layer=yes site_states=yes,MAXDEPTH=1,NODEPAR=4,SEED=1,MEM=48000" ;;
    *) echo "unknown dataset $ds"; exit 2 ;;
  esac
  test ! -e $RUN || { echo "$RUN exists, skipped"; continue; }
  # --export splits at commas, so the keys go through the environment instead
  env_assign() { while IFS= read -r kv; do export "$kv"; done < <(echo "$X" | sed 's/,\([A-Z][A-Z]*=\)/\n\1/g'); }
  unset PROJ FIX DSKEYS MAXDEPTH NODEPAR SEED MEM; env_assign
  [ -n "${SEEDX:-}" ] && SEED=$SEEDX
  jid=$(sbatch --parsable -J mrw_${ARM}_$ds --export=ALL,RUN=$RUN,BUILD=$BUILD $M/sst_gate.sbatch)
  echo "$ds -> $RUN job $jid (keys: $DSKEYS; maxdepth $MAXDEPTH nodepar $NODEPAR seed $SEED)"
done
