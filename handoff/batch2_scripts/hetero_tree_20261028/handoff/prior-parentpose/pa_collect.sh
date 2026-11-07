#!/bin/bash
# one block per node: the label readouts (production, each replay), agreement with production, drift (in -> B), r131 map readout
A=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/anchor
for key in $(ls $A/replays | sed -E 's/_(def|anc)_s[0-9]+$//' | sort -u); do
  echo "########## $key"
  first=1
  for r in $(ls -d $A/replays/${key}_*); do
    n=$(basename $r); arm=${n#${key}_}
    [ -s $r/readout/labels.txt ] || { echo "  $arm: (no readout yet)"; continue; }
    if [ $first = 1 ]; then grep '^production' $r/readout/labels.txt | sed 's/^/  /'; first=0; fi
    grep -E "^${arm} " $r/readout/labels.txt | sed 's/^/  /'
    grep -E "^production vs" $r/readout/labels.txt | sed 's/^/    /'
    grep -A3 '^in->B' $r/readout/drift.txt | sed -n '1,4p' | sed 's/^/    drift /'
    [ -s $r/readout/arm.txt ] && grep -v '^#' $r/readout/arm.txt | awk '{n=split($1,p,"/"); print "    map " p[n], $2, $3, $4}'
  done
done
