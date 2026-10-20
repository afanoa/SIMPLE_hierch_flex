#!/bin/bash
cd /mnt/beegfs/elmlund/afan/routes/toy55/emdb
for id in 74131 74201 74209 74230 74231 74246 74259 74132 74133 74074 ; do
  f=emd_$id.map.gz; [ -s $f ] && continue
  curl -s -L -o $f.part "https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-$id/map/emd_$id.map.gz" && mv $f.part $f && echo "$(date +%T) $f $(stat -c %s $f)" || echo "$(date +%T) FAILED $id"
done
echo "DOWNLOAD DONE"
