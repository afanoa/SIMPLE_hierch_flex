#!/bin/bash
cd /mnt/beegfs/elmlund/afan/routes/toy55/emdb
for id in 74131 74140 74141 74195 74196 74197 74198 74199 74200 74201 74202 74203 74205 74204 74206 74207 74208 74209 74213 74215 74223 74224 74226 74227 74228 74229 74230 74231 74232 74235 74237 74238 74239 74240 74243 74244 74246 74247 74248 74249 74250 74251 74252 74254 74257 74259 74261 74132 74133 74255 74074 74263 74260 74264 74137 ; do
  f=emd_$id.map.gz; [ -s $f ] && continue
  curl -s -L -o $f.part "https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-$id/map/emd_$id.map.gz" && mv $f.part $f && echo "$(date +%T) $f $(stat -c %s $f)" || echo "$(date +%T) FAILED $id"
done
echo "DOWNLOAD DONE"
