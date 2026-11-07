#!/bin/bash
#SBATCH -J klem_lp10
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 00:15:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# 10 A low-passed, unsharpened copies (SIMPLE filter, lp=10, no B-factor) of motion-series maps; env SITES (dirs holding
# bin1..5.mrc); writes <dir>/lp10/bin1..5.mrc at the input box and pixel size and prints both headers
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_stree/snap/${SNAP:-m8c52t}; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_EMAIL=nil OMP_NUM_THREADS=4
PY=/mnt/nasapps/production/python/3.12/bin/python3
hdr(){ $PY -c "import sys, numpy as np; h = np.frombuffer(open(sys.argv[1], 'rb').read(1024), np.int32); f = h.view(np.float32); print(int(h[0]), int(h[1]), int(h[2]), '%.6f' % (f[10] / h[0]))" $1; }
for d in ${SITES:?}; do mkdir -p $d/lp10; cd $d/lp10
  for k in 1 2 3 4 5; do read nx ny nz smpd <<< "$(hdr $d/bin$k.mrc)"
    simple_exec prg=filter vol1=$d/bin$k.mrc outvol=$d/lp10/bin$k.mrc smpd=$smpd lp=10 nthr=4 mkdir=no > filter_bin$k.log 2>&1 || echo "FAILED $d bin$k: $(grep -a -m2 ERROR filter_bin$k.log)"
    echo "$d bin$k: in $nx $ny $nz $smpd -> out $(hdr $d/lp10/bin$k.mrc 2>/dev/null || echo missing)"; done; done
