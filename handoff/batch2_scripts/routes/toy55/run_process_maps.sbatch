#!/bin/bash
#SBATCH -J toy55_maps
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# wait for the 55 downloads, resample/mirror to our lattice, dock every conformer into the consensus frame (reference: the
# already-docked conformer 1), write the image-count-weighted mean as the toy consensus.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$PATH OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; cd $W; log(){ echo "[$(date '+%F %T')] $*"; }
n=0; while [ $n -lt 180 ] && ! grep -q "DOWNLOAD DONE" $W/emdb/dl.log 2>/dev/null; do sleep 60; n=$((n+1)); done
log "downloads: $(ls $W/emdb/*.map.gz 2>/dev/null | wc -l) files, failed: $(grep -c FAILED $W/emdb/dl.log)"
$PY $W/process_maps.py $W/${CSV:-rnasep_conformers10.csv} $W/emdb $W/maps160 || { log "process FAILED"; exit 1; }
REF=$A/rnasep_coupled/degenhardt/emd_74131_m160_docked.mrc; mkdir -p $W/maps160/docked
for m in $W/maps160/conf*_m160.mrc; do b=$(basename $m .mrc)
  [ -s $W/maps160/docked/$b.mrc ] && continue
  simple_exec prg=dock_volpair vol1=$REF vol2=$m smpd=1.464 lp=10 hp=120 mskdiam=200 trs=12 outvol=$W/maps160/docked/$b.mrc outfile=$W/maps160/docked/$b.txt nthr=8 mkdir=no > $W/maps160/docked/$b.log 2>&1 \
    || { log "dock FAILED $b"; continue; }
  echo "$b $(grep '^cc=' $W/maps160/docked/$b.txt) $(grep '^euler=' $W/maps160/docked/$b.txt) $(grep '^shift=' $W/maps160/docked/$b.txt)"
done
$PY - <<PYEOF
import numpy as np, csv, glob, os
def rd(fn):
    h = np.fromfile(fn, np.int32, 256); nx, ny, nz = h[:3]; return np.fromfile(fn, np.float32, offset=1024 + int(h[23]), count=nx*ny*nz).reshape(nz, ny, nx).astype(np.float64)
def wr(fn, v, smpd=1.464):
    v = np.ascontiguousarray(v, np.float32); nz, ny, nx = v.shape; h = np.zeros(256, np.int32); hf = h.view(np.float32)
    h[0:3] = (nx, ny, nz); h[3] = 2; h[7:10] = (nx, ny, nz); hf[10:13] = (nx*smpd, ny*smpd, nz*smpd); hf[13:16] = 90.; h[16:19] = (1, 2, 3)
    hf[19:22] = (v.min(), v.max(), v.mean()); h[52] = int.from_bytes(b'MAP ', 'little'); h[53] = int.from_bytes(bytes([68, 65, 0, 0]), 'little'); hf[54] = v.std()
    open(fn, 'wb').write(h.tobytes() + v.tobytes())
rows = list(csv.DictReader(open("$W/${CSV:-rnasep_conformers10.csv}"))); tot = 0.; acc = None; cc = []
for r in rows:
    fn = glob.glob('$W/maps160/docked/conf%02d_*.mrc' % int(r['conformer']))
    if not fn: continue
    v = rd(fn[0]); w = float(r['images']); acc = v * w if acc is None else acc + v * w; tot += w
    cc.append(float(open(fn[0][:-4] + '.txt').read().split('cc=')[1].split()[0]))
wr('$W/maps160/toy_consensus.mrc', acc / tot); print('consensus from %d docked conformers (image-weighted); docking cc min %.3f median %.3f' % (len(cc), min(cc), np.median(cc)))
ref = rd('$REF'); m = np.hypot(np.hypot(*np.meshgrid(np.arange(160)-80, np.arange(160)-80, indexing='ij')), 0) < 68
X = (acc / tot); print('toy consensus vs docked conformer 1: corr %.3f' % np.corrcoef(X[m].ravel() if X.ndim == 2 else X.ravel(), ref.ravel())[0, 1])
PYEOF
log DONE
