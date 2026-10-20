#!/bin/bash
#SBATCH -J klem_server_test
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# board server test on the toy: board A (parent config) then board B (same particles, focus mask) on the SAME 4 workers,
# versus board B run alone; the cells of B must agree (cross-run corr), and the wall time shows the saving.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
FB=$HOME/SIMPLE_klem2/build; export SIMPLE_PATH=$FB PATH=$FB/bin:$FB/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8 TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; T=$A/rnasep_klem/toyah_tree; R=$A/rnasep_klem/server_test; rm -rf $R; mkdir -p $R/A $R/B $R/Bsolo; log(){ echo "[$(date '+%F %T') srv] $*"; }
# configs: A = the toy parent board (nfree 2, 20..10 A), B = same particles with the parent's region mask as focus and the region schedule 12..8 A
sed -e "s#^subset .*#subset $T/subset_all.txt#" $T/parent/board.txt | grep -v "^lps\|^sigmas\|^nthr_part\|^focus_mask\|^init_file\|^boards_more" > $R/A/board.txt
echo "lps 20 14 10" >> $R/A/board.txt; echo "sigmas 1.5 1.0 0.7" >> $R/A/board.txt; echo "init_file $T/parent/init_positions.txt" >> $R/A/board.txt; echo "nthr_part 8" >> $R/A/board.txt
cp $R/A/board.txt $R/B/board.txt; sed -i "s#^lps .*#lps 12 10 8#; s#^sigmas .*#sigmas 1.0 0.8 0.7#" $R/B/board.txt; echo "focus_mask $T/regions_parent/region_1.mrc" >> $R/B/board.txt
cp $R/B/board.txt $R/Bsolo/board.txt; echo "$R/B" > $R/A/boards_more.txt; echo "boards_more $R/A/boards_more.txt" >> $R/A/board.txt
cd $R; t0=$(date +%s); log "server: boards A then B on 4 workers"
simple_exec prg=klem_board projfile=$T/src.simple infile=$R/A/board.txt dir=$R/A nthr=8 nparts=4 mkdir=no > $R/server.log 2>&1 || { log "server FAILED: $(grep -a -m3 'ERROR\|STOP\|rror' $R/server.log | cut -c1-200)"; }
log "server done in $(( $(date +%s) - t0 )) s"; grep -a "^>>> board\|new board\|written to\|^>>> it" $R/server.log | cut -c1-140
cd $R/Bsolo; t0=$(date +%s); log "board B alone on 4 workers"
simple_exec prg=klem_board projfile=$T/src.simple infile=$R/Bsolo/board.txt dir=$R/Bsolo nthr=8 nparts=4 mkdir=no > $R/solo.log 2>&1 || log "solo FAILED"
log "solo done in $(( $(date +%s) - t0 )) s"
$PY - <<PY
import numpy as np, glob, sys
sys.path.insert(0, '$W'); from model_class import read_mrc
a = sorted(glob.glob('$R/B/cells/cell_*.mrc')); b = sorted(glob.glob('$R/Bsolo/cells/cell_*.mrc')); cc = []
for fa, fb in zip(a, b):
    x = read_mrc(fa).ravel(); y = read_mrc(fb).ravel(); x -= x.mean(); y -= y.mean(); cc.append(float(x @ y / max(np.linalg.norm(x) * np.linalg.norm(y), 1e-30)))
print('board B cells, server vs alone: %d cells, corr median %.3f min %.3f' % (len(cc), np.median(cc), np.min(cc)))
la = np.loadtxt('$R/B/klem_latent.txt', comments='#'); lb = np.loadtxt('$R/Bsolo/klem_latent.txt', comments='#'); print('argmax agreement %.3f' % np.mean(la[:, 3] == lb[:, 3]))
PY
log DONE
