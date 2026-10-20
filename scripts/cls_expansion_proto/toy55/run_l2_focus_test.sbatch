#!/bin/bash
#SBATCH -J klem_l2_focus_test
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# sequential (host) focused board timing after the parallel focus: the toy level-2 board of region_1 (14k particles, 16 threads)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; T=$A/rnasep_klem/toyah_srv/fam1; L=$T/region_1/level2_s1; R=$A/rnasep_klem/l2_focus_test; rm -rf $R; mkdir -p $R
PS=$(ls -td $T/refine/*_refine3D_states | head -1)/p.simple; t0=$(date +%s)
FB=$HOME/SIMPLE_klem2/build NTHR=16 SERVER_LIST= NPARTS= FOCUSMASK=$T/regions/region_1.mrc META=$T/meta_fb.npz bash $W/fboard.sh $R/board $PS $L/subset.txt "12 12 10 10 8 8 8 8 8" "1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7"
echo "total $(( $(date +%s) - t0 )) s (was 1193 s)"; grep -a "^>>> it" $R/board/board.log | cut -c1-60,150-200
