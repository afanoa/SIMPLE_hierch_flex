#!/bin/bash
#SBATCH -J klem_loadtest
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# loader profile: the real 100k subset through the distributed board (10 workers) at the fine band (12->8 A, kmax 29)
# with the instrumented loader; the worker logs print stack-read vs FFT/extraction seconds
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/loadtest; mkdir -p $R/klem2_bin; cp -r $HOME/SIMPLE_klem2/build/bin $R/klem2_bin/ 2>/dev/null
export FB=$R/klem2_bin NTHR=8 NPARTS=10 CTFABS=1 PSCALE=1; t0=$(date +%s)
NFREE=2 META=$A/rnasep_coupled/basis_span_all/meta_all.npz bash $W/fboard.sh $R/board $A/rnasep_coupled/full692k_2d/rnasep.simple $A/rnasep_klem/nb_real/subset100k.txt "12 10 8" "1.0 0.8 0.7" || echo "board FAILED"
echo "total $(( $(date +%s) - t0 )) s"; grep -a -h "klem_load_particles" $R/board/klem_worker_*.log | head -10; grep -a "^>>> it\|^>>> band" $R/board/board.log | cut -c1-160
