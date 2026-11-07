#!/bin/bash
#SBATCH -J land_real
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_real.%j.out
# real RNase P (f692st tree): per-particle coordinate inside one node, held-out window unmixing on the rigid swing and the deposited conformers
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; H=$T/simrnasep692
NODE=${NODE:-r1}; R=$K_/rflex_f692st; O=$K_/landscape/real_f692st_$NODE; DUMP=${DUMP:-$K_/lsamp/f692st_ps/dump_db61decd1108.txt}
$PY -u $T/land_real.py --tree $R --node $NODE --dump $DUMP --sims ${H}_A,${H}_B --mr ${H}_A/meta_rows.txt,${H}_B/meta_rows.txt --cache $K_/lsamp/caches/rnasep692 \
   --swing $T/maps160/swing60 --confs $T/maps160/docked_all,$T/maps160/holo_docked --out $O ${ARGS:-} 2>&1 | grep --line-buffered -v "Warning\|warn("
echo "[$(date '+%F %T')] land_real $NODE DONE"
