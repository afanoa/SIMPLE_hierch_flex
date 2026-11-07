#!/bin/bash
#SBATCH -J klem_lsamp_wait
#SBATCH -p csbdevel
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 24:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# waits (read-only) for the f692st tree of the state-tree session to finish (leaves_list.txt + final_labels.txt with every master row, and
# its master job WAITJOB gone), then submits the landscape-sampling chain on it (lsamp_chain.sbatch -> rnasep_klem/lsamp/f692st_ps,
# maps in view_rflex/lsamp_f692st). Touches nothing of f692st. env WAITJOB (54050336)
set -uo pipefail
A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; R=$K_/rflex_f692st; J=${WAITJOB:-54050336}
for i in $(seq 1 140); do
  if [ -s $R/leaves_list.txt ] && [ -s $R/final_labels.txt ] && [ "$(wc -l < $R/final_labels.txt)" -eq 932496 ] && ! squeue -h -j $J 2>/dev/null | grep -q .; then
    sleep 300; echo "[$(date '+%F %T')] f692st finished ($(wc -l < $R/leaves_list.txt) leaves): submitting the landscape chain"
    cd $K_/lsamp && TREE=$R OUT=$K_/lsamp/f692st_ps SIM=$T/simrnasep692_A,$T/simrnasep692_B MR=$T/simrnasep692_A/meta_rows.txt,$T/simrnasep692_B/meta_rows.txt \
      CP=$K_/lsamp/caches/rnasep692 CTFABS=1 MSKDIAM=232 VIEWBASE=$K_/view_rflex/lsamp_f692st sbatch --export=ALL $T/lsamp_chain.sbatch
    exit 0
  fi
  sleep 600
done
echo "[$(date '+%F %T')] f692st not finished within the wait"
