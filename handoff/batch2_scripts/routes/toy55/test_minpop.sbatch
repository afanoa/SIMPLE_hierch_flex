#!/bin/bash
#SBATCH -J klem_minpop_test
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# node_states2.py with and without the merge floor on finished boards: real family 1 and 2 region boards (932,496 rows,
# floor 20k) and the toy family child boards (floor 8k). Prints the state tables side by side; writes nothing into the runs.
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; O=$A/rnasep_klem/minpop_test; mkdir -p $O
for spec in "f692/fam1/region_1/board 932496 20000" "f692/fam2/region_1/board 932496 20000" "f692/fam1/klem_fb 932496 20000" "toyah_tree/fam1/child 0 8000" "toyah_tree/fam2/child 0 8000"; do
  set -- $spec; b=$A/rnasep_klem/$1; tag=$(echo $1 | tr / _)
  for mp in off on; do echo "=================== $1  support $mp"; $PY $W/node_states2.py $b $O/${tag}_$mp --nrows $2 --sim 0.7 --nmax 30000 $([ $mp = on ] && echo --support) 2>&1 | grep -v Warning; done
done
echo "ALL DONE"
