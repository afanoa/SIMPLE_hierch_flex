#!/bin/bash
#SBATCH -J klem_site_comp2
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 00:40:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# reproducible components of an existing site_<TAG> group set, support-restricted and amplitude-matched (site_components2.py)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; O=$K_/site_${TAG:?}
$PY $T/site_components2.py $O/comps2 --list $O/list.txt --mskdiam ${MSKDIAM:-220} ${TRUTH:+--truth $TRUTH} 2>&1 | grep -v "Warning\|lpf = "
