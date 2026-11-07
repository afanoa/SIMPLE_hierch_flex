#!/bin/bash
#SBATCH -J klem_site_reposed
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=160G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# BASE's site J EM at the poses of PSRC (rows checked against SRC0) from the start SEEDF, compared with REF (site_reseed.sbatch); env RT NODE PSRC SRC0 BASE J REF SEEDF FLEXRANK OUT MSKDIAM
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local
A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; H=${H:-$T/simrnasep692}; MRA=${MRA:-${H}_A/meta_rows.txt}; MRB=${MRB:-${H}_B/meta_rows.txt}
BASE=${BASE:?}; J=${J:-2}; NB=${OUT:?}; PSRC=${PSRC:?}; SRC0=${SRC0:?}; mkdir -p $NB/picked $NB/site_$J $NB/reseed; cd $NB; log(){ echo "[$(date '+%F %T') site reposed] $*"; }
if [ ! -s $NB/poses_B.txt ]; then
  simple_exec prg=print_project_field projfile=$PSRC oritype=ptcl3D mkdir=no 2>/dev/null > $NB/pm.txt; simple_exec prg=print_project_field projfile=$SRC0 oritype=ptcl3D mkdir=no 2>/dev/null > $NB/p0.txt
  key(){ awk '{k = ""; for (i = 1; i <= NF; i++) if ($i ~ /^(stkind|indstk|pind)=/) k = k " " $i; print k}' $1; }
  [ $(grep -c . $NB/pm.txt) = $(grep -c . $NB/p0.txt) ] && cmp -s <(key $NB/pm.txt) <(key $NB/p0.txt) || { log "PSRC rows do not match SRC0 row for row"; exit 1; }
  awk 'NR == FNR {u[FNR] = $1; next} u[FNR] > 0 {n++; for (i = 1; i <= NF; i++) if ($i ~ /^state=/) {split($i, a, "="); if (a[2] + 0 == 0) z++}} END {printf "rows in BASE rows_use %d, state 0 at PSRC %d\n", n, z}' $BASE/rows_use.txt $NB/pm.txt
  awk 'NR == FNR {l[FNR] = $0; next} {print l[$1 + 0]}' $NB/pm.txt $MRA > $NB/poses_A.txt; awk 'NR == FNR {l[FNR] = $0; next} {print l[$1 + 0]}' $NB/pm.txt $MRB > $NB/poses_B.txt; rm -f $NB/pm.txt $NB/p0.txt; fi
ln -sfn $BASE/rows_use.txt $NB/rows_use.txt; ln -sfn $(readlink -f $BASE/picked/site_$(printf %02d $J).mrc) $NB/picked/site_$(printf %02d $J).mrc; ln -sfn ${REF:?} $NB/site_$J/post_master.txt
[ -n "${FLEXRANK:-}" ] && cp -f $FLEXRANK $NB/reseed/flexrank.txt
log "site $J of $BASE at the poses of $PSRC, start $SEEDF, reference ranking $REF"
RT=${RT:?} NODE=${NODE:-r} SRC=$PSRC BASE=$NB J=$J MSKDIAM=${MSKDIAM:-220} OUT=$NB/reseed SEEDS="reposed=${SEEDF:?}@1" PARS=1 bash $T/site_reseed.sbatch
