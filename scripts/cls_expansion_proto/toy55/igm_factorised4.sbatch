#!/bin/bash
#SBATCH -J igm_factorised4
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 08:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Truth-free factorised route, merged stage v4 (truth only in readouts): the 12 candidate regions + boards of
# igm_factorised2/3 merged with the gap rule at a GUARD-sd guard (region_merge.py gap:GUARD) -> NSTART random starts per
# merged region as one SLURM array (igm_mstart_board.sbatch, CONC at a time) -> per region the most confident start
# (start_select.py, mean max responsibility) -> factorised_readout.py. env MODE (armsigm1), GUARD (1.5), NSTART (4), CONC (16)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; KL=$A/rnasep_klem; MODE=${MODE:-armsigm1}
O1=$W/factorised_$MODE; O2=$W/factorised2_$MODE; O3=$W/factorised3_$MODE; O=$O3/merged4; mkdir -p $O; cd $O; NS=${NSTART:-4}
log(){ echo "[$(date '+%F %T') factorised4 $MODE] $*"; }
CANDS=""; for f in $(ls $O2/regions/region_*.mrc); do r=$(basename $f .mrc); r=${r#region_}; [ -s $O2/arm_$r/board/klem_latent.txt ] && CANDS="$CANDS $f:$O2/arm_$r/board"; done
for f in $(ls $O3/resid/regions/region_*.mrc); do r=$(basename $f .mrc); r=${r#region_}; [ -s $O3/resid/arm_$r/board/klem_latent.txt ] && CANDS="$CANDS $f:$O3/resid/arm_$r/board"; done
[ -s $O/regions/region_1.mrc ] || $PY $W/region_merge.py $O/regions gap:${GUARD:-1.5} $CANDS 2>&1 | grep -v Warning | grep -v "^  "
NR=$(ls $O/regions/region_*.mrc | wc -l); log "$NR merged regions -> $NS starts each ($((NR * NS)) boards, ${CONC:-16} at a time)"
J=$(sbatch --parsable --array=1-$((NR * NS))%${CONC:-16} --export=ALL,RUN=$O,SUB=$O1/clean_subset.txt,PS=$KL/synth_$MODE/src.simple,NSTART=$NS $W/igm_mstart_board.sbatch)
while squeue -h -j $J 2>/dev/null | grep -q .; do sleep 30; done; log "starts done ($(cat $W/igm_mstart.${J}_*.out 2>/dev/null | grep -c done) of $((NR * NS)))"
echo "== start selection (mean max responsibility; * = kept; truth only in brackets)"; $PY $W/start_select.py $O --gt $W/sim$MODE/gt_particles.txt --truth $W/truth_$MODE.txt 2>&1 | grep -v Warning
echo "== merged stage readout (truth)"; $PY $W/factorised_readout.py $O $W/sim$MODE/gt_particles.txt $W/truth_$MODE.txt $W/maps160/synth_$MODE/masks 2>&1 | grep -v Warning
log DONE
