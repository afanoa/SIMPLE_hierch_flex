#!/bin/bash
#SBATCH -J klem_smerge
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# bottom-up organisation of a multi-state refinement's states by gold-standard distinguishability (state_merge.py).
# env RD (dir with recvol_stateXX_even/odd), LABELS (per-row final states), OUT, GT (truth) or ARM (arm_states.txt), HP (40), LP (10), ZTHR (5), UNFIL (1), SMPD (header), NCUT, BY (d2), LINK (average), BUDGET (particles)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55
$PY $T/state_merge.py ${OUT:?} --rdir ${RD:?} --labels ${LABELS:?} ${GT:+--gt $GT} --hp ${HP:-40} --lp ${LP:-10} --zthr ${ZTHR:-5} ${SMPD:+--smpd $SMPD} ${ARM:+--arm $ARM} ${NCUT:+--ncut $NCUT} ${BY:+--by $BY} ${LINK:+--link $LINK} ${BUDGET:+--budget $BUDGET} $([ "${UNFIL:-1}" = 1 ] && echo --unfil) 2>&1 | grep -v Warning
