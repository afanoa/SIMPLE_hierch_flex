#!/bin/bash
#SBATCH -J klem_dk_merge
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# view-matched merge test (env RTAG, MODE): (a) every leaf re-reconstructed from a VIEW-MATCHED subset of its particles
# (identity labels, VIEWMATCH=1) -> consolidated_12_vm; (b) the data-derived merge (fixed 12 A and band-scan auto) run
# on THOSE half maps (envelope mask); (c) truth purity / junk of the merged states
set -uo pipefail; PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; RTAG=${RTAG:?}; MODE=${MODE:?}; R=$A/rnasep_klem/$RTAG; GT=$W/sim$MODE/gt_particles.txt
export DOCK=1 SPECMATCH=${SPECMATCH:-0} RTAG LEAFGLOB="region*" FLOOR=12 MASK=region WHITEN=0 UNCAP=1 VIEWMATCH=1 META=$R/meta_toy.npz TAG=vm CUT=12 CTAG=_vm MINP=1000 TREE=$R/merge_rule_12_12.npz
[ -s $TREE ] || { echo "no identity tree $TREE"; exit 1; }
echo "== (a) view-matched leaf half maps already present"
C=$R/consolidated_12_vm; ls $C/fam0/regionall/level1/maps/recvol_state??_even.mrc | wc -l | sed 's/^/view-matched leaf half maps: /'
for FL in 12 auto; do echo "== (b) merge of the view-matched leaves, floor $FL"; FLOOR=$FL $PY $W/merge_data_derived.py $C $W/maps160/toy_consensus.mrc $C/merge_dk_$FL 2>&1 | grep -v Warning | grep -a "leaf states\|merges ->\|remaining\|MISSING\|dock " | head -4 | cut -c1-200
  $PY - <<PYEOF
import numpy as np
M = np.load('$C/merge_dk_$FL.npz', allow_pickle=True); names, lab = list(M['names']), M['labels']; st = np.loadtxt('$C/states.txt', dtype=int); out = np.zeros_like(st)
for i, nm in enumerate(names):
    s = int(nm.split('_s')[-1]); out[st == s] = lab[i]
np.savetxt('$C/states_merged_dk_$FL.txt', out, fmt='%d'); print('%d leaf states -> %d merged (view-matched, docked)' % (len(names), lab.max()))
PYEOF
  mkdir -p $C/mdk_$FL; cp -f $C/states_merged_dk_$FL.txt $C/mdk_$FL/states.txt; echo "== (c) truth of the merged states, floor $FL"; $PY $W/consolidated_purity.py $C/mdk_$FL $GT 2>&1 | grep -v Warning; done
echo "DK MERGE DONE"
