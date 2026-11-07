#!/bin/bash
#SBATCH -J refstates_nmi2
#SBATCH -p csbdevel
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# truth readout of a refinement's re-classified states (refine3D_states state field) vs its input labels (env O = dir with A/ B/ refined/poses_full.txt + states_{A,B}.txt, MODE)
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; MODE=${MODE:-55_700k}; O=${O:?}; GT=$W/sim${MODE}_A/gt_particles.txt,$W/sim${MODE}_B/gt_particles.txt; TR="$W/truthtree_55_200k.npz 0.45"
for h in A B; do awk 'NF {v = 0; for (i = 1; i <= NF; i++) if ($i ~ /^state=/) {sub(/^state=/, "", $i); v = int($i + 0.5)}; print v}' $O/$h/refined/poses_full.txt > $O/refstates_$h.txt; done
$PY $W/final_nmi.py input_labels $O/states_A.txt,$O/states_B.txt $GT $TR 2>&1 | grep -v Warning
$PY $W/final_nmi.py refined_states $O/refstates_A.txt,$O/refstates_B.txt $GT $TR 2>&1 | grep -v Warning
for h in A B; do echo "== half $h refined-state composition"; $PY $W/dtree_compose.py $O/refstates_$h.txt $W/sim${MODE}_$h/gt_particles.txt 2>&1 | grep -v Warning; done
$PY - <<PYEOF
import numpy as np
for h in 'AB':
    a = np.loadtxt('$O/states_%s.txt' % h).astype(int); b = np.loadtxt('$O/refstates_%s.txt' % h).astype(int); m = a > 0
    print('half %s: %.1f %% of the labelled particles changed state in the refinement; sizes in %s -> out %s' % (h, 100 * np.mean(a[m] != b[m]), np.bincount(a[m])[1:].tolist(), np.bincount(b[m])[1:].tolist()))
PYEOF
