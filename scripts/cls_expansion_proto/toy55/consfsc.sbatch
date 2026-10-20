#!/bin/bash
#SBATCH -J klem_consfsc
#SBATCH -p csbdevel
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 00:05:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
/mnt/nasapps/production/python/3.12/bin/python3 - <<PYEOF
import numpy as np, glob
for f in glob.glob("/mnt/beegfs/elmlund/afan/rnasep_states/25_refine3D_states/fsc_state0*.bin"):
    v = np.fromfile(f, np.float32); n = len(v); k = [i for i, x in enumerate(v) if x < 0.143 and i > 2]; kk = k[0] if k else n
    print("%s: %d shells, FSC 0.143 at shell %d -> %.2f A (box 160 @ 1.464); FSC at 8/6/5/4 A: %s" % (f.split("/")[-1], n, kk, 160 * 1.464 / kk, " ".join("%.2f" % v[int(round(160 * 1.464 / r))] for r in (8, 6, 5, 4))))
print("CONSFSC DONE")
PYEOF
