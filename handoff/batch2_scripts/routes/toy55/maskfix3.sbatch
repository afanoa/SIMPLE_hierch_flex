#!/bin/bash
#SBATCH -J klem_maskfix3
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# soft-edge (falloff) node masks on unit 7's cheap children: arm-tip weight per radius and total weight
K=/mnt/beegfs/elmlund/afan/rnasep_klem; T=/mnt/beegfs/elmlund/afan/routes/toy55; M=$T/maps160; P=/mnt/nasapps/production/python/3.12/bin/python3
L=$K/rflex_f692u7/node_r/cheap/mask_list.txt; O=$K/maskfix; cd $O
for spec in 0.2:3:6 0.2:3:9 0.1:3:6 0.1:2:8; do IFS=: read fr dl sf <<< "$spec"
  $P $T/arm_mask2.py $O/u7_f${fr}_d${dl}_s${sf}.mrc --list $L --mskdiam 232 --frac $fr --dil $dl --soft $sf 2>&1 | grep -v Warn | tail -1; done
$P - <<PYEOF
import numpy as np
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
for nm in ['u7_f0.1_d6', 'u7_f0.2_d3_s6', 'u7_f0.2_d3_s9', 'u7_f0.1_d3_s6', 'u7_f0.1_d2_s8']:
    m = rd('$O/%s.mrc' % nm); print('%-16s total weight %8.0f voxels, voxels > 0.5 %7d, > 0.05 %7d' % (nm, m.sum(), (m > 0.5).sum(), (m > 0.05).sum()))
PYEOF
$P $T/radial_extent.py $M/toy_consensus.mrc $M/docked/conf18_emd74209_m160.mrc --masks f0.1_d6_s3:$O/u7_f0.1_d6.mrc f0.2_d3_s6:$O/u7_f0.2_d3_s6.mrc f0.2_d3_s9:$O/u7_f0.2_d3_s9.mrc f0.1_d3_s6:$O/u7_f0.1_d3_s6.mrc f0.1_d2_s8:$O/u7_f0.1_d2_s8.mrc --mskdiam 232 2>&1 | grep -v Warn | sed -n 2,11p
