#!/bin/bash
#SBATCH -J merge_test
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# dry run of region_merge.py on the v2 region boards (scratch output)
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; O2=$W/factorised2_${MODE:-armsigm1}
C=""; for f in $(ls $O2/regions/region_*.mrc); do r=$(basename $f .mrc); r=${r#region_}; C="$C $f:$O2/arm_$r/board"; done
$PY $W/region_merge.py $O2/merge_test ${THR:-auto} $C; $PY $W/region_merge.py $O2/merge_test_nuis ${THR:-auto} --nuis $W/lat_${MODE:-armsigm1}/ptcl3D_r0.txt $C 2>&1 | grep -v Warning; echo MERGE_TEST DONE
