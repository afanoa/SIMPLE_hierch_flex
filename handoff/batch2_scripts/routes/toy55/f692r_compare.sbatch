#!/bin/bash
#SBATCH -J f692r_compare
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=24G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# stage-by-stage comparison of the f692e reproduction (env NEW, default f692r) with the Oct 17 f692e run
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; O=$K/f692e; N=$K/${NEW:-f692r}
for f in subset_all.txt expand/family1_subset.txt expand/family2_subset.txt; do a=$(md5sum < $O/$f | cut -c1-12); b=$(md5sum < $N/$f 2>/dev/null | cut -c1-12); echo "$f: $([ "$a" = "$b" ] && echo identical || echo DIFFERENT)"; done
for d in parent fam1/child fam2/child fam1/klem_fb fam2/klem_fb $(cd $O; ls -d fam?/region*/board 2>/dev/null); do echo "== board $d"; [ -s $N/$d/klem_latent.txt ] || { echo "   missing in the new run"; continue; }; $PY $W/repro_compare.py $O/$d $N/$d 2>&1 | grep -v Warning | sed 's/^/   /'; done
for j in 1 2; do echo "== family $j refined poses"; for k in e1 e2 e3 x y; do a=$(grep -o " $k=[^ ]*" $O/fam$j/refined_ptcl3D.txt | md5sum | cut -c1-12); b=$(grep -o " $k=[^ ]*" $N/fam$j/refined_ptcl3D.txt 2>/dev/null | md5sum | cut -c1-12); echo "   $k $([ "$a" = "$b" ] && echo identical || echo DIFFERENT)"; done; done
