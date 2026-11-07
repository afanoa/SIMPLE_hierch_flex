#!/bin/bash
#SBATCH -J klem_eostates
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# even/odd composition of every region level-1 state and zoom sub-state of run RTAG (eo per row cached in famF/eo_rows.txt)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH
R=/mnt/beegfs/elmlund/afan/rnasep_klem/${RTAG:?}; T=$(mktemp -d /tmp/eost.XXXXXX); cd $T
for F in 1 2; do D=$R/fam$F; PS=$(ls -td $D/refine/*_refine3D_states | head -1)/p.simple; EO=$D/eo_rows.txt
  [ -s $EO ] || { simple_exec prg=print_project_field projfile=$PS oritype=ptcl3D mkdir=no 2>/dev/null | awk '{e = -1; for (i = 1; i <= NF; i++) if ($i ~ /^eo=/) {split($i, a, "="); e = a[2] + 0}; print e}' > $EO.$$ && mv -f $EO.$$ $EO; }
  echo "== $RTAG fam$F: $(grep -c . $EO) rows, even $(grep -c '^0$' $EO), odd $(grep -c '^1$' $EO)"
  for L in $D/region*/level1 $D/zoomregion*/level1; do sc=$L/states_clusters.txt; [ -s $sc ] || continue
    awk -v t=$(basename $(dirname $L)) 'FILENAME == ARGV[1] {e[FNR] = $1; next} $1 > 0 {n[$1]++; if (e[FNR] == 0) ne[$1]++} END {printf "%s:", t; for (s in n) {f = ne[s] / n[s]; printf " s%d %d %.2f", s, n[s], f; if (n[s] >= 1000 && (f < 0.35 || f > 0.65)) bad++; if (n[s] >= 1000) tot++}; printf "  | skewed(>=1000, outside 35-65 %%): %d/%d\n", bad, tot}' $EO $sc
  done
done
cd /; rm -rf $T
