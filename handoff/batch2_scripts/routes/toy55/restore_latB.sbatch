#!/bin/bash
#SBATCH -J restore_latB
#SBATCH -p csbdevel
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# lat_55_700k_B/ptcl3D_r0.txt was overwritten (02:30) through a fallback symlink by oracle_refine 54013658: keep that result
# as the oracle's own file, restore the original from the 01:00 copy (final_xem_cons/poses_B.txt), verify against src.simple B
W=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; L=$W/lat_55_700k_B/ptcl3D_r0.txt; C=$K/dtree_55_700k_sixc8r/final_xem_cons/poses_B.txt
OO=$K/oracle_55_700k_sixc8/B/refined; OP=$K/polish_55_700k_sixc8r/B/refined
ls -la --time-style=full-iso $L $C
rm -f $OO/poses_full.txt && cp $L $OO/poses_full.txt && echo "oracle B result kept: $(grep -c . $OO/poses_full.txt) rows"
rm -f $OP/poses_full.txt; echo "polish B result was overwritten: symlink removed (rerun needed)"
cp $C $L.tmp && mv $L.tmp $L && echo "restored $L from $C"
for k in e1 e2 e3 x y dfx; do echo "$k $(grep -o " $k=[^ ]*" $L | md5sum | cut -c1-12)"; done
echo "expected (src.simple B, projcheck 54013583): e1 af686e543c86 e2 8c7a09163c99 e3 729cba2cfde4 x d7320c631312 y 4c5e46127f98 dfx 30cacf8047ff"
