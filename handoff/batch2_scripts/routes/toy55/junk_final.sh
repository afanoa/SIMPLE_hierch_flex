#!/bin/bash
# final-rule (xd3) reruns of junk_score.py on saved half labels / poses -> rnasep_klem/factorial/junk_test/final/<name>; arg: toys | real
K=/mnt/beegfs/elmlund/afan/rnasep_klem; T=/mnt/beegfs/elmlund/afan/routes/toy55; JT=$K/factorial/junk_test; PY=/mnt/nasapps/production/python/3.12/bin/python3; F=$JT/final; mkdir -p $F
run(){ nm=$1; shift; mkdir -p $F/$nm; $PY -u $T/junk_score.py --out $F/$nm --tag $nm "$@" 2>&1 | grep -v Warning > $F/$nm/junk.log; echo "#### $nm"; grep -v "^  guard a[1-5]_ph" $F/$nm/junk.log; }
SS="--sim $T/simsites_A,$T/simsites_B --mr $T/simsites_A/meta_rows.txt,$T/simsites_B/meta_rows.txt --nmaster 199996 --cache $T/benchcache_sites --ctfabs 0 --truth $JT/simsites_truth_jt.txt --jtype jtype --guard C_closed,A_domain"
RN="--sim $T/simrnasep692_A,$T/simrnasep692_B --mr $T/simrnasep692_A/meta_rows.txt,$T/simrnasep692_B/meta_rows.txt --nmaster 932496 --cache $T/benchcache_rnasep692 --ctfabs 1"
if [ "$1" = toys ]; then
  run simsites_K1 $SS --poses $JT/simsites/poses_A.txt,$JT/simsites/poses_B.txt --labels $JT/simsites/K1/lab_A.txt,$JT/simsites/K1/lab_B.txt
  run simsites_convC $SS --poses $K/factorial/simsites_sites_d1/poses_A.txt,$K/factorial/simsites_sites_d1/poses_B.txt --labels $JT/simsites_convC/convC/lab_A.txt,$JT/simsites_convC/convC/lab_B.txt
  run igm_K1 --sim $T/simarmsigm1_A,$T/simarmsigm1_B --mr $K/armsigm1_src/meta_rows_A.txt,$K/armsigm1_src/meta_rows_B.txt --nmaster 200000 --cache $T/benchcache_armsigm1 --ctfabs 0 \
    --poses $JT/igm/poses_A.txt,$JT/igm/poses_B.txt --labels $JT/igm/K1/lab_A.txt,$JT/igm/K1/lab_B.txt --truth $K/factorial/armsigm1_truth/arms_truth.txt --guardtail a1_th,a2_th,a3_th,a4_th,a5_th
else
  run r1conv_c1 $RN --poses $JT/r1conv/poses_A.txt,$JT/r1conv/poses_B.txt --labels $JT/r1conv/c1/lab_A.txt,$JT/r1conv/c1/lab_B.txt
  run f692st_r_states $RN --poses $JT/f692st_r/poses_A.txt,$JT/f692st_r/poses_B.txt --labels $JT/f692st_r/r/lab_A.txt,$JT/f692st_r/r/lab_B.txt
fi
