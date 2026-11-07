#!/bin/bash
#SBATCH -J rnasep_apo_rest
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# The remaining deposited apo conformers: list every EMDB entry EMD-74070..74300 through the API (title "conformer N"),
# keep the apo ones (no protein in the title/sample) not yet docked, download + verify, resample/mirror (process_maps.py),
# dock into the apo consensus frame (dock_volpair vs docked conformer 1) into maps160/docked_all (with the 10 we had),
# then the identification readout of the real run's family states (refined) and region level-1 states against ALL.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$PATH OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; cd $W; log(){ echo "[$(date '+%F %T') apo_rest] $*"; }
LIST=$W/rnasep_emdb_all.csv
if [ ! -s $LIST ]; then log "listing EMDB entries 74070..74300"
$PY - <<PY
import json, re, subprocess, csv, time
rows = []
for i in range(74070, 74301):
    try: j = json.loads(subprocess.run(['curl', '-s', '-m', '30', 'https://www.ebi.ac.uk/emdb/api/entry/EMD-%d' % i], capture_output=True, text=True).stdout)
    except Exception: continue
    if not isinstance(j, dict) or 'admin' not in j: continue
    t = j['admin'].get('title', ''); m = re.search(r'conformer\s+(\d+)', t, re.I)
    if 'RNase P' not in t and 'Ribonuclease P' not in t: continue
    sample = json.dumps(j.get('sample', {})).lower(); holo = ('rnpa' in sample or 'protein' in sample or 'rnpa' in t.lower() or 'protein' in t.lower())
    try: rec = j['structure_determination_list']['structure_determination'][0]['image_processing'][0]['final_reconstruction']; res = rec['resolution']['valueOf_']; nimg = rec.get('number_images_used', 0)
    except Exception: res, nimg = '', 0
    try: px = j['map']['pixel_spacing']['x']['valueOf_']; box = j['map']['dimensions']['col']
    except Exception: px, box = 0.732, 480
    rows.append(dict(emdb_id='EMD-%d' % i, conformer=int(m.group(1)) if m else -1, images=nimg, resolution=res, pixel=px, box=box, file='emd_%d.map.gz' % i, holo=int(holo), title=t)); time.sleep(0.2)
w = csv.DictWriter(open('$LIST', 'w'), fieldnames=['emdb_id', 'conformer', 'images', 'resolution', 'pixel', 'box', 'file', 'holo', 'title']); w.writeheader(); [w.writerow(r) for r in rows]
print('%d RNase P entries, %d apo conformers, %d holo' % (len(rows), sum(1 for r in rows if not r['holo']), sum(1 for r in rows if r['holo'])))
PY
fi
# apo conformers not yet docked (the 10 in rnasep_conformers10.csv)
HAVE=$(tail -n +2 $W/rnasep_conformers10.csv | cut -d, -f2 | tr '\n' ' '); CSV=$W/rnasep_apo_rest.csv
$PY - <<PY
import csv
have = set(int(v) for v in '$HAVE'.split()); rows = [r for r in csv.DictReader(open('$LIST')) if r['holo'] == '0' and int(r['conformer']) > 0 and int(r['conformer']) not in have]
w = csv.DictWriter(open('$CSV', 'w'), fieldnames=['emdb_id', 'conformer', 'images', 'resolution', 'pixel', 'box', 'file']); w.writeheader(); [w.writerow({k: r[k] for k in w.fieldnames}) for r in rows]; print('%d apo conformers to add: %s' % (len(rows), ' '.join(r['conformer'] for r in rows)))
PY
mkdir -p $W/emdb $W/maps160/apo_rest $W/maps160/docked_all
for id in $(tail -n +2 $CSV | cut -d, -f1 | sed 's/EMD-//'); do f=$W/emdb/emd_$id.map.gz
  if [ -s $f ] && gzip -t $f 2>/dev/null; then continue; else rm -f $f; fi
  for try in 1 2 3 4; do wget -q -c -O $f https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-$id/map/emd_$id.map.gz && gzip -t $f 2>/dev/null && break; sleep 30; done
  gzip -t $f 2>/dev/null || { log "download of emd_$id.map.gz FAILED verification"; rm -f $f; }
done; log "downloads done"
$PY $W/process_maps.py $CSV $W/emdb $W/maps160/apo_rest || { log "process FAILED"; exit 1; }
REF=$A/rnasep_coupled/degenhardt/emd_74131_m160_docked.mrc; cp -n $W/maps160/docked/conf??_*_m160.mrc $W/maps160/docked_all/ 2>/dev/null
for m in $W/maps160/apo_rest/conf*_m160.mrc; do b=$(basename $m .mrc); [ -s $W/maps160/docked_all/$b.mrc ] && continue
  simple_exec prg=dock_volpair vol1=$REF vol2=$m smpd=1.464 lp=10 hp=120 mskdiam=200 trs=12 outvol=$W/maps160/docked_all/$b.mrc outfile=$W/maps160/docked_all/$b.txt nthr=8 mkdir=no > $W/maps160/docked_all/$b.log 2>&1
  echo "$b $(grep '^cc=' $W/maps160/docked_all/$b.txt)"
done; log "docking done: $(ls $W/maps160/docked_all/conf*_m160.mrc | wc -l) docked apo maps"
R=$A/rnasep_klem/f692; CONS=$W/maps160/toy_consensus.mrc
for f in 1 2; do echo "== family $f refined states vs ALL deposited apo conformers"; $PY $W/deposited_match.py $R/fam$f/states_dir $W/maps160/docked_all $CONS --save $R/fam$f/state_matches_all.txt 2>&1 | grep -v Warning | tail -60
  for d in $R/fam$f/region_1 $R/fam$f/regionmotion_1; do [ -d $d/level1/states_dir ] && { echo "== family $f $(basename $d) level-1 states vs ALL"; $PY $W/deposited_match.py $d/level1/states_dir $W/maps160/docked_all $CONS 2>&1 | grep -v Warning | grep "identified\|  [0-9]" | tail -60; }; done; done
log DONE
