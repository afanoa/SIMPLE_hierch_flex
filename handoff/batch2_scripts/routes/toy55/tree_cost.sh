#!/bin/bash
# Core-hours of a run (login-node safe: sacct only): the master job ids given + every SIMPLE worker job named simple_<PROJ> that started
# inside the masters' time span (allocated cores x elapsed). usage: tree_cost.sh PROJ JOBID [JOBID ...]
PROJ=$1; shift; IDS=$(echo $* | tr ' ' ',')
read S E < <(sacct -n -X -j $IDS -o Start,End -P | awk -F'|' 'NR == 1 {s = $1; e = $2} {if ($1 < s) s = $1; if ($2 > e || $2 == "Unknown") e = $2} END {if (e == "Unknown") e = "now"; print s, e}')
M=$(sacct -n -X -j $IDS -o CPUTimeRAW -P | awk '{s += $1} END {printf "%.1f", s / 3600}')
W=$(sacct -n -X -u $USER -S $S -E $E --name=simple_$PROJ -o CPUTimeRAW,Start -P | awk -F'|' -v s=$S '$2 >= s {c += $1; n++} END {printf "%.1f %d", c / 3600, n}')
echo "run $PROJ ($IDS, $S .. $E): masters $M core-h, workers simple_$PROJ ${W% *} core-h in ${W#* } jobs"
