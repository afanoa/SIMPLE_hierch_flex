#!/bin/bash
# read-only job monitor helper (login node: squeue / sacct / scontrol / grep on logs only). env START (sacct start),
# NAMES (egrep of my orchestration job names). Emits Q| (queue), S| (finished jobs), E| (log events), A| (log age, min).
K=/mnt/beegfs/elmlund/afan/rnasep_klem; M=$K/monitor; mkdir -p $M; L=$M/logs.txt; touch $L
PATS='FAILED|ERROR STOP|Error termination|Traceback|falling back|no refinement|Killed|oom-kill|OUT_OF_MEMORY|out of memory|empty even/odd|Segmentation|DUE TO TIME|TIME LIMIT|families finished|ALL DONE|refined project|refined poses for|zooms finished|regions stage|parent subset from|== LANDSCAPE|POOLED DONE|leaf state maps|distinct conformers|CHUNK [0-9]* DONE|board server died|timed out'
squeue -u $USER -h -o '%i|%j|%T|%M|%C|%P|%r|%V' > $M/squeue.txt
awk -F'|' -v n="$NAMES" '$2 ~ n {print $1}' $M/squeue.txt | while read j; do grep -q "\.$j\.out" $L || scontrol show job $j 2>/dev/null | grep -o 'StdOut=[^ ]*' | sed 's/StdOut=//' >> $L; done
sed 's/^/Q|/' $M/squeue.txt
sacct -u $USER -S ${START:-now-1days} -X -n -P -o JobID,JobName,State,ExitCode,Elapsed 2>/dev/null | awk -F'|' -v n="$NAMES" '$2 ~ n && $3 !~ /RUNNING|PENDING/' | sed 's/^/S|/'
now=$(date +%s); sort -u $L | while read f; do [ -f "$f" ] || continue; grep -a -n -E "$PATS" "$f" 2>/dev/null | grep -a -v 'modulefile\|grep -v' | cut -c1-220 | sed "s#^#E|$f|#"; echo "A|$f|$(( (now - $(stat -c %Y "$f")) / 60 ))"; done
