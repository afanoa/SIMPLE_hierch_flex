#!/bin/bash
# move the even/odd half maps of consolidation dirs out of the way: DIR/fam0/regionall/level1/maps/*_{even,odd}.mrc -> maps/halves/
# (the user opens whole folders of maps; readers that need halves look in maps/halves/ too). usage: tidy_halves.sh DIR ...
for C in "$@"; do p=$C/fam0/regionall/level1/maps; [ -d $p ] || continue; mkdir -p $p/halves
  ls $p/*_even.mrc $p/*_odd.mrc > /dev/null 2>&1 && mv -n $p/*_even.mrc $p/*_odd.mrc $p/halves/; echo "tidied $p: $(ls $p | grep -c '\.mrc$') full maps, $(ls $p/halves | wc -l) half maps in halves/"; done
