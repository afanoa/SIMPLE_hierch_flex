---
name: map-cavgs-selection-ranked-stack-trap
description: "map_cavgs_selection matches by image CORRELATION, so a 'top N' selection must be cut from the *_ranked.mrc stack; cutting from the unranked stack silently selects classes 1..N by index"
metadata:
  type: feedback
---

`map_cavgs_selection` (`simple_commanders_cavgs.f90:1084`) does NOT select by index. It
correlates every image in `stk2` against every image in `stk` and marks the best match:

```fortran
call calc_cartesian_corrmat(imgs_sel, imgs_all, correlations)
do isel=1,nsel
    loc = maxloc(correlations(isel,:)); states(loc(1)) = 1
end do
```

So the selection is only as good as the stack `stk2` was cut from. Cutting the first N
images out of `cavgs_iterNNN.mrc` (class order) selects **classes 1..N by index** —
quality-blind — while looking exactly like a "top N" selection on the command line.
To select by rank, cut from **`cavgs_iterNNN_ranked.mrc`**.

**Why:** this silently destroyed the EMPIAR-10028 campaign of 2026-08-21. "First 240 of
27_ranked_cavgs" was cut from the unranked stack, so the run discarded **class 285,
rank 1 of 298** (3.89 A, pop 511) plus 15 more of the top 60, and **kept 46 of the
bottom 60** including classes at 68, 80, 96 and 121 A. Rejected classes had median rank
132 where a rank-based cut would give ~268. Every flex_pca result on that project was
computed on a quality-scrambled 83,626-particle set containing genuine junk.

**How to apply:** after any `map_cavgs_selection`, verify before building on it — join the
project's cls2D `state` against `classdoc_ranked.txt` and check that the rejected classes'
ranks are the BOTTOM of the range, not the middle:

    rejected ranks: min 261 median 279 max 298   <- correct
    rejected ranks: min 1   median 132 max 298   <- silently wrong

A resolution cut is easier to sanity-check than a count: on 10028, `res <= 12 A` keeps
260 classes / 100,667 ptcls (96%), and the genuine junk is only ~16 classes at 40-121 A.
See [[flex-pca-ab-testing-traps]] for other silent invalidators.
