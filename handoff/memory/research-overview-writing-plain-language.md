---
name: research-overview-writing-plain-language
description: "overviews of the heterogeneity work are for a non-technical audience - no specific numbers or tests, 3DVA-style plain framing, order why -> flex_pca -> flex2D -> hierarchy, hierarchy as concept only"
metadata:
  node_type: memory
  type: feedback
  originSessionId: feefe332-1f43-4b38-9d2c-f28b2e30a3db
  modified: 2026-10-19T15:31:02.083Z
---

2026-10-19, on the 1.5-page overview doc (Claude Doc "Heterogeneity analysis in SIMPLE: current work", export
doc/implementation_notes/heterogeneity_overview_2026_10_19.docx): "We dont need these specific numbers and tests",
"take some wording/ideas from cryosparc 3DVA as some of this needs to be aimed at a not as technical audience", and the
current flex model goes BEFORE flex2D. On the hierarchical work: explain only the concept (fit -> refine -> fit at each
node, because poses matter a lot), say it works well for many states and large body motions, no details (work in
progress). Streaming flex2D is partly implemented (not on the flex2D branch's src/main/stream).

**Why:** the audience is not specialist; run metrics and toy tests read as lab notes, not an overview.

**How to apply:** for overviews/summaries of the work: plain language (family of structures, modes = where density is
gained or lost, per-particle coordinate, landscape of clusters vs continuous spread), headline claims without metrics,
keep "I developed flex_pca" and the RNase P headline, without naming cryoSPARC ([[flex-pca-authorship-and-headline-result]]).

2026-10-19 (3-slide deck, heterogeneity_3slide_2026_10_19*.pptx): never name cryoSPARC, 3DVA, 3DFlex,
cryoDRGN or any other method on slides; say "existing methods". Do name flex_pca's model as probabilistic PCA
(one phrase, no maths). Cuts that lose the thread are worse than a little text: give each slide a one-line lead.
The user edits decks in PowerPoint in place: read their saved copy, edit THAT XML, write a new _vN file.
