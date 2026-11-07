---
name: selection-without-mkdir-no-edits-a-copy
description: "TRAP: simple_exec prg=selection (and other project editors) default mkdir=yes and write the selected project into N_selection/, leaving the input projfile unchanged"
metadata:
  node_type: memory
  type: reference
  originSessionId: ea91e751-00ee-44f7-a73f-361aba072d8d
  modified: 2026-11-03T17:31:48.421Z
---

`simple_exec prg=selection projfile=X.simple ...` without `mkdir=no` creates `N_selection/X.simple` and applies the states there; `X.simple` itself keeps every particle active. Downstream runs on `X.simple` then use the unselected set, and nothing errors.

**Why:** 2026-11-03, IgM 3D tests. The prep job ran four selections (cls2D infile, ptcl2D nran, ptcl2D infile) without `mkdir=no`. All four test projects still had 2,400,639 of 2,400,639 particles active. The test masters were cancelled before any worker started.

**How to apply:** pass `mkdir=no` to selection, update_project and any other step that must edit a project in place. After each selection, check the active count (`print_project_field oritype=ptcl3D | grep -c ' state=1\.'`). Make the job exit before it submits anything when the count is wrong. Related: [[persist-project-lineage-no-new-projects]], [[check-all-program-keys-before-launch]], [[import-particles-resets-state-and-eo]].
