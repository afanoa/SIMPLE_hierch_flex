---
name: persist-project-lineage-no-new-projects
description: "run every SIMPLE step on the previous step's output project (numbered mkdir=yes dirs in one chain dir); never new_project + import_cavgs/import_particles side projects or subset scripts - user 2026-10-15: \"my directories are full of your nonsense\""
metadata:
  node_type: memory
  type: feedback
  originSessionId: 950a20f7-af3a-4db4-93f0-bde563ac38a8
  modified: 2026-10-19T14:36:31.924Z
---

2026-10-15, NOT expanded averages: I rebuilt the states3 recipe as new_project + import_cavgs on a python-cut subset
stack (34/35/36_* dirs). Those projects held averages only, no particles, so solve3D_cavgs could not map poses/states
back to the particles and nothing downstream (refine3D_states_cavgs, flex_pca on particles) could continue. User:
"so we didnt copy the last project file. you just did some nonsense and made a new project file", "my directories are
full of your nonsense. The chain should be take the selected class averages --> cls expansion --> solve3D_cavgs",
"please commit this to memory to persist project files".

**Why:** SIMPLE workflows carry state in the project file (cls2D/cls3D/ptcl2D/ptcl3D, stacks, maps). Each program run
on the previous output project keeps the particle lineage (solve3D_cavgs map2ptcls, refine3D_states_cavgs hand-back).

**How to apply:**
- Chain = one working dir, each step `simple_exec prg=<next> projfile=<previous N_prg>/<proj>.simple mkdir=yes`, giving
  1_cls_expansion, 2_solve3D_cavgs, 3_... next to each other. Copy the previous project file, never re-import.
- Subsetting goes through `selection` on that project (cls2D/ptcl states), not a python-cut stack + import.
- No new_project/import_cavgs/import_particles side projects, no extra NN_* experiment dirs, unless the user asks.
- An old recipe that used import_cavgs is not a reason to repeat it; say it breaks the lineage and use the project.
- 2026-10-19 (ribosome_10509, "make the directories clean... no nesting jobs. it should be clear what the workflow
  was"): new dataset dir = README.md (inputs, parameters + sources, build, step table, dated log) + movies.txt +
  `<proj>.simple` (new_project projname= dir=<the dir>) + flat SIMPLE step dirs run from the top + `scripts/` holding
  the sbatch scripts, their logs and diagnostics. Diagnostics run in node-local scratch, never as extra dirs.
Related: [[do-the-requested-job-not-diagnostic-sprawl]], [[keep-test-loops-minimal]].
