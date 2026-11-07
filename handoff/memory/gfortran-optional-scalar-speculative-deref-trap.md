---
name: gfortran-optional-scalar-speculative-deref-trap
description: "TRAP (batch2 gcc 15.2 release build): an absent OPTIONAL real dummy used as `if( present(x) ) k = ... / x` was dereferenced unconditionally (if-converted select) -> NULL deref SIGSEGV in every caller omitting it; state_split/state_mergeup crashed from build 23 until em_lim became a required argument (build 32)"
metadata:
  node_type: memory
  type: project
  originSessionId: 58de2c8f-ad53-44fb-b8a2-05aadff16ec4
  modified: 2026-10-22T07:51:40.678Z
---

2026-10-22: state_split segfaulted (exit 139) within seconds in state_mergeup%new for every caller that omitted the optional
`real, optional, intent(in) :: em_lim`; site_rank (which passes em_lim=15.) worked. Bisection over snapshots: stree14-17 fine,
stree18 (build 23, the build that added em_lim) onward broken. Disassembly: the caller zeroes rcx (NULL for absent), new() saves rcx
at entry and later does `vdivss (%rax)` through it with no present() test: the compiler turned
`kb = .../EM_LIM; if( present(em_lim) ) kb = .../em_lim` into an unconditional load.
**Why:** the failure is silent at build time, data-independent, and hides behind "worked with the other program".
**How to apply:**
- Do not read an optional scalar dummy inside arithmetic, even under `if( present(...) )`. Pass the value explicitly as a required
  argument (here: `new(params, build, em_lim)` with the exported default EM_LIM), or copy it under a separate branch with no arithmetic.
- A SIGSEGV with no output and an addr2line hit on a constructor: disassemble around the fault
  (`objdump -d --start-address`) before suspecting data; `x/6i $pc` from `gdb -batch` inside sbatch gives the instruction.
Related: [[state-tree-mergeup-port-2026-10-21]].
