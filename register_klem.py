#!/usr/bin/env python3
"""Register the klem_estep program in a SIMPLE worktree: exec router (simple_exec_refine3D.f90) and UI
(simple_ui_refine3D.f90). Idempotent. usage: register_klem.py WORKTREE"""
import sys, re
wt = sys.argv[1]
# ---- exec router
p = wt + '/src/main/exec/simple_exec_refine3D.f90'; s = open(p).read()
if 'commander_klem_estep' not in s:
    s = s.replace("use simple_commander_refine3D_pose_cont, only: commander_refine3D_pose_cont\n", "use simple_commander_refine3D_pose_cont, only: commander_refine3D_pose_cont\nuse simple_commanders_klem,      only: commander_klem_estep\n")
    s = s.replace("type(commander_refine3D)        :: xrefine3D\n", "type(commander_refine3D)        :: xrefine3D\ntype(commander_klem_estep)      :: xklem_estep\n")
    m = re.search(r'\n(\s*)case\s*\(?\s*default\s*\)?', s, flags=re.I)
    assert m, 'no case default in the exec router'
    ind = m.group(1)
    s = s[:m.start()] + "\n%scase( 'klem_estep' )\n%s    call xklem_estep%%execute(cline)" % (ind, ind) + s[m.start():]
    open(p, 'w').write(s); print('exec router patched')
else: print('exec router already patched')
# ---- UI
p = wt + '/src/main/ui/simple/simple_ui_refine3D.f90'; s = open(p).read()
if 'klem_estep' not in s:
    s = s.replace("type(ui_program), target :: refine3D_auto\n", "type(ui_program), target :: refine3D_auto\ntype(ui_program), target :: klem_estep\n", 1)
    s = s.replace("        call new_refine3D_auto(prgtab)\n", "        call new_refine3D_auto(prgtab)\n        call new_klem_estep(prgtab)\n", 1)
    # registration line template: the last prgtab line of new_refine3D_auto
    blk = s[s.index('subroutine new_refine3D_auto'):s.index('end subroutine new_refine3D_auto')]
    reg = [l for l in blk.split('\n') if 'prgtab' in l and 'intent' not in l][-1]
    reg_k = reg.replace('refine3D_auto', 'klem_estep')
    spec = '''
    subroutine new_klem_estep( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call klem_estep%new(&
        &'klem_estep',&                                                                            ! name
        &'Kernel-latent EM E-step: cross-half log-likelihoods against cell volumes',&               ! summary
        &'is a program scoring every particle against K cell volumes of the other data half at its own pose, band-limited (the klem E-step)',& ! help
        &'simple_exec',&                                                                           ! executable
        &.true.,&                                                                                  ! requires sp_project
        &visibility=UI_VIS_ADVANCED, display_name='klem E-step')
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call klem_estep%add_input(UI_FILE, 'dir', 'dir', 'Cell volumes directory', 'Directory with cell_NN_even.mrc and cell_NN_odd.mrc on the coarse grid', 'e.g. cells_half', .true., '')
        call klem_estep%add_input(UI_FILE, 'infile', 'file', 'Noise/scale table', 'Per-shell noise variance and reference scale: shell sig2 scale', 'e.g. sig2_scale.txt', .true., '')
        call klem_estep%add_input(UI_FILE, 'outfile', 'file', 'Output stream file', 'Binary stream with the per-particle cell log-likelihoods and comparison samples', 'e.g. klem_estep.bin', .true., '')
        ! parameter input/output
        call klem_estep%add_input(UI_PARM, 'nstates', 'num', 'Number of cells', 'Number of cell volumes K', '# cells', .true., 0.0)
        call klem_estep%add_input(UI_PARM, 'lp', 'num', 'Band limit', 'Low-pass limit of the band in Angstroms', 'in A', .true., 0.0)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call klem_estep%add_input(UI_COMP, nthr)
''' + reg_k + '\n    end subroutine new_klem_estep\n'
    s = s.replace('\nend module', spec + '\nend module', 1) if '\nend module' in s else s + spec
    open(p, 'w').write(s); print('ui patched; registration line: ' + reg_k.strip())
else: print('ui already patched')

# ---- second program: klem_board (same files)
p = wt + '/src/main/exec/simple_exec_refine3D.f90'; s = open(p).read()
if 'commander_klem_board' not in s:
    i = s.index('implicit none'); s = s[:i] + 'use simple_commanders_klem_board, only: commander_klem_board\n' + s[i:]
    i = s.index('contains'); s = s[:i] + 'type(commander_klem_board)      :: xklem_board\n' + s[i:]
    m = re.search(r'\n(\s*)case\s*\(?\s*default\s*\)?', s, flags=re.I); ind = m.group(1)
    s = s[:m.start()] + "\n%scase( 'klem_board' )\n%s    call xklem_board%%execute(cline)" % (ind, ind) + s[m.start():]
    open(p, 'w').write(s); print('exec router: klem_board added')
p = wt + '/src/main/ui/simple/simple_ui_refine3D.f90'; s = open(p).read()
if 'klem_board' not in s:
    s = s.replace("type(ui_program), target :: klem_estep\n", "type(ui_program), target :: klem_estep\ntype(ui_program), target :: klem_board\n", 1)
    s = s.replace("        call new_klem_estep(prgtab)\n", "        call new_klem_estep(prgtab)\n        call new_klem_board(prgtab)\n", 1)
    reg = [l for l in s[s.index('subroutine new_klem_estep'):s.index('end subroutine new_klem_estep')].split('\n') if 'add_ui_program' in l][-1].replace('klem_estep', 'klem_board')
    spec = '''
    subroutine new_klem_board( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call klem_board%new(&
        &'klem_board',&                                                                            ! name
        &'Kernel-latent EM board: lattice of cell volumes by cross-half EM',&                      ! summary
        &'is a program fitting a lattice of free cell volumes plus junk components to particles at fixed poses by cross-half EM over a band schedule (the klem board)',& ! help
        &'simple_exec',&                                                                           ! executable
        &.true.,&                                                                                  ! requires sp_project
        &visibility=UI_VIS_ADVANCED, display_name='klem board')
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call klem_board%add_input(UI_FILE, 'infile', 'file', 'Board config', 'Config file with key value lines (n1, n2, nfree, lps, sigmas, shrink, balance, focus_mask, subset, nsub, seed, ss, fs, conv, sign, ctfabs)', 'e.g. board.txt', .true., '')
        call klem_board%add_input(UI_FILE, 'dir', 'dir', 'Output directory', 'Directory for the board outputs (klem_latent.txt, klem_states.txt, klem_resp.bin, klem_halves.bin, cells/)', 'e.g. board', .true., '')
        ! parameter input/output
        ! <empty>
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call klem_board%add_input(UI_COMP, nthr)
''' + reg + '\n    end subroutine new_klem_board\n'
    s = s.replace('\nend module', spec + '\nend module', 1)
    open(p, 'w').write(s); print('ui: klem_board added')

# ---- third program: klem_board_part (the distributed board's worker) + nparts on klem_board
p = wt + '/src/main/exec/simple_exec_refine3D.f90'; s = open(p).read()
if 'commander_klem_board_part' not in s:
    s = s.replace('use simple_commanders_klem_board, only: commander_klem_board\n', 'use simple_commanders_klem_board, only: commander_klem_board, commander_klem_board_part\n', 1)
    i = s.index('contains'); s = s[:i] + 'type(commander_klem_board_part) :: xklem_board_part\n' + s[i:]
    m = re.search(r'\n(\s*)case\s*\(?\s*default\s*\)?', s, flags=re.I); ind = m.group(1)
    s = s[:m.start()] + "\n%scase( 'klem_board_part' )\n%s    call xklem_board_part%%execute(cline)" % (ind, ind) + s[m.start():]
    open(p, 'w').write(s); print('exec router: klem_board_part added')
p = wt + '/src/main/ui/simple/simple_ui_refine3D.f90'; s = open(p).read()
if 'klem_board%add_input(UI_COMP, nparts)' not in s:
    s = s.replace('        call klem_board%add_input(UI_COMP, nthr)\n', '        call klem_board%add_input(UI_COMP, nparts)\n        call klem_board%add_input(UI_COMP, nthr)\n', 1)
    open(p, 'w').write(s); print('ui: nparts added to klem_board')
if 'klem_board_part' not in s:
    s = s.replace("type(ui_program), target :: klem_board\n", "type(ui_program), target :: klem_board\ntype(ui_program), target :: klem_board_part\n", 1)
    s = s.replace("        call new_klem_board(prgtab)\n", "        call new_klem_board(prgtab)\n        call new_klem_board_part(prgtab)\n", 1)
    reg = [l for l in s[s.index('subroutine new_klem_board('):s.index('end subroutine new_klem_board')].split('\n') if 'add_ui_program' in l][-1].replace('klem_board', 'klem_board_part')
    spec = '''
    subroutine new_klem_board_part( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call klem_board_part%new(&
        &'klem_board_part',&                                                                       ! name
        &'Kernel-latent EM board: one persistent worker of the distributed board',&                ! summary
        &'is a program holding one part of the particles for the distributed klem board and answering the master process round by round (launched by klem_board with nparts > 1)',& ! help
        &'simple_exec',&                                                                           ! executable
        &.true.,&                                                                                  ! requires sp_project
        &visibility=UI_VIS_ADVANCED, display_name='klem board worker')
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call klem_board_part%add_input(UI_FILE, 'infile', 'file', 'Board config', 'The board config file of the master', 'e.g. board.txt', .true., '')
        call klem_board_part%add_input(UI_FILE, 'dir', 'dir', 'Exchange directory', 'The master board directory holding distr/', 'e.g. board', .true., '')
        ! parameter input/output
        call klem_board_part%add_input(UI_PARM, 'part', 'num', 'Part index', 'Index of this worker (1..nparts)', '# part', .true., 1.0)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call klem_board_part%add_input(UI_COMP, nparts)
        call klem_board_part%add_input(UI_COMP, nthr)
''' + reg + '\n    end subroutine new_klem_board_part\n'
    s = s.replace('\nend module', spec + '\nend module', 1)
    open(p, 'w').write(s); print('ui: klem_board_part added')

# ---- fourth program: klem_tree (the tree driver)
p = wt + '/src/main/exec/simple_exec_refine3D.f90'; s = open(p).read()
if 'commander_klem_tree' not in s:
    i = s.index('implicit none'); s = s[:i] + 'use simple_commanders_klem_tree, only: commander_klem_tree\n' + s[i:]
    i = s.index('contains'); s = s[:i] + 'type(commander_klem_tree)       :: xklem_tree\n' + s[i:]
    m = re.search(r'\n(\s*)case\s*\(?\s*default\s*\)?', s, flags=re.I); ind = m.group(1)
    s = s[:m.start()] + "\n%scase( 'klem_tree' )\n%s    call xklem_tree%%execute(cline)" % (ind, ind) + s[m.start():]
    open(p, 'w').write(s); print('exec router: klem_tree added')
p = wt + '/src/main/ui/simple/simple_ui_refine3D.f90'; s = open(p).read()
if 'klem_tree' not in s:
    s = s.replace("type(ui_program), target :: klem_board_part\n", "type(ui_program), target :: klem_board_part\ntype(ui_program), target :: klem_tree\n", 1)
    s = s.replace("        call new_klem_board_part(prgtab)\n", "        call new_klem_board_part(prgtab)\n        call new_klem_tree(prgtab)\n", 1)
    reg = [l for l in s[s.index('subroutine new_klem_board('):s.index('end subroutine new_klem_board')].split('\n') if 'add_ui_program' in l][-1].replace('klem_board', 'klem_tree')
    spec = '''
    subroutine new_klem_tree( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call klem_tree%new(&
        &'klem_tree',&                                                                             ! name
        &'Kernel-latent EM tree: the full landscape workflow around the klem board',&              ! summary
        &'is a program running the kernel-latent EM tree (parent board, families, pose refinement, feedback boards, data-derived regions, states, level 2, landscape count) as one resumable workflow; stages run as jobs of this program',& ! help
        &'simple_exec',&                                                                           ! executable
        &.true.,&                                                                                  ! requires sp_project
        &visibility=UI_VIS_ADVANCED, display_name='klem tree')
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call klem_tree%add_input(UI_FILE, 'infile', 'file', 'Tree config', 'Config file with key value lines (run, stage, nparts, nthr_part, np_refine, maxreg, nl2, parmax, skip_refine, norec, docked_list, cons, meta_init, ...)', 'e.g. tree.txt', .true., '')
        ! parameter input/output
        ! <empty>
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call klem_tree%add_input(UI_COMP, nthr)
''' + reg + '\n    end subroutine new_klem_tree\n'
    s = s.replace('\nend module', spec + '\nend module', 1)
    open(p, 'w').write(s); print('ui: klem_tree added')

# ---- private executable: the distributed board's worker and the tree's sub-jobs run through simple_private_exec
# (the public executable appends a job record to the project at exit; concurrent worker exits corrupt it)
p = wt + '/production/simple_private_exec_driver.f90'; s = open(p).read()
if 'commander_klem_board_part' not in s:
    s = s.replace('use simple_private_exec_api\n', 'use simple_private_exec_api\nuse simple_commanders_klem_board, only: commander_klem_board_part\nuse simple_commanders_klem_tree,  only: commander_klem_tree\n', 1)
    s = s.replace('        type(commander_preprocess)              :: xpreprocess\n', '        type(commander_preprocess)              :: xpreprocess\n        type(commander_klem_board_part)         :: xklem_board_part\n        type(commander_klem_tree)               :: xklem_tree\n', 1)
    m = re.search(r'\n(\s*)case DEFAULT\s*\n\s*THROW_HARD\(.prg=', s); assert m, 'no case DEFAULT in the private driver'
    ind = m.group(1)
    s = s[:m.start()] + "\n%scase( 'klem_board_part' )\n%s    call xklem_board_part%%execute(cline)\n%scase( 'klem_tree' )\n%s    call xklem_tree%%execute(cline)" % (ind, ind, ind, ind) + s[m.start():]
    open(p, 'w').write(s); print('private exec driver: klem_board_part and klem_tree added')
else: print('private exec driver already patched')
