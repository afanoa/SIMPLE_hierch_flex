!@descr: the kernel-latent EM tree as one program (prg=klem_tree): the stages of the bash chain in Fortran, resumable
!! from existing outputs. stage=all: parent board, families, parent regions, one family job per family (async),
!! the landscape. stage=family: child board, k-means, refine3D_states (unless skip_refine), feedback board, regions,
!! one region job per region (async, parmax at a time). stage=region: focal board, node states, reconstructions
!! (unless norec), level 2 along the node's own axis. stage=node: the docked/closed node enumeration. Config file
!! (infile): "key value" lines; see klem_tree_cfg below. Boards run through prg=klem_board (distributed when nparts > 1).
module simple_commanders_klem_tree
use simple_core_module_api
use simple_commander_base, only: commander_base
use simple_cmdline,        only: cmdline
use simple_parameters,     only: parameters
use simple_sp_project,     only: sp_project
use simple_string,         only: string
use simple_qsys_env,       only: qsys_env
use simple_commanders_klem_board, only: commander_klem_board
use simple_commanders_refine3D,   only: commander_refine3D_states
use simple_commanders_rec,        only: commander_rec3D
use simple_commanders_project_core, only: commander_selection
use simple_klem_tree
implicit none
public :: commander_klem_tree
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_klem_tree
  contains
    procedure :: execute => exec_klem_tree
end type commander_klem_tree

type klem_tree_cfg
    character(len=512) :: run = 'klem_tree', meta_init = 'none', subset = 'none', docked_list = 'none', cons = 'none', stage = 'all', region_mask = 'none', region_tag = 'r'
    integer :: nfam = 2, nparts = 1, nthr_part = 8, np_refine = 10, maxits = 6, nfree_parent = 2, maxreg = 6, nl2 = 3, parmax = 2, kgroups = 5, nmax = 30000, n2 = 10000, maxk = 8, fam = 1, nsub = 0
    integer :: skip_refine = 0, norec = 0, nprec = 4, qmax = 4
    character(len=8) :: rank = 'pc1'
    real    :: balance = 0.5, devw = 0.05, sim = 0.7, frac = 0.8, minmass = 2000., mskdiam = 232., mskdiam_refine = 220., region_thr = 0.5, cmerge = 0.8
    character(len=256) :: lps_parent = '20 20 20 14 14 14 10 10 10', sig_parent = '1.5 1.5 1.2 1.0 1.0 0.8 0.7 0.7 0.7'
    character(len=256) :: lps_child = '14 14 12 12 10 10 8 8 8', sig_child = '1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7'
    character(len=256) :: lps_region = '12 12 10 10 8 8 8 8 8', sig_region = '1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7'
end type klem_tree_cfg

contains

    subroutine exec_klem_tree( self, cline )
        class(commander_klem_tree), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(parameters)    :: params
        type(klem_tree_cfg) :: cfg
        type(sp_project)    :: spproj
        type(string)        :: cwd0
        character(len=512)  :: R, src, line, tok(64)
        integer :: u, io, ntok, box, nrows
        real    :: smpd
        call params%new(cline)
        call read_tree_cfg(params%infile%to_char(), cfg)
        call simple_getcwd(cwd0); R = trim(cfg%run); call simple_mkdir(trim(R))
        src = params%projfile%to_char(); box = params%box; smpd = params%smpd
        call spproj%read(params%projfile); nrows = spproj%os_ptcl3D%get_noris(); call spproj%kill
        write(logfhandle,'(A)') '>>> KLEM TREE stage '//trim(cfg%stage)//' in '//trim(R)//' on '//trim(src)
        select case(trim(cfg%stage))
            case('all');    call stage_all(params, cfg, R, src, box, smpd, nrows)
            case('family'); call stage_family(params, cfg, R, src, box, smpd, nrows, cfg%fam)
            case('region'); call stage_region(params, cfg, R, src, box, smpd, nrows, cfg%fam, trim(cfg%region_mask), trim(cfg%region_tag))
            case('node');   call stage_node(params, cfg, R, src, box, smpd, nrows, cfg%fam)
            case default;   THROW_HARD('unknown stage '//trim(cfg%stage))
        end select
        call simple_chdir(cwd0)
        call simple_end('**** SIMPLE_KLEM_TREE NORMAL STOP ****')
    end subroutine exec_klem_tree

    ! ------------------------------------------------------------------ config
    subroutine read_tree_cfg( fname, c )
        character(len=*),    intent(in)    :: fname
        type(klem_tree_cfg), intent(inout) :: c
        character(len=1024) :: line
        character(len=512)  :: key, val
        integer :: u, io, p
        open(newunit=u, file=trim(fname), status='old', action='read')
        do
            read(u,'(A)',iostat=io) line
            if( io /= 0 ) exit
            if( len_trim(line) == 0 .or. line(1:1) == '#' ) cycle
            p = index(trim(line), ' '); if( p <= 0 ) cycle
            key = line(1:p-1); val = adjustl(line(p+1:))
            select case(trim(key))
                case('run');            c%run = val
                case('meta_init');      c%meta_init = val
                case('subset');         c%subset = val
                case('docked_list');    c%docked_list = val
                case('cons');           c%cons = val
                case('stage');          c%stage = val
                case('region_mask');    c%region_mask = val
                case('region_tag');     c%region_tag = val
                case('nfam');           read(val,*) c%nfam
                case('nparts');         read(val,*) c%nparts
                case('nthr_part');      read(val,*) c%nthr_part
                case('np_refine');      read(val,*) c%np_refine
                case('maxits');         read(val,*) c%maxits
                case('nfree_parent');   read(val,*) c%nfree_parent
                case('maxreg');         read(val,*) c%maxreg
                case('nl2');            read(val,*) c%nl2
                case('parmax');         read(val,*) c%parmax
                case('kgroups');        read(val,*) c%kgroups
                case('nmax');           read(val,*) c%nmax
                case('n2');             read(val,*) c%n2
                case('maxk');           read(val,*) c%maxk
                case('fam');            read(val,*) c%fam
                case('nsub');           read(val,*) c%nsub
                case('skip_refine');    read(val,*) c%skip_refine
                case('norec');          read(val,*) c%norec
                case('rank');           c%rank = val
                case('qmax');           read(val,*) c%qmax
                case('nprec');          read(val,*) c%nprec
                case('balance');        read(val,*) c%balance
                case('devw');           read(val,*) c%devw
                case('sim');            read(val,*) c%sim
                case('frac');           read(val,*) c%frac
                case('minmass');        read(val,*) c%minmass
                case('mskdiam');        read(val,*) c%mskdiam
                case('mskdiam_refine'); read(val,*) c%mskdiam_refine
                case('region_thr');     read(val,*) c%region_thr
                case('cmerge');         read(val,*) c%cmerge
                case('lps_parent');     c%lps_parent = val
                case('sig_parent');     c%sig_parent = val
                case('lps_child');      c%lps_child = val
                case('sig_child');      c%sig_child = val
                case('lps_region');     c%lps_region = val
                case('sig_region');     c%sig_region = val
            end select
        enddo
        close(u)
    end subroutine read_tree_cfg

    ! ------------------------------------------------------------------ one board (prg=klem_board in-process)
    subroutine run_board( params, cfg, out, proj, subset, lps, sigmas, nfree, focus, init_file, pose_every )
        class(parameters),   intent(inout) :: params
        type(klem_tree_cfg), intent(in)    :: cfg
        character(len=*),    intent(in)    :: out, proj, subset, lps, sigmas, focus, init_file
        integer,             intent(in)    :: nfree, pose_every
        type(commander_klem_board) :: xboard
        type(cmdline) :: cl
        type(string)  :: s
        integer :: u
        if( file_exists(trim(out)//'/klem_latent.txt') )then
            write(logfhandle,'(A)') '>>> board exists, reused: '//trim(out); return
        endif
        call simple_mkdir(trim(out))
        open(newunit=u, file=trim(out)//'/board.txt', status='replace', action='write')
        write(u,'(A)') 'n1 7'; write(u,'(A)') 'n2 7'; write(u,'(A,I0)') 'nfree ', nfree
        write(u,'(A)') 'lps '//trim(lps); write(u,'(A)') 'sigmas '//trim(sigmas)
        write(u,'(A)') 'shrink 0.02'; write(u,'(A)') 'eps0 0.05'; write(u,'(A,F6.3)') 'balance ', cfg%balance
        write(u,'(A)') 'subset '//trim(subset); write(u,'(A)') 'nsub 0'; write(u,'(A)') 'seed 1'
        write(u,'(A)') 'conv 0'; write(u,'(A)') 'ss -1'; write(u,'(A)') 'fs 1'; write(u,'(A)') 'sign 1'; write(u,'(A)') 'ctfabs 0'; write(u,'(A)') 'write_halves 0'
        write(u,'(A,F6.3)') 'devw ', cfg%devw; write(u,'(A)') 'pnoise 1'
        if( trim(focus) /= 'none' ) write(u,'(A)') 'focus_mask '//trim(focus)
        if( trim(init_file) /= 'none' ) write(u,'(A)') 'init_file '//trim(init_file)
        if( pose_every > 0 )then
            write(u,'(A,I0)') 'pose_every ', pose_every; write(u,'(A)') 'pose_step 3'; write(u,'(A)') 'pose_prior 4'
        endif
        if( cfg%nparts > 1 ) write(u,'(A,I0)') 'nthr_part ', cfg%nthr_part
        close(u)
        call cl%set('prg', 'klem_board'); s = trim(proj); call cl%set('projfile', s); s = trim(out)//'/board.txt'; call cl%set('infile', s)
        s = trim(out); call cl%set('dir', s); call cl%set('nthr', merge(cfg%nthr_part, params%nthr, cfg%nparts > 1)); call cl%set('mkdir', 'no')
        if( cfg%nparts > 1 ) call cl%set('nparts', cfg%nparts)
        write(logfhandle,'(A)') '>>> board -> '//trim(out); flush(logfhandle)
        call xboard%execute(cl)
    end subroutine run_board

    !> the project's rows that belong to a subset file (1-based) as a states file (nrows lines, state or 0)
    subroutine write_states_from_subset( subset, nrows, fname, state )
        character(len=*), intent(in) :: subset, fname
        integer,          intent(in) :: nrows, state
        integer, allocatable :: st(:)
        integer :: u, io, i
        allocate(st(nrows), source=0)
        open(newunit=u, file=trim(subset), status='old', action='read')
        do
            read(u,*,iostat=io) i
            if( io /= 0 ) exit
            if( i >= 1 .and. i <= nrows ) st(i) = state
        enddo
        close(u)
        open(newunit=u, file=trim(fname), status='replace', action='write')
        do i = 1,nrows
            write(u,'(I0)') st(i)
        enddo
        close(u); deallocate(st)
    end subroutine write_states_from_subset

    !> selection + reconstruct3D of the states in a states file at the poses of PROJ, maps copied to OUT/maps
    subroutine reconstruct_states( params, cfg, proj, states_file, nstates, out )
        class(parameters),   intent(inout) :: params
        type(klem_tree_cfg), intent(in)    :: cfg
        character(len=*),    intent(in)    :: proj, states_file, out
        integer,             intent(in)    :: nstates
        type(commander_selection) :: xsel
        type(commander_rec3D)     :: xrec
        type(cmdline) :: cl
        type(string)  :: s, cwd
        integer :: st
        character(len=512) :: rdir
        if( file_exists(trim(out)//'/maps/recvol_state01.mrc') )then
            write(logfhandle,'(A)') '>>> reconstructions exist, reused: '//trim(out); return
        endif
        rdir = trim(out)//'/rec'; call simple_mkdir(trim(rdir)); call simple_mkdir(trim(out)//'/maps'); call simple_mkdir(trim(out)//'/states_dir')
        call simple_getcwd(cwd); call simple_chdir(trim(rdir))
        call copyf(trim(proj), trim(rdir)//'/p.simple')
        call cl%set('prg', 'selection'); s = 'p.simple'; call cl%set('projfile', s); call cl%set('oritype', 'ptcl3D'); s = trim(states_file); call cl%set('infile', s); call cl%set('mkdir', 'no')
        call xsel%execute(cl); call cl%kill
        call cl%set('prg', 'reconstruct3D'); s = 'p.simple'; call cl%set('projfile', s); call cl%set('pgrp', 'c1'); call cl%set('mskdiam', cfg%mskdiam); call cl%set('nstates', nstates)
        call cl%set('nparts', cfg%nprec); call cl%set('nthr', 8); call cl%set('mkdir', 'yes')
        call xrec%execute(cl); call cl%kill
        call simple_chdir(cwd)
        do st = 1,nstates
            call copy_latest(trim(rdir), 'recvol_state'//int2str_pad(st,2)//'.mrc', trim(out)//'/maps/recvol_state'//int2str_pad(st,2)//'.mrc')
            call copy_latest(trim(rdir), 'recvol_state'//int2str_pad(st,2)//'_even.mrc', trim(out)//'/maps/recvol_state'//int2str_pad(st,2)//'_even.mrc')
            call copy_latest(trim(rdir), 'recvol_state'//int2str_pad(st,2)//'_odd.mrc', trim(out)//'/maps/recvol_state'//int2str_pad(st,2)//'_odd.mrc')
            call copy_latest(trim(rdir), 'recvol_state'//int2str_pad(st,2)//'.mrc', trim(out)//'/states_dir/state_'//int2str_pad(st,2)//'.mrc')
        enddo
    end subroutine reconstruct_states

    !> copy FILE from the most recent *_reconstruct3D subdirectory of DIR (as the chain's ls -td | head -1)
    subroutine copy_latest( dir, file, dest )
        character(len=*), intent(in) :: dir, file, dest
        type(string), allocatable :: list(:)
        integer :: i, best
        character(len=512) :: cand
        call simple_list_files(trim(dir)//'/*_reconstruct3D', list)
        if( .not. allocated(list) ) return
        best = 0
        do i = 1,size(list)
            cand = list(i)%to_char()
            if( file_exists(trim(cand)//'/'//trim(file)) ) best = i
        enddo
        if( best > 0 ) call copyf(list(best)%to_char()//'/'//trim(file), trim(dest))
    end subroutine copy_latest

    ! ------------------------------------------------------------------ stage: all
    subroutine stage_all( params, cfg, R, src, box, smpd, nrows )
        class(parameters),   intent(inout) :: params
        type(klem_tree_cfg), intent(in)    :: cfg
        character(len=*),    intent(in)    :: R, src
        integer,             intent(in)    :: box, nrows
        real,                intent(in)    :: smpd
        type(board_t)   :: B
        type(sp_project):: spproj
        type(qsys_env)  :: qenv
        type(cmdline)   :: cl
        type(string)    :: eb, scr, outf
        integer, allocatable :: nmem(:)
        character(len=512), allocatable :: leaves(:)
        character(len=512) :: psub
        integer :: u, i, j, nreg, nleaf, waited
        real    :: sims(4)
        ! (1) parent subset (all particles with state > 0, or the given subset) and the parent board
        if( .not. file_exists(trim(R)//'/subset_all.txt') )then
            call spproj%read(params%projfile)
            open(newunit=u, file=trim(R)//'/subset_all.txt', status='replace', action='write')
            do i = 1,nrows
                if( spproj%os_ptcl3D%get_state(i) > 0 ) write(u,'(I0)') i
            enddo
            close(u); call spproj%kill
        endif
        psub = trim(R)//'/subset_all.txt'; if( trim(cfg%subset) /= 'none' ) psub = trim(cfg%subset)
        call run_board(params, cfg, trim(R)//'/parent', src, trim(psub), cfg%lps_parent, cfg%sig_parent, cfg%nfree_parent, 'none', cfg%meta_init, 0)
        call read_board(trim(R)//'/parent', B)
        if( trim(cfg%docked_list) /= 'none' ) call tree_deposited_match(B, box, smpd, cfg%docked_list, cfg%cons, 60., 20., trim(R)//'/parent/cell_matches.txt')
        ! (2) families and the parent's own moving region
        call simple_mkdir(trim(R)//'/expand'); allocate(nmem(cfg%nfam))
        if( .not. file_exists(trim(R)//'/expand/family1_subset.txt') ) call tree_families(B, cfg%nfam, trim(R)//'/expand', nmem)
        if( .not. file_exists(trim(R)//'/regions_parent/region_1.mrc') ) call tree_regions(B, box, smpd, 1, cfg%region_thr, 40, 2, 2, cfg%kgroups, cfg%cmerge, 0.5, trim(R)//'/regions_parent', nreg)
        call kill_board(B)
        ! (3) the families as concurrent jobs of this program (stage family)
        call qenv%new(params, cfg%nfam, qsys_nthr=max(cfg%nthr_part, 8)); call simple_mkdir('stderrout')   ! sub-jobs via simple_private_exec (no project job records)
        do j = 1,cfg%nfam
            if( file_exists(trim(R)//'/fam'//int2str(j)//'/DONE') ) cycle
            call write_stage_cfg(cfg, trim(R)//'/tree_fam'//int2str(j)//'.txt', 'family', j, 'none', 'r')
            call cl%kill; call cl%set('prg', 'klem_tree'); call cl%set('projfile', params%projfile); scr = trim(R)//'/tree_fam'//int2str(j)//'.txt'; call cl%set('infile', scr); call cl%set('nthr', max(cfg%nthr_part, 8)); call cl%set('mkdir', 'no')
            scr = 'klem_tree_fam'//int2str(j); outf = 'klem_tree_fam'//int2str(j)//'.log'; call qenv%exec_simple_prg_in_queue_async(cl, scr, outf)
            write(logfhandle,'(A,I0,A)') '>>> family ', j, ' submitted as a job'; flush(logfhandle)
        enddo
        waited = 0
        do
            if( all([(file_exists(trim(R)//'/fam'//int2str(j)//'/DONE'), j=1,cfg%nfam)]) ) exit
            call sleep(60); waited = waited + 60
            if( waited > 86400 ) THROW_HARD('family jobs did not finish within a day')
        enddo
        write(logfhandle,'(A)') '>>> families finished'
        ! (4) the landscape over the leaves (feedback boards + region boards)
        nleaf = 0; allocate(leaves(2 * cfg%nfam * (cfg%maxreg + 1)))
        do j = 1,cfg%nfam
            if( file_exists(trim(R)//'/fam'//int2str(j)//'/klem_fb/klem_halves.bin') )then
                nleaf = nleaf + 1; leaves(nleaf) = trim(R)//'/fam'//int2str(j)//'/klem_fb'
            endif
            do i = 1,cfg%maxreg
                if( file_exists(trim(R)//'/fam'//int2str(j)//'/region_'//int2str(i)//'/board/klem_halves.bin') )then
                    nleaf = nleaf + 1; leaves(nleaf) = trim(R)//'/fam'//int2str(j)//'/region_'//int2str(i)//'/board'
                endif
                if( file_exists(trim(R)//'/fam'//int2str(j)//'/regionmotion_'//int2str(i)//'/board/klem_halves.bin') )then
                    nleaf = nleaf + 1; leaves(nleaf) = trim(R)//'/fam'//int2str(j)//'/regionmotion_'//int2str(i)//'/board'
                endif
            enddo
        enddo
        sims = [0.6, 0.7, 0.8, 0.9]
        if( nleaf > 0 ) call tree_count_states(leaves(1:nleaf), nleaf, 0.5, sims, 4, box)
        write(logfhandle,'(A)') '>>> ALL DONE ('//trim(R)//')'
        call qenv%kill
    end subroutine stage_all

    !> a stage config = the run config with stage/fam/region set
    subroutine write_stage_cfg( cfg, fname, stage, fam, mask, tag )
        type(klem_tree_cfg), intent(in) :: cfg
        character(len=*),    intent(in) :: fname, stage, mask, tag
        integer,             intent(in) :: fam
        integer :: u
        open(newunit=u, file=trim(fname), status='replace', action='write')
        write(u,'(A)') 'run '//trim(cfg%run); write(u,'(A)') 'stage '//trim(stage); write(u,'(A,I0)') 'fam ', fam
        write(u,'(A)') 'region_mask '//trim(mask); write(u,'(A)') 'region_tag '//trim(tag)
        write(u,'(A)') 'meta_init '//trim(cfg%meta_init); write(u,'(A)') 'docked_list '//trim(cfg%docked_list); write(u,'(A)') 'cons '//trim(cfg%cons)
        write(u,'(A,I0)') 'nfam ', cfg%nfam; write(u,'(A,I0)') 'nparts ', cfg%nparts; write(u,'(A,I0)') 'nthr_part ', cfg%nthr_part; write(u,'(A,I0)') 'np_refine ', cfg%np_refine
        write(u,'(A,I0)') 'maxits ', cfg%maxits; write(u,'(A,I0)') 'maxreg ', cfg%maxreg; write(u,'(A,I0)') 'nl2 ', cfg%nl2; write(u,'(A,I0)') 'parmax ', cfg%parmax; write(u,'(A,I0)') 'kgroups ', cfg%kgroups
        write(u,'(A,I0)') 'nmax ', cfg%nmax; write(u,'(A,I0)') 'n2 ', cfg%n2; write(u,'(A,I0)') 'maxk ', cfg%maxk; write(u,'(A,I0)') 'skip_refine ', cfg%skip_refine; write(u,'(A,I0)') 'norec ', cfg%norec; write(u,'(A,I0)') 'nprec ', cfg%nprec; write(u,'(A)') 'rank '//trim(cfg%rank); write(u,'(A,I0)') 'qmax ', cfg%qmax
        write(u,'(A,F7.3)') 'balance ', cfg%balance; write(u,'(A,F7.3)') 'devw ', cfg%devw; write(u,'(A,F7.3)') 'sim ', cfg%sim; write(u,'(A,F7.3)') 'frac ', cfg%frac; write(u,'(A,F9.1)') 'minmass ', cfg%minmass
        write(u,'(A,F7.1)') 'mskdiam ', cfg%mskdiam; write(u,'(A,F7.1)') 'mskdiam_refine ', cfg%mskdiam_refine; write(u,'(A,F6.2)') 'region_thr ', cfg%region_thr; write(u,'(A,F6.2)') 'cmerge ', cfg%cmerge
        write(u,'(A)') 'lps_parent '//trim(cfg%lps_parent); write(u,'(A)') 'sig_parent '//trim(cfg%sig_parent); write(u,'(A)') 'lps_child '//trim(cfg%lps_child); write(u,'(A)') 'sig_child '//trim(cfg%sig_child)
        write(u,'(A)') 'lps_region '//trim(cfg%lps_region); write(u,'(A)') 'sig_region '//trim(cfg%sig_region)
        close(u)
    end subroutine write_stage_cfg

    ! ------------------------------------------------------------------ stage: family
    subroutine stage_family( params, cfg, R, src, box, smpd, nrows, fam )
        class(parameters),   intent(inout) :: params
        type(klem_tree_cfg), intent(in)    :: cfg
        character(len=*),    intent(in)    :: R, src
        integer,             intent(in)    :: box, nrows, fam
        real,                intent(in)    :: smpd
        type(board_t) :: B
        type(commander_selection)       :: xsel
        type(commander_refine3D_states) :: xref
        type(qsys_env) :: qenv
        type(cmdline)  :: cl
        type(string)   :: s, cwd, eb, scr, outf
        character(len=512) :: D, sub, ps, focus
        character(len=512), allocatable :: masks(:), tags(:)
        real,    allocatable :: X(:,:), cent(:,:)
        integer, allocatable :: lab(:), st(:)
        integer :: u, i, j, k, non, nreg, nreg2, nm, nrun, waited, K10
        D = trim(R)//'/fam'//int2str(fam); sub = trim(R)//'/expand/family'//int2str(fam)//'_subset.txt'; call simple_mkdir(trim(D))
        ! (1) child board at the consensus poses, focus = the parent's own moving region
        focus = 'none'; if( file_exists(trim(R)//'/regions_parent/region_1.mrc') ) focus = trim(R)//'/regions_parent/region_1.mrc'
        call run_board(params, cfg, trim(D)//'/child', src, sub, cfg%lps_child, cfg%sig_child, 0, focus, cfg%meta_init, 0)
        call read_board(trim(D)//'/child', B)
        if( trim(cfg%docked_list) /= 'none' ) call tree_deposited_match(B, box, smpd, cfg%docked_list, cfg%cons, 60., 20., trim(D)//'/child/cell_matches.txt')
        ! (2) k-means K = 10 on the child's posterior means -> states.txt
        K10 = 10
        if( .not. file_exists(trim(D)//'/states.txt') )then
            non = count(B%gm > 0.5); allocate(X(non,2), lab(non), cent(K10,2)); j = 0
            do i = 1,B%n
                if( B%gm(i) > 0.5 )then
                    j = j + 1; X(j,1) = B%z1(i); X(j,2) = B%z2(i)
                endif
            enddo
            call tree_kmeans(X, non, 2, K10, 10, lab, cent); allocate(st(nrows), source=0); j = 0
            do i = 1,B%n
                if( B%gm(i) > 0.5 )then
                    j = j + 1; st(B%pinds(i)) = lab(j)
                endif
            enddo
            open(newunit=u, file=trim(D)//'/states.txt', status='replace', action='write')
            do i = 1,nrows
                write(u,'(I0)') st(i)
            enddo
            close(u); deallocate(X, lab, cent, st)
            write(logfhandle,'(A,I0,A,I0,A)') '>>> family ', fam, ': k-means K=', K10, ' on the child board posterior means'
        endif
        ! (3) refine3D_states (or the consensus poses as the refined project)
        ps = trim(D)//'/refine/1_refine3D_states/p.simple'
        if( .not. file_exists(trim(ps)) )then
            call simple_mkdir(trim(D)//'/refine'); call copyf(trim(src), trim(D)//'/refine/p.simple'); call simple_getcwd(cwd); call simple_chdir(trim(D)//'/refine')
            call cl%set('prg', 'selection'); s = 'p.simple'; call cl%set('projfile', s); call cl%set('oritype', 'ptcl3D'); s = trim(D)//'/states.txt'; call cl%set('infile', s); call cl%set('mkdir', 'no')
            call xsel%execute(cl); call cl%kill
            if( cfg%skip_refine /= 0 )then
                call simple_mkdir(trim(D)//'/refine/1_refine3D_states'); call copyf(trim(D)//'/refine/p.simple', trim(ps))
                write(logfhandle,'(A)') '>>> skip_refine: the consensus poses stand in as the refined project'
            else
                call cl%set('prg', 'refine3D_states'); s = 'p.simple'; call cl%set('projfile', s); call cl%set('nstates', K10); call cl%set('pose_policy', 'local')
                call cl%set('lpstart', 20.); call cl%set('lpstop', 8.); call cl%set('maxits', cfg%maxits); call cl%set('mskdiam', cfg%mskdiam_refine); call cl%set('pgrp', 'c1')
                call cl%set('nparts', cfg%np_refine); call cl%set('nthr', 8)
                write(logfhandle,'(A,I0,A,I0,A)') '>>> refine3D_states K=', K10, ' local 20->8 A maxits ', cfg%maxits, ' (distributed)'; flush(logfhandle)
                call xref%execute(cl); call cl%kill
            endif
            call simple_chdir(cwd)
        endif
        if( .not. file_exists(trim(ps)) ) THROW_HARD('no refined project: '//trim(ps))
        ! (4) feedback board at the refined poses, focus = the child board's own moving region
        if( .not. file_exists(trim(D)//'/regions_child/region_1.mrc') ) call tree_regions(B, box, smpd, 1, cfg%region_thr, 40, 2, 2, cfg%kgroups, cfg%cmerge, 0.5, trim(D)//'/regions_child', nreg)
        call kill_board(B)
        focus = 'none'; if( file_exists(trim(D)//'/regions_child/region_1.mrc') ) focus = trim(D)//'/regions_child/region_1.mrc'
        call run_board(params, cfg, trim(D)//'/klem_fb', ps, sub, cfg%lps_child, cfg%sig_child, 0, focus, cfg%meta_init, 0)
        call read_board(trim(D)//'/klem_fb', B)
        if( trim(cfg%docked_list) /= 'none' ) call tree_deposited_match(B, box, smpd, cfg%docked_list, cfg%cons, 60., 20., trim(D)//'/klem_fb/cell_matches.txt')
        ! (5) regions of the feedback board: spatial (up to 2) then motion groups
        if( .not. file_exists(trim(D)//'/regions/region_1.mrc') ) call tree_regions(B, box, smpd, 1, cfg%region_thr, 40, 2, 2, cfg%kgroups, cfg%cmerge, 0.5, trim(D)//'/regions', nreg)
        if( .not. file_exists(trim(D)//'/regions_motion/region_1.mrc') ) call tree_regions(B, box, smpd, 2, cfg%region_thr, 40, 2, 6, cfg%kgroups, cfg%cmerge, 0.5, trim(D)//'/regions_motion', nreg2)
        call kill_board(B)
        allocate(masks(cfg%maxreg), tags(cfg%maxreg)); nm = 0
        do i = 1,2
            if( nm < cfg%maxreg .and. file_exists(trim(D)//'/regions/region_'//int2str(i)//'.mrc') )then
                nm = nm + 1; masks(nm) = trim(D)//'/regions/region_'//int2str(i)//'.mrc'; tags(nm) = '_'//int2str(i)
            endif
        enddo
        do i = 1,6
            if( nm < cfg%maxreg .and. file_exists(trim(D)//'/regions_motion/region_'//int2str(i)//'.mrc') )then
                nm = nm + 1; masks(nm) = trim(D)//'/regions_motion/region_'//int2str(i)//'.mrc'; tags(nm) = 'motion_'//int2str(i)
            endif
        enddo
        ! (6) region jobs, parmax at a time
        call qenv%new(params, max(cfg%parmax,1), qsys_nthr=max(cfg%nthr_part, 8)); call simple_mkdir('stderrout'); nrun = 0
        do i = 1,nm
            if( file_exists(trim(D)//'/region'//trim(tags(i))//'/DONE') ) cycle
            do while( count_running_regions(D, masks, tags, nm) >= cfg%parmax )
                call sleep(60)
            enddo
            call write_stage_cfg(cfg, trim(D)//'/tree_region'//trim(tags(i))//'.txt', 'region', fam, trim(masks(i)), trim(tags(i)))
            call cl%kill; call cl%set('prg', 'klem_tree'); call cl%set('projfile', params%projfile); s = trim(D)//'/tree_region'//trim(tags(i))//'.txt'; call cl%set('infile', s); call cl%set('nthr', max(cfg%nthr_part, 8)); call cl%set('mkdir', 'no')
            call simple_mkdir(trim(D)//'/region'//trim(tags(i))); call touch(trim(D)//'/region'//trim(tags(i))//'/RUNNING')
            scr = 'klem_tree_fam'//int2str(fam)//'_region'//trim(tags(i)); outf = trim(scr%to_char())//'.log'; call qenv%exec_simple_prg_in_queue_async(cl, scr, outf)
            write(logfhandle,'(A)') '>>> region '//trim(tags(i))//' submitted as a job'; flush(logfhandle)
        enddo
        waited = 0
        do
            if( all([(file_exists(trim(D)//'/region'//trim(tags(i))//'/DONE'), i=1,nm)]) ) exit
            call sleep(60); waited = waited + 60
            if( waited > 43200 ) exit
        enddo
        call touch(trim(D)//'/DONE'); write(logfhandle,'(A,I0,A)') '>>> family ', fam, ' DONE'
        call qenv%kill
    end subroutine stage_family

    function count_running_regions( D, masks, tags, nm ) result( c )
        character(len=*), intent(in) :: D, masks(:), tags(:)
        integer,          intent(in) :: nm
        integer :: c, i
        c = 0
        do i = 1,nm
            if( file_exists(trim(D)//'/region'//trim(tags(i))//'/RUNNING') .and. .not. file_exists(trim(D)//'/region'//trim(tags(i))//'/DONE') ) c = c + 1
        enddo
    end function count_running_regions

    ! ------------------------------------------------------------------ stage: region
    subroutine stage_region( params, cfg, R, src, box, smpd, nrows, fam, mask, tag )
        class(parameters),   intent(inout) :: params
        type(klem_tree_cfg), intent(in)    :: cfg
        character(len=*),    intent(in)    :: R, src, mask, tag
        integer,             intent(in)    :: box, nrows, fam
        real,                intent(in)    :: smpd
        type(board_t) :: B
        character(len=512) :: D, RB, sub, ps, L
        integer :: nst, s, n, n2, u, io, i
        D = trim(R)//'/fam'//int2str(fam); RB = trim(D)//'/region'//trim(tag); sub = trim(R)//'/expand/family'//int2str(fam)//'_subset.txt'
        ps = trim(D)//'/refine/1_refine3D_states/p.simple'; call simple_mkdir(trim(RB))
        ! (1) focal board at the refined poses
        call run_board(params, cfg, trim(RB)//'/board', ps, sub, cfg%lps_region, cfg%sig_region, 0, mask, cfg%meta_init, 0)
        call read_board(trim(RB)//'/board', B)
        if( trim(cfg%docked_list) /= 'none' ) call tree_deposited_match(B, box, smpd, cfg%docked_list, cfg%cons, 60., 20., trim(RB)//'/board/cell_matches.txt')
        ! (2) states from the reproducible cell clusters
        call tree_node_states(B, nrows, cfg%sim, cfg%nmax, cfg%frac, cfg%minmass, cfg%maxk, 0.5, trim(RB)//'/level1', nst)
        call kill_board(B)
        if( nst < 1 .or. cfg%norec /= 0 )then
            call touch(trim(RB)//'/DONE'); return
        endif
        call reconstruct_states(params, cfg, ps, trim(RB)//'/level1/states_clusters.txt', nst, trim(RB)//'/level1')
        ! (3) level 2 for the largest states: within-state board with the inherited region, ranked along its own axis
        do s = 1,min(nst, cfg%nl2)
            L = trim(RB)//'/level2_s'//int2str(s); call simple_mkdir(trim(L))
            ! subset = rows with state s
            open(newunit=u, file=trim(RB)//'/level1/states_clusters.txt', status='old', action='read'); n = 0
            open(newunit=io, file=trim(L)//'/subset.txt', status='replace', action='write')
            do i = 1,nrows
                read(u,*) n2
                if( n2 == s )then
                    write(io,'(I0)') i; n = n + 1
                endif
            enddo
            close(u); close(io)
            if( n < 4000 ) cycle
            call run_board(params, cfg, trim(L)//'/board', ps, trim(L)//'/subset.txt', cfg%lps_region, cfg%sig_region, 0, mask, cfg%meta_init, 0)
            call read_board(trim(L)//'/board', B); n2 = min(cfg%n2, n * 4 / 10); if( trim(cfg%rank) == 'ppca' )then; call tree_rank_ppca(B, n2, nrows, trim(L), cfg%qmax); else; call tree_rank_pc1(B, n2, nrows, trim(L)); endif; call kill_board(B)
            call reconstruct_states(params, cfg, ps, trim(L)//'/states_top.txt', 2, trim(L))
        enddo
        call touch(trim(RB)//'/DONE')
    end subroutine stage_region

    ! ------------------------------------------------------------------ stage: node (the docked node as a multi-state node)
    subroutine stage_node( params, cfg, R, src, box, smpd, nrows, fam )
        class(parameters),   intent(inout) :: params
        type(klem_tree_cfg), intent(in)    :: cfg
        character(len=*),    intent(in)    :: R, src
        integer,             intent(in)    :: box, nrows, fam
        real,                intent(in)    :: smpd
        type(board_t) :: B
        character(len=512) :: D, O, ps
        integer :: nst
        D = trim(R)//'/fam'//int2str(fam); O = trim(D)//'/node'; ps = trim(D)//'/refine/1_refine3D_states/p.simple'; call simple_mkdir(trim(O))
        if( trim(cfg%subset) == 'none' .or. trim(cfg%region_mask) == 'none' ) THROW_HARD('stage node needs subset (the node particles) and region_mask (the inherited region)')
        call run_board(params, cfg, trim(O)//'/board', ps, cfg%subset, cfg%lps_region, cfg%sig_region, 0, cfg%region_mask, cfg%meta_init, 0)
        call read_board(trim(O)//'/board', B)
        if( trim(cfg%docked_list) /= 'none' ) call tree_deposited_match(B, box, smpd, cfg%docked_list, cfg%cons, 60., 20., trim(O)//'/board/cell_matches.txt')
        call tree_node_states(B, nrows, cfg%sim, cfg%nmax, cfg%frac, cfg%minmass, cfg%maxk, 0.5, trim(O)//'/states', nst); call kill_board(B)
        if( nst >= 1 .and. cfg%norec == 0 ) call reconstruct_states(params, cfg, ps, trim(O)//'/states/states_clusters.txt', nst, trim(O)//'/states')
        call touch(trim(O)//'/DONE')
    end subroutine stage_node

    subroutine copyf( a, b )
        character(len=*), intent(in) :: a, b
        type(string) :: sa, sb
        sa = trim(a); sb = trim(b); call simple_copy_file(sa, sb)
    end subroutine copyf

    subroutine touch( fname )
        character(len=*), intent(in) :: fname
        integer :: u
        open(newunit=u, file=trim(fname), status='replace', action='write'); write(u,'(A)') 'ok'; close(u)
    end subroutine touch

end module simple_commanders_klem_tree
