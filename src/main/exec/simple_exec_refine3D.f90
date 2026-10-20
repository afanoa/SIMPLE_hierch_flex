!@descr: execution of 3D refinement commanders
module simple_exec_refine3D
use simple_cmdline,              only: cmdline
use simple_string,               only: string
use simple_exec_helpers,         only: restarted_exec
use simple_commanders_mask,      only: commander_automask
use simple_commanders_volops,    only: commander_postprocess
use simple_commanders_rec,       only: commander_rec3D
use simple_commanders_refine3D,  only: commander_refine3D, commander_refine3D_auto, commander_refine3D_states, &
    &commander_refine3D_states_cavgs, commander_classify3D_refs, commander_bootstrap_rec3D
use simple_commanders_klem,      only: commander_klem_estep
use simple_commanders_klem_board, only: commander_klem_board, commander_klem_board_part
use simple_commanders_klem_tree, only: commander_klem_tree
implicit none

public :: exec_refine3D_commander
private

type(commander_automask)              :: xautomask
type(commander_postprocess)           :: xpostprocess
type(commander_rec3D)                 :: xrec3D
type(commander_bootstrap_rec3D)       :: xbootstrap_rec3D
type(commander_refine3D_auto)         :: xrefine3D_auto
type(commander_refine3D_states)       :: xrefine3D_states
type(commander_refine3D_states_cavgs) :: xrefine3D_states_cavgs
type(commander_classify3D_refs)       :: xclassify3D_refs
type(commander_refine3D)              :: xrefine3D

type(commander_klem_estep)      :: xklem_estep
type(commander_klem_board)      :: xklem_board
type(commander_klem_board_part) :: xklem_board_part
type(commander_klem_tree)       :: xklem_tree
contains

    subroutine exec_refine3D_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'automask' )
                call xautomask%execute(cline)
            case( 'postprocess' )
                call xpostprocess%execute(cline)
            case( 'reconstruct3D' )
                call xrec3D%execute(cline)
            case( 'bootstrap_rec3D' )
                call xbootstrap_rec3D%execute(cline)
            case( 'refine3D' )
                if( cline%defined('nrestarts') )then
                    call restarted_exec(cline, string('refine3D'), string('simple_exec'))
                else
                    call xrefine3D%execute(cline)
                endif
            case( 'refine3D_auto' )
                if( cline%defined('nrestarts') )then
                    call restarted_exec(cline, string('refine3D_auto'), string('simple_exec'))
                else
                    call xrefine3D_auto%execute(cline)
                endif
            case( 'refine3D_states' )
                if( cline%defined('nrestarts') )then
                    call restarted_exec(cline, string('refine3D_states'), string('simple_exec'))
                else
                    call xrefine3D_states%execute(cline)
                endif
            case( 'refine3D_states_cavgs' )
                if( cline%defined('nrestarts') )then
                    call restarted_exec(cline, string('refine3D_states_cavgs'), string('simple_exec'))
                else
                    call xrefine3D_states_cavgs%execute(cline)
                endif
            case( 'classify3D_refs' )
                if( cline%defined('nrestarts') )then
                    call restarted_exec(cline, string('classify3D_refs'), string('simple_exec'))
                else
                    call xclassify3D_refs%execute(cline)
                endif
            case( 'klem_estep' )
                call xklem_estep%execute(cline)
            case( 'klem_board' )
                call xklem_board%execute(cline)
            case( 'klem_board_part' )
                call xklem_board_part%execute(cline)
            case( 'klem_tree' )
                call xklem_tree%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_refine3D_commander

end module simple_exec_refine3D
