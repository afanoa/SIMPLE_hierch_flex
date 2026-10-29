!@descr: serial export of native observation statistics in reference complex-variance units
module simple_flex_pca_observation_export
use simple_core_module_api, only: dp, string, fopen
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
private
public :: open_observation_export, write_observation_export, close_observation_export
contains

    !> The caller owns source/model/mean/basis/whitening hash bindings and row-role validation.
    !! Reference complex variance one defines units only; it is not an estimated noise level.
    subroutine open_observation_export( filename, model_id, role, rank, nrows, unit, status )
        character(len=*), intent(in) :: filename, model_id, role
        integer, intent(in) :: rank, nrows
        integer, intent(out) :: unit, status
        unit=0; status=1
        if( rank<1 .or. nrows<1 .or. len_trim(model_id)<1 ) return
        if( scan(model_id,achar(10)//achar(13))/=0 ) return
        if( role/='heldout' .and. role/='inference' ) return
        call fopen(unit,file=string(filename),status='NEW',action='WRITE',iostat=status)
        if( status/=0 ) return
        write(unit,'(A)',iostat=status) 'SCHEMA SIMPLE_FLEX_OBSERVATION_STATISTICS_V1'
        if( status==0 ) write(unit,'(A)',iostat=status) 'UNITS reference_complex_variance_1_real_coordinate_variance_half'
        if( status==0 ) write(unit,'(2A)',iostat=status) 'ROLE ',trim(role)
        if( status==0 ) write(unit,'(2A)',iostat=status) 'MODEL_ID ',trim(model_id)
        if( status==0 ) write(unit,'(A,2(1X,I0))',iostat=status) 'DIMENSIONS',rank,nrows
        if( status==0 ) write(unit,'(A)',iostat=status) 'CONVENTION fixed_pose_unit_amplitude_assumed_transformed_noise'
        if( status/=0 ) call finish_observation_export(unit,status)
    end subroutine open_observation_export

    !> Called serially from the shared native row former; this routine never fits or scores.
    !! G,b,s supplied at complex variance nu become Gref,bref,sref by multiplying by nu.
    subroutine write_observation_export( unit, pind, ndof, nu, gram, rhs, energy, status )
        integer, intent(in) :: unit, pind, ndof
        real(dp), intent(in) :: nu, gram(:,:), rhs(:), energy
        integer, intent(out) :: status
        real(dp) :: reference_gram(size(rhs),size(rhs)), reference_rhs(size(rhs)), reference_energy
        integer :: n, j
        status=1; n=size(rhs)
        if( n<1 .or. pind<1 .or. ndof<1 ) return
        if( any(shape(gram)/=[n,n]) ) return
        if( .not. ieee_is_finite(nu) .or. nu<=0._dp ) return
        if( .not. all(ieee_is_finite(gram)) .or. .not. all(ieee_is_finite(rhs)) ) return
        if( .not. ieee_is_finite(energy) .or. energy<0._dp ) return
        reference_gram=nu*gram; reference_rhs=nu*rhs; reference_energy=nu*energy
        if( .not. all(ieee_is_finite(reference_gram)) .or. .not. all(ieee_is_finite(reference_rhs)) ) return
        if( .not. ieee_is_finite(reference_energy) ) return
        write(unit,'(A,2(1X,I0),2(1X,ES26.17E3))',iostat=status) 'ROW',pind,ndof,nu,reference_energy
        if( status/=0 ) return
        write(unit,'(A,*(1X,ES26.17E3))',iostat=status) 'RHS',reference_rhs
        if( status/=0 ) return
        do j=1,n
            write(unit,'(A,1X,I0,*(1X,ES26.17E3))',iostat=status) 'GRAM_COLUMN',j,reference_gram(:,j)
            if( status/=0 ) return
        enddo
        write(unit,'(A,1X,I0)',iostat=status) 'ROW_END',pind
    end subroutine write_observation_export

    !> Every requested row gets a final validity record, including rows skipped by geometry checks.
    subroutine close_observation_export( unit, rows, ndof, valid, status )
        integer, intent(in) :: unit, rows(:), ndof(:)
        logical, intent(in) :: valid(:)
        integer, intent(out) :: status
        integer :: i
        status=1
        if( size(rows)==size(ndof) .and. size(rows)==size(valid) )then
            status=0
            do i=1,size(rows)
                write(unit,'(A,2(1X,I0),1X,L1)',iostat=status) 'ROW_STATUS',rows(i),ndof(i),valid(i)
                if( status/=0 ) exit
            enddo
            if( status==0 ) write(unit,'(A,2(1X,I0))',iostat=status) 'COMPLETE',size(rows),count(valid)
        endif
        call finish_observation_export(unit,status)
    end subroutine close_observation_export

    !> Attempt both operations, preserving the first write/validation/flush/close failure.
    !! Return status so the caller controls when to report failure.
    subroutine finish_observation_export( unit, status )
        integer, intent(in) :: unit
        integer, intent(inout) :: status
        integer :: io_status
        flush(unit,iostat=io_status)
        if( status==0 ) status=io_status
        close(unit,iostat=io_status)
        if( status==0 ) status=io_status
    end subroutine finish_observation_export
end module simple_flex_pca_observation_export
