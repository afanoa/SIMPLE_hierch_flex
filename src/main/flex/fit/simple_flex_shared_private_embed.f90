!@descr: frozen shared/private model scoring through the native observation evaluator
module simple_flex_shared_private_embed
use simple_core_module_api,     only: dp, simple_exception, string, fopen, fclose, fileiochk
use simple_parameters,          only: parameters
use simple_builder,             only: builder
use simple_flex_pca_records,     only: flex_selection
use simple_flex_pca_fit_types,   only: flex_fit
use simple_flex_pca_embed,       only: evaluate_predictive_rows, infer_fixed_rows
use simple_flex_shared_private, only: sp_posterior
use simple_flex_pca_observation_export, only: open_observation_export, write_observation_export, close_observation_export
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
implicit none
private
public :: score_shared_private_fit, evaluate_shared_private_fit, export_shared_private_inference
#include "simple_local_flags.inc"
contains

    !> The same immutable full prior used by fitting, with all image normalization terms.
    !! No union covariance inverse, covariance floor, merged prior or diagonal fallback.
    subroutine score_shared_private_fit( fit, gram, rhs, energy, logdet_noise, ndof, &
        &logp, valid, responsibility, posterior_mean, posterior_second, component_mean, component_second )
        class(flex_fit), intent(in) :: fit
        real(dp), intent(in) :: gram(:,:), rhs(:), energy, logdet_noise
        integer, intent(in) :: ndof
        real(dp), intent(out) :: logp, responsibility(:), posterior_mean(:), posterior_second(:,:)
        real(dp), optional, intent(out) :: component_mean(:,:), component_second(:,:,:)
        logical, intent(out) :: valid
        integer :: n, d, k, status
        real(dp) :: norm0
        valid=.false.; logp=ieee_value(0._dp,ieee_quiet_nan)
        responsibility=0._dp; posterior_mean=0._dp; posterior_second=0._dp
        if( present(component_mean) ) component_mean=0._dp
        if( present(component_second) ) component_second=0._dp
        if( .not. fit%spec%l_sp_req ) return
        if( .not. fit%history%l_mix_active .or. .not. fit%history%l_mix_used ) return
        if( .not. allocated(fit%history%sp_loading) ) return
        if( .not. allocated(fit%history%sp_sigma) ) return
        if( .not. allocated(fit%history%mix_xi) ) return
        if( .not. allocated(fit%history%mix_pi) ) return
        n=fit%model%ncomp; k=fit%spec%kmix
        d=fit%spec%sp_shared+fit%spec%sp_private
        if( n<1 .or. k<2 .or. d<1 .or. ndof<1 ) return
        if( any(shape(gram)/=[n,n]) .or. size(rhs)/=n ) return
        if( size(responsibility)/=k .or. size(posterior_mean)/=n ) return
        if( any(shape(posterior_second)/=[n,n]) ) return
        if( any(shape(fit%history%sp_loading)/=[n,d,k]) ) return
        if( any(shape(fit%history%sp_sigma)/=[d,d,k]) ) return
        if( any(shape(fit%history%mix_xi)/=[n,k]) .or. size(fit%history%mix_pi)/=k ) return
        if( present(component_mean) )then
            if( any(shape(component_mean)/=[n,k]) ) return
        endif
        if( present(component_second) )then
            if( any(shape(component_second)/=[n,n,k]) ) return
        endif
        norm0=real(ndof,dp)*log(2._dp*acos(-1._dp))+logdet_noise
        block
            real(dp) :: means(n,k), seconds(n,n,k)
            call sp_posterior(n,d,k,gram,rhs,energy,norm0,fit%history%sp_loading,fit%history%mix_xi, &
                &fit%history%sp_sigma,fit%history%mix_pi,logp,responsibility,posterior_mean,posterior_second, &
                &means,seconds,status)
            if( status==0 )then
                if( present(component_mean) ) component_mean=means
                if( present(component_second) ) component_second=seconds
            endif
        end block
        valid=status==0
        if( valid ) valid=ieee_is_finite(logp)
        if( .not. valid ) logp=ieee_value(0._dp,ieee_quiet_nan)
    end subroutine score_shared_private_fit

    !> New observations get a fresh E-step in the final fit frame through the shared row former.
    !! Responsibilities are diagnostic memberships; physical-state admission is a caller decision.
    subroutine evaluate_shared_private_fit( params, build, fit, train_sel, eval_pinds, fixed_sig2, &
        &logp, ndof, valid, responsibility, export_status, posterior_mean, posterior_second )
        class(parameters), intent(inout) :: params
        type(builder), intent(inout) :: build
        class(flex_fit), intent(in) :: fit
        type(flex_selection), intent(in) :: train_sel
        integer, intent(in) :: eval_pinds(:)
        real(dp), intent(in) :: fixed_sig2
        real(dp), allocatable, intent(out) :: logp(:), responsibility(:,:)
        integer, allocatable, intent(out) :: ndof(:)
        logical, allocatable, intent(out) :: valid(:)
        integer, intent(out) :: export_status
        real(dp), allocatable, optional, intent(out) :: posterior_mean(:,:), posterior_second(:,:,:)
        integer, allocatable :: row_slot(:)
        integer :: n, k, nrows, i, naudit, statistics_unit, open_status, write_status, close_status
        integer :: statistics_length, statistics_env_status
        real(dp), allocatable :: zero_gram(:,:), zero_rhs(:), check_resp(:), check_mean(:), check_second(:,:)
        real(dp) :: check_logp, missing
        real(dp), allocatable :: audit_mean(:,:), audit_second(:,:,:), audit_component_mean(:,:,:)
        real(dp), allocatable :: audit_component_second(:,:,:,:)
        character(len=:), allocatable :: audit_prefix, model_identity
        character(len=3) :: statistics_env
        logical :: check_valid, l_statistics
        export_status=0; write_status=0
        n=fit%model%ncomp; k=fit%spec%kmix
        if( n<1 .or. k<2 ) THROW_HARD('private evaluation requires a fitted mixture')
        allocate(zero_gram(n,n),zero_rhs(n),check_resp(k),check_mean(n),check_second(n,n))
        zero_gram=0._dp; zero_rhs=0._dp
        call score_shared_private_fit(fit,zero_gram,zero_rhs,0._dp,0._dp,1, &
            &check_logp,check_valid,check_resp,check_mean,check_second)
        if( .not. check_valid ) THROW_HARD('private evaluation requires the complete valid frozen prior')
        nrows=build%spproj_field%get_noris()
        if( size(eval_pinds)<1 .or. nrows<1 ) THROW_HARD('private evaluation requires nonempty source rows')
        if( any(eval_pinds<1) .or. any(eval_pinds>nrows) ) THROW_HARD('private evaluation row is outside source project')
        allocate(row_slot(nrows),source=0)
        do i=1,size(eval_pinds)
            if( row_slot(eval_pinds(i))/=0 ) THROW_HARD('duplicate private evaluation row')
            row_slot(eval_pinds(i))=i
        enddo
        missing=ieee_value(0._dp,ieee_quiet_nan)
        allocate(responsibility(k,size(eval_pinds)),source=missing)
        if( present(posterior_mean) ) allocate(posterior_mean(n,size(eval_pinds)),source=missing)
        if( present(posterior_second) ) allocate(posterior_second(n,n,size(eval_pinds)),source=missing)
        naudit=min(2,size(eval_pinds))
        allocate(audit_mean(n,naudit),audit_second(n,n,naudit),source=missing)
        allocate(audit_component_mean(n,k,naudit),audit_component_second(n,n,k,naudit),source=missing)
        audit_prefix=fit%spec%fprefix%to_char()//'_shared_private'
        model_identity=audit_prefix//'_model.txt'
        ! Match the native two-row audit convention: only the exact value yes enables export.
        call get_environment_variable('SIMPLE_PCA_EVAL_STATS',statistics_env, &
            &length=statistics_length,status=statistics_env_status)
        l_statistics=.false.
        if( statistics_env_status==0 ) l_statistics=statistics_length==3 .and. statistics_env=='yes'
        if( l_statistics )then
            call open_observation_export(audit_prefix//'_unit_statistics.txt',model_identity,'heldout', &
                &n,size(eval_pinds),statistics_unit,open_status)
            if( open_status/=0 ) THROW_HARD('cannot open private reference-statistics export')
            call evaluate_predictive_rows(params,build,fit%model,train_sel,eval_pinds,fixed_sig2, &
                &logp,ndof,valid,score_row=score_one,audit_prefix=audit_prefix,score_model_id=model_identity, &
                &statistics_row=observe_one,statistics_status=export_status)
            call close_observation_export(statistics_unit,eval_pinds,ndof,valid,close_status)
            if( export_status==0 ) export_status=write_status
            if( export_status==0 ) export_status=close_status
        else
            call evaluate_predictive_rows(params,build,fit%model,train_sel,eval_pinds,fixed_sig2, &
                &logp,ndof,valid,score_row=score_one,audit_prefix=audit_prefix,score_model_id=model_identity)
        endif
        call write_prediction_audit(fit,eval_pinds(:naudit),logp(:naudit),ndof(:naudit),valid(:naudit), &
            &responsibility(:,:naudit),audit_mean,audit_second,audit_component_mean,audit_component_second)
    contains
        subroutine observe_one( pind, half, gram, rhs, energy, nu, logdet_noise, row_ndof, status )
            integer, intent(in) :: pind, half, row_ndof
            real(dp), intent(in) :: gram(:,:), rhs(:), energy, nu, logdet_noise
            integer, intent(out) :: status
            real(dp) :: expected_logdet, tolerance
            status=write_status
            if( status/=0 ) return
            ! Observation export is independent of the numerical score and validity.
            if( half<0 .or. half>1 .or. row_ndof<2 .or. mod(row_ndof,2)/=0 )then
                write_status=1; status=write_status
                return
            endif
            if( .not. ieee_is_finite(nu) .or. nu<=0._dp .or. .not. ieee_is_finite(logdet_noise) )then
                write_status=1; status=write_status
                return
            endif
            expected_logdet=real(row_ndof,dp)*(log(nu)-log(2._dp))
            tolerance=64._dp*epsilon(1._dp)*max(1._dp,abs(expected_logdet))
            if( abs(logdet_noise-expected_logdet)>tolerance )then
                write_status=1; status=write_status
                return
            endif
            call write_observation_export(statistics_unit,pind,row_ndof,nu,gram,rhs,energy,status)
            if( status/=0 ) write_status=status
        end subroutine observe_one
        subroutine score_one( pind, gram, rhs, energy, logdet_noise, row_ndof, row_logp, row_valid )
            integer, intent(in) :: pind, row_ndof
            real(dp), intent(in) :: gram(:,:), rhs(:), energy, logdet_noise
            real(dp), intent(out) :: row_logp
            logical, intent(out) :: row_valid
            real(dp) :: row_resp(k), row_mean(n), row_second(n,n)
            real(dp) :: row_component_mean(n,k), row_component_second(n,n,k)
            integer :: slot
            row_valid=.false.; row_logp=missing
            if( pind<1 .or. pind>size(row_slot) ) return
            slot=row_slot(pind)
            if( slot<1 ) return
            call score_shared_private_fit(fit,gram,rhs,energy,logdet_noise,row_ndof, &
                &row_logp,row_valid,row_resp,row_mean,row_second,row_component_mean,row_component_second)
            if( .not. row_valid ) return
            responsibility(:,slot)=row_resp
            if( present(posterior_mean) ) posterior_mean(:,slot)=row_mean
            if( present(posterior_second) ) posterior_second(:,:,slot)=row_second
            if( slot<=naudit )then
                audit_mean(:,slot)=row_mean; audit_second(:,:,slot)=row_second
                audit_component_mean(:,:,slot)=row_component_mean
                audit_component_second(:,:,:,slot)=row_component_second
            endif
        end subroutine score_one
    end subroutine evaluate_shared_private_fit

    subroutine write_prediction_audit( fit, rows, logp, ndof, valid, resp, mean, second, component_mean, component_second )
        class(flex_fit), intent(in) :: fit
        integer, intent(in) :: rows(:), ndof(:)
        logical, intent(in) :: valid(:)
        real(dp), intent(in) :: logp(:), resp(:,:), mean(:,:), second(:,:,:)
        real(dp), intent(in) :: component_mean(:,:,:), component_second(:,:,:,:)
        character(len=:), allocatable :: prefix
        integer :: unit, ios, b, j, i
        prefix=fit%spec%fprefix%to_char()//'_shared_private'
        call fopen(unit,file=string(prefix//'_posterior_audit.txt'),status='NEW',action='WRITE',iostat=ios)
        call fileiochk('private posterior audit',ios)
        write(unit,'(A)') 'simple_flex_shared_private_posterior 1'
        write(unit,'(A)') 'score_mode conditional_full_prior_gaussian_mixture'
        write(unit,'(A)') 'convention fixed_pose_unit_amplitude_assumed_transformed_noise'
        write(unit,'(A,A)') 'model_receipt ',prefix//'_model.txt'
        write(unit,'(A,4(1X,I0))') 'dimensions',fit%model%ncomp,fit%spec%kmix, &
            &fit%spec%sp_shared+fit%spec%sp_private,size(rows)
        write(unit,'(A,*(1X,ES25.16))') 'pi',fit%history%mix_pi
        do b=1,fit%spec%kmix
            write(unit,'(A,1X,I0,*(1X,ES25.16))') 'xi',b,fit%history%mix_xi(:,b)
            do j=1,size(fit%history%sp_loading,2)
                write(unit,'(A,2(1X,I0),*(1X,ES25.16))') 'loading_column',b,j,fit%history%sp_loading(:,j,b)
                write(unit,'(A,2(1X,I0),*(1X,ES25.16))') 'sigma_column',b,j,fit%history%sp_sigma(:,j,b)
            enddo
        enddo
        do i=1,size(rows)
            write(unit,'(A,2(1X,I0),1X,L1,1X,ES25.16)') 'row',rows(i),ndof(i),valid(i),logp(i)
            if( .not. valid(i) ) cycle
            write(unit,'(A,*(1X,ES25.16))') 'responsibility',resp(:,i)
            write(unit,'(A,*(1X,ES25.16))') 'mean',mean(:,i)
            do j=1,fit%model%ncomp
                write(unit,'(A,1X,I0,*(1X,ES25.16))') 'second_column',j,second(:,j,i)
            enddo
            do b=1,fit%spec%kmix
                write(unit,'(A,1X,I0,*(1X,ES25.16))') 'component_mean',b,component_mean(:,b,i)
                do j=1,fit%model%ncomp
                    write(unit,'(A,2(1X,I0),*(1X,ES25.16))') 'component_second_column',b,j,component_second(:,j,b,i)
                enddo
            enddo
        enddo
        call fclose(unit)
    end subroutine write_prediction_audit

    !> Unit-reference statistics for every explicitly requested node row, including F. Never a held-out score.
    !! Final-frame priors are read only; the separately frozen C scalar is applied by the offline consumer.
    subroutine export_shared_private_inference( params, build, fit, train_sel, rows )
        use simple_defs, only: logfhandle
        class(parameters), intent(inout) :: params
        type(builder), intent(inout) :: build
        class(flex_fit), intent(in) :: fit
        type(flex_selection), intent(in) :: train_sel
        integer, intent(in) :: rows(:)
        real(dp), allocatable :: logp(:)
        integer, allocatable :: ndof(:)
        logical, allocatable :: valid(:)
        integer :: unit, status, export_status, close_status
        character(len=:), allocatable :: prefix, model_id
        prefix=fit%spec%fprefix%to_char()//'_shared_private_inference'
        model_id=fit%spec%fprefix%to_char()//'_shared_private_model.txt'
        call open_observation_export(prefix//'_unit_statistics.txt',model_id,'inference', &
            &fit%model%ncomp,size(rows),unit,status)
        if( status/=0 ) THROW_HARD('cannot open non-heldout inference statistics')
        call infer_fixed_rows(params,build,fit%model,train_sel,rows,1._dp,logp,ndof,valid, &
            &score_row=score_one_inference,audit_prefix=prefix,score_model_id=model_id, &
            &statistics_row=observe_inference,statistics_status=export_status)
        call close_observation_export(unit,rows,ndof,valid,close_status)
        if( export_status==0 ) export_status=close_status
        write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA NON-HELDOUT inference export valid=',count(valid),' / ',size(rows)
        call flush(logfhandle)
        if( export_status/=0 ) THROW_HARD('non-heldout inference statistics export failed')
        if( .not. all(valid) ) THROW_HARD('non-heldout inference requires every requested physical row valid')
    contains
        subroutine observe_inference( pind, half, gram, rhs, energy, nu, logdet_noise, dof, row_status )
            integer, intent(in) :: pind, half, dof
            real(dp), intent(in) :: gram(:,:), rhs(:), energy, nu, logdet_noise
            integer, intent(out) :: row_status
            row_status=1
            if( half/=0 .and. half/=1 ) return
            if( dof<2 .or. mod(dof,2)/=0 ) return
            if( nu/=1._dp .or. .not. ieee_is_finite(logdet_noise) ) return
            if( abs(logdet_noise+real(dof,dp)*log(2._dp))>64*epsilon(1._dp)*real(dof,dp) ) return
            call write_observation_export(unit,pind,dof,nu,gram,rhs,energy,row_status)
        end subroutine observe_inference
        subroutine score_one_inference( pind, gram, rhs, energy, logdet_noise, dof, value, ok )
            integer, intent(in) :: pind, dof
            real(dp), intent(in) :: gram(:,:), rhs(:), energy, logdet_noise
            real(dp), intent(out) :: value
            logical, intent(out) :: ok
            real(dp) :: resp(fit%spec%kmix), mean(fit%model%ncomp), second(fit%model%ncomp,fit%model%ncomp)
            call score_shared_private_fit(fit,gram,rhs,energy,logdet_noise,dof,value,ok,resp,mean,second)
        end subroutine score_one_inference
    end subroutine export_shared_private_inference
end module simple_flex_shared_private_embed
