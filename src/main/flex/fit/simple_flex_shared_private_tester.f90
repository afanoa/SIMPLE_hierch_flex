!@descr: shared/private full-covariance posterior, update and frame regression checks
module simple_flex_shared_private_tester
use simple_core_module_api, only: dp
use simple_flex_shared_private, only: sp_posterior, sp_update, sp_transport
use simple_flex_shared_private_embed, only: score_shared_private_fit
use simple_flex_pca_fit_types, only: flex_fit
use simple_flex_pca_posterior, only: probe_solve_mix
use simple_test_utils, only: assert_true, assert_int
use simple_flex_pca_observation_export, only: write_observation_export, close_observation_export
use omp_lib, only: omp_get_num_threads
implicit none
private
public :: test_flex_shared_private
contains
    subroutine test_flex_shared_private()
        integer, parameter :: N=3, D=2, K=2
        real(dp) :: g(N,N), h(N), loading(N,D,K), xi(N,K), sigma(D,D,K), pi(K)
        real(dp) :: resp(K), mean(N), second(N,N), cm(N,K), cs(N,N,K), logp
        real(dp) :: expected_resp(K), expected_mean(N), expected_second(N,N)
        real(dp) :: transform(N,N), inverse(N,N), moved(N,D,K), shifted(N,K), moved_g(N,N), moved_h(N)
        real(dp) :: old_loading(N,D,K), old_xi(N,K), old_sigma(D,D,K), count(K), first(N,K), moments(N,N,K)
        real(dp) :: expected_logp, energy, norm0
        integer :: status, b, i, j
        g = reshape([ &
            1.39999999999999991E+00_dp, 2.00000000000000011E-01_dp, -1.00000000000000006E-01_dp, &
            2.00000000000000011E-01_dp, 1.10000000000000009E+00_dp, 2.99999999999999989E-01_dp, &
            -1.00000000000000006E-01_dp, 2.99999999999999989E-01_dp, 9.00000000000000022E-01_dp ], shape(g))
        h = [ &
            2.99999999999999989E-01_dp, -8.00000000000000044E-01_dp, 5.00000000000000000E-01_dp ]
        loading = reshape([ &
            1.00000000000000000E+00_dp, 0.00000000000000000E+00_dp, 2.99999999999999989E-01_dp, &
            2.00000000000000011E-01_dp, 1.00000000000000000E+00_dp, -1.00000000000000006E-01_dp, &
            1.00000000000000006E-01_dp, 1.00000000000000000E+00_dp, 4.00000000000000022E-01_dp, &
            1.00000000000000000E+00_dp, 0.00000000000000000E+00_dp, 2.99999999999999989E-01_dp ], shape(loading))
        xi = reshape([ &
            1.00000000000000006E-01_dp, -2.00000000000000011E-01_dp, 6.99999999999999956E-01_dp, &
            -2.99999999999999989E-01_dp, 4.00000000000000022E-01_dp, -5.00000000000000000E-01_dp ], shape(xi))
        sigma = reshape([ &
            9.00000000000000022E-01_dp, 2.00000000000000011E-01_dp, 2.00000000000000011E-01_dp, &
            1.39999999999999991E+00_dp, 4.20000000000000018E+00_dp, -2.99999999999999989E-01_dp, &
            -2.99999999999999989E-01_dp, 1.49999999999999994E-01_dp ], shape(sigma))
        pi = [ &
            6.99999999999999956E-01_dp, 2.99999999999999989E-01_dp ]
        expected_resp = [ &
            8.97572302567503422E-01_dp, 1.02427697432496634E-01_dp ]
        expected_mean = [ &
            1.37742999669938804E-01_dp, -5.61917159541375777E-01_dp, 6.44453018425136248E-01_dp ]
        expected_second = reshape([ &
            3.94098589290687251E-01_dp, -7.40943433557153391E-02_dp, 2.42898999686641193E-01_dp, &
            -7.40943433557153253E-02_dp, 8.87394649313295814E-01_dp, -4.85325723448859936E-01_dp, &
            2.42898999686641193E-01_dp, -4.85325723448859936E-01_dp, 6.52622369889442888E-01_dp ], shape(expected_second))
        transform = reshape([ &
            1.19999999999999996E+00_dp, 0.00000000000000000E+00_dp, 2.99999999999999989E-01_dp, &
            1.00000000000000006E-01_dp, 8.00000000000000044E-01_dp, 0.00000000000000000E+00_dp, &
            2.00000000000000011E-01_dp, -1.00000000000000006E-01_dp, 1.10000000000000009E+00_dp ], shape(transform))
        inverse = reshape([ &
            8.75621890547263693E-01_dp, -2.98507462686567158E-02_dp, -2.38805970149253727E-01_dp, &
            -1.09452736318407962E-01_dp, 1.25373134328358193E+00_dp, 2.98507462686567158E-02_dp, &
            -1.69154228855721400E-01_dp, 1.19402985074626863E-01_dp, 9.55223880597014907E-01_dp ], shape(inverse))
        energy = 2.70000000000000018e+00_dp
        norm0 = 5.51363119922803602e+00_dp
        expected_logp = -4.71762756903871239e+00_dp
        call sp_posterior(N,D,K,g,h,energy,norm0,loading,xi,sigma,pi,logp,resp,mean,second,cm,cs,status)
        call assert_int(status,0,'shared/private posterior accepts full active covariances')
        call assert_true(abs(logp-expected_logp)<1.d-11,'shared/private normalized image density')
        call assert_true(maxval(abs(resp-expected_resp))<1.d-11,'unequal covariance determinants affect responsibilities')
        call assert_true(maxval(abs(mean-expected_mean))<1.d-11,'global posterior mean includes component offset')
        call assert_true(maxval(abs(second-expected_second))<1.d-11,'global second moment includes posterior covariance')
        call assert_true(abs(sum(resp)-1.d0)<1.d-13,'one image contributes unit total responsibility')
        call test_probe_caller(g,h,loading,xi,sigma,pi,expected_logp,energy,norm0,expected_mean,expected_second)
        moved=loading; shifted=xi
        call sp_transport(N,D,K,transform,moved,shifted,status)
        call assert_int(status,0,'invertible nonorthogonal frame transport accepted')
        moved_g=matmul(matmul(transpose(inverse),g),inverse)
        moved_h=matmul(transpose(inverse),h)
        call sp_posterior(N,D,K,moved_g,moved_h,energy,norm0,moved,shifted,sigma,pi,logp,resp,mean,second,cm,cs,status)
        call assert_int(status,0,'posterior in transported frame valid')
        call assert_true(abs(logp-expected_logp)<1.d-11,'full nonorthogonal frame density invariant')
        call assert_true(maxval(abs(mean-matmul(transform,expected_mean)))<1.d-11,'posterior mean transports')
        call assert_true(maxval(abs(second-matmul(matmul(transform,expected_second),transpose(transform))))<1.d-10, &
            &'full posterior second moment transports')
        old_loading=loading; old_xi=xi; old_sigma=sigma
        count=[3.d0,2.d0]
        do b=1,K
            first(:,b)=count(b)*xi(:,b)
            moments(:,:,b)=matmul(matmul(loading(:,:,b),sigma(:,:,b)),transpose(loading(:,:,b)))
            do j=1,N
                do i=1,N
                    moments(i,j,b)=count(b)*(moments(i,j,b)+xi(i,b)*xi(j,b))
                enddo
            enddo
        enddo
        call sp_update(N,D,K,count,first,moments,loading,xi,sigma,pi,status)
        call assert_int(status,0,'full-moment shared/private update')
        call assert_true(maxval(abs(xi-old_xi))<1.d-12,'updated component means')
        call assert_true(maxval(abs(sigma-old_sigma))<1.d-11,'updated full active covariance')
        call assert_true(maxval(abs(pi-[.6d0,.4d0]))<1.d-13,'component masses normalized once')
        loading(:,2,1)=loading(:,1,1)
        call sp_posterior(N,D,K,g,h,energy,norm0,loading,old_xi,old_sigma,pi,logp,resp,mean,second,cm,cs,status)
        call assert_int(status,0,'deficient loading still defines a finite image density')
        call sp_update(N,D,K,count,first,moments,loading,xi,sigma,pi,status)
        call assert_true(status/=0,'inverse-loading update rejects deficient representation')
        loading=old_loading; sigma=old_sigma; sigma(1,1,1)=-1.d0
        call sp_posterior(N,D,K,g,h,energy,norm0,loading,old_xi,sigma,pi,logp,resp,mean,second,cm,cs,status)
        call assert_true(status/=0,'nonpositive active covariance rejected without hidden ridge')
    end subroutine test_flex_shared_private

    subroutine test_probe_caller( g, h, loading, xi, sigma, pi, logp, energy, norm0, mean, second )
        integer, parameter :: N=3, D=2, K=2, NT=4
        real(dp), intent(in) :: g(N,N), h(N), loading(N,D,K), xi(N,K), sigma(D,D,K), pi(K)
        real(dp), intent(in) :: logp, energy, norm0, mean(N), second(N,N)
        class(flex_fit), allocatable :: fit
        real(dp) :: ld(NT), nll(NT), expected_nll, serial_mean(N), serial_second(N,N)
        real(dp) :: serial_count(K), serial_first(N,K), serial_moments(N,N,K)
        real(dp) :: embedded_logp, embedded_resp(K), embedded_mean(N), embedded_second(N,N)
        real(dp), parameter :: AMPLITUDE=1.7d0, NOISE=0.6d0
        logical :: valid(NT), embedded_valid
        integer :: t, team(NT), update_status
        allocate(fit)
        fit%model%ncomp=N; fit%model%sig2=NOISE
        fit%spec%l_sp_req=.true.; fit%spec%sp_shared=1; fit%spec%sp_private=1; fit%spec%kmix=K
        fit%history%l_mix_active=.true.; fit%history%l_mix_used=.true.
        fit%history%sp_loading=loading; fit%history%sp_sigma=sigma
        fit%history%mix_xi=xi; fit%history%mix_pi=pi
        allocate(fit%iter%Gth(N,N,NT),fit%iter%bth(N,NT),fit%iter%cth(N,NT),fit%iter%zth(N,NT))
        allocate(fit%iter%Ainvth(N,N,NT),fit%iter%dens(N,N,NT))
        allocate(fit%history%mxa_sr(K,NT),source=0.d0)
        allocate(fit%history%mxa_sm(N,K,NT),source=0.d0)
        allocate(fit%history%mxa_smm(N,N,K,NT),source=0.d0)
        do t=1,NT
            fit%iter%Gth(:,:,t)=NOISE*g/AMPLITUDE**2
            fit%iter%cth(:,t)=[.2d0,-.3d0,.1d0]
            fit%iter%bth(:,t)=NOISE*h/AMPLITUDE+AMPLITUDE*fit%iter%cth(:,t)
        enddo
        call probe_solve_mix(fit,1,1,AMPLITUDE,ld(1),valid(1),nll(1))
        expected_nll=-2.d0*logp-energy-norm0
        call assert_true(valid(1),'actual shared/private probe caller succeeds')
        call assert_true(abs(nll(1)-expected_nll)<1.d-11,'probe adds component determinants exactly once')
        call assert_true(maxval(abs(fit%iter%zth(:,1)-mean))<1.d-11,'probe applies contrast and noise once')
        call assert_true(maxval(abs(fit%iter%dens(:,:,1)-second))<1.d-11,'probe retains full second moment')
        call test_reference_export(fit,g,h,energy,logp)
        serial_mean=fit%iter%zth(:,1); serial_second=fit%iter%dens(:,:,1)
        serial_count=fit%history%mxa_sr(:,1); serial_first=fit%history%mxa_sm(:,:,1)
        serial_moments=fit%history%mxa_smm(:,:,:,1)
        fit%history%mxa_sr=0.d0; fit%history%mxa_sm=0.d0; fit%history%mxa_smm=0.d0
        !$omp parallel do num_threads(NT) default(none) &
        !$omp shared(fit,ld,valid,nll,team) private(t)
        do t=1,NT
            team(t)=omp_get_num_threads()
            call probe_solve_mix(fit,t,t,AMPLITUDE,ld(t),valid(t),nll(t))
        enddo
        !$omp end parallel do
        call assert_true(all(team>1),'probe regression actually uses a parallel team')
        call assert_true(all(valid),'parallel probe calls succeed')
        do t=1,NT
            call assert_true(abs(nll(t)-expected_nll)<1.d-11,'parallel likelihood matches serial')
            call assert_true(maxval(abs(fit%iter%zth(:,t)-serial_mean))<1.d-12,'parallel posterior mean matches')
            call assert_true(maxval(abs(fit%iter%dens(:,:,t)-serial_second))<1.d-12,'parallel full moment matches')
            call assert_true(maxval(abs(fit%history%mxa_sr(:,t)-serial_count))<1.d-12,'thread counts isolated')
            call assert_true(maxval(abs(fit%history%mxa_sm(:,:,t)-serial_first))<1.d-12,'thread first moments isolated')
            call assert_true(maxval(abs(fit%history%mxa_smm(:,:,:,t)-serial_moments))<1.d-12, &
                &'thread full component second moments isolated')
        enddo
        call assert_true(abs(sum(fit%history%mxa_sr)-real(NT,dp))<1.d-12,'one count per image across threads')
        call score_shared_private_fit(fit,g,h,energy,norm0-real(N,dp)*log(2.d0*acos(-1.d0)),N, &
            &embedded_logp,embedded_valid,embedded_resp,embedded_mean,embedded_second)
        call assert_true(embedded_valid,'full fitted prior accepted by embedding consumer')
        call assert_true(abs(embedded_logp-logp)<1.d-11,'embedding consumer includes image normalization exactly once')
        call assert_true(maxval(abs(embedded_resp-serial_count))<1.d-12,'fit and embedding responsibilities agree')
        call assert_true(maxval(abs(embedded_mean-serial_mean))<1.d-12,'fit and embedding means agree')
        call assert_true(maxval(abs(embedded_second-serial_second))<1.d-12,'fit and embedding full moments agree')
        call sp_update(N,D,K,serial_count,serial_first,serial_moments,fit%history%sp_loading, &
            &fit%history%mix_xi,fit%history%sp_sigma,fit%history%mix_pi,update_status)
        call assert_int(update_status,0,'component model update retains full covariance')
        fit%history%mxa_sr=0.d0; fit%history%mxa_sm=0.d0; fit%history%mxa_smm=0.d0
        call probe_solve_mix(fit,1,1,AMPLITUDE,ld(1),valid(1),nll(1))
        call score_shared_private_fit(fit,g,h,energy,norm0-real(N,dp)*log(2.d0*acos(-1.d0)),N, &
            &embedded_logp,embedded_valid,embedded_resp,embedded_mean,embedded_second)
        call assert_true(valid(1).and.embedded_valid,'updated fitted model supports a fresh embedding E-step')
        call assert_true(abs(nll(1)+energy+norm0+2.d0*embedded_logp)<1.d-11,'updated model score equivalence')
        call assert_true(maxval(abs(embedded_mean-fit%iter%zth(:,1)))<1.d-12,'updated model posterior equivalence')
        call assert_true(maxval(abs(embedded_second-fit%iter%dens(:,:,1)))<1.d-12,'updated model full moment equivalence')
        serial_count=fit%history%mxa_sr(:,1)
        fit%history%sp_sigma(1,1,1)=-1.d0
        call probe_solve_mix(fit,1,1,AMPLITUDE,ld(1),valid(1),nll(1))
        call assert_true(.not.valid(1),'invalid private model reports failure without throwing inside parallel solve')
        call assert_true(maxval(abs(fit%history%mxa_sr(:,1)-serial_count))<1.d-12, &
            &'rejected model contributes no responsibility')
        call score_shared_private_fit(fit,g,h,energy,0.d0,N, &
            &embedded_logp,embedded_valid,embedded_resp,embedded_mean,embedded_second)
        call assert_true(.not.embedded_valid,'invalid full prior has no diagonal embedding fallback')
        call fit%kill
        deallocate(fit)
    end subroutine test_probe_caller

    subroutine test_reference_export( fit, g, h, energy, expected_logp )
        class(flex_fit), intent(in) :: fit
        real(dp), intent(in) :: g(3,3), h(3), energy, expected_logp
        real(dp), parameter :: SUPPLIED_NU=.37d0, TARGET_NU=2.d0
        real(dp) :: reference_g(3,3), reference_h(3), reference_energy, supplied_read
        real(dp) :: score, resp(2), mean(3), second(3,3)
        integer :: unit, status, row, dof, column, j
        logical :: valid
        character(len=32) :: tag
        open(newunit=unit,status='scratch',form='formatted',action='readwrite')
        call write_observation_export(unit,71,3,SUPPLIED_NU, &
            &TARGET_NU*g/SUPPLIED_NU,TARGET_NU*h/SUPPLIED_NU,TARGET_NU*energy/SUPPLIED_NU,status)
        call assert_int(status,0,'reference statistics export accepts nonunit supplied variance')
        rewind(unit)
        read(unit,*) tag,row,dof,supplied_read,reference_energy
        call assert_true(trim(tag)=='ROW','reference export row record')
        call assert_int(row,71,'reference export preserves original project row')
        call assert_int(dof,3,'reference export preserves real degrees of freedom')
        call assert_true(abs(supplied_read-SUPPLIED_NU)<1.d-15,'reference export declares supplied variance')
        read(unit,*) tag,reference_h
        call assert_true(trim(tag)=='RHS','reference export right hand side')
        do j=1,3
            read(unit,*) tag,column,reference_g(:,j)
            call assert_true(trim(tag)=='GRAM_COLUMN','reference export matrix orientation')
            call assert_int(column,j,'reference export column index')
        enddo
        read(unit,*) tag,row
        call assert_true(trim(tag)=='ROW_END','reference export complete row')
        call assert_true(maxval(abs(reference_g-TARGET_NU*g))<1.d-12,'Gram converted to reference units once')
        call assert_true(maxval(abs(reference_h-TARGET_NU*h))<1.d-12,'right hand side reference conversion')
        call assert_true(abs(reference_energy-TARGET_NU*energy)<1.d-12,'energy reference conversion')
        call score_shared_private_fit(fit,reference_g/TARGET_NU,reference_h/TARGET_NU, &
            &reference_energy/TARGET_NU,real(dof,dp)*log(TARGET_NU/2.d0),dof,score,valid,resp,mean,second)
        call assert_true(valid,'exported statistics support a calibrated full-prior posterior')
        call assert_true(abs(score-expected_logp)<1.d-11,'reference rescaling preserves normalized dense golden density')
        call write_observation_export(unit,71,3,0.d0,g,h,energy,status)
        call assert_true(status/=0,'reference export rejects zero variance without a floor')
        call close_observation_export(unit,[71],[3],[.true.],status)
        call assert_int(status,0,'reference export final row-status close')
    end subroutine test_reference_export
end module simple_flex_shared_private_tester
