!@descr: shared/private mixture inference from one image's noise-scaled sufficient statistics
module simple_flex_shared_private
use simple_core_module_api, only: dp
use simple_linalg,          only: cholesky, chol_forward, chol_backward, jacobi
use simple_kmeans,          only: kmeans
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
private
public :: sp_posterior, sp_update, sp_transport, sp_initialize
contains

    !> Full active covariance may imply a singular union-frame covariance; only active space is inverted.
    !! e is whitened residual energy; norm0 is image_dim*log(2*pi)+logdet(noise).
    !! Passing both zero returns the varying term for a caller which adds them once.
    subroutine sp_posterior( n, d, k, g, h, e, norm0, loading, xi, sigma, pi, &
        &logdensity, resp, zmean, zsecond, component_mean, component_second, status )
        integer, intent(in) :: n, d, k
        real(dp), intent(in) :: g(n,n), h(n), e, norm0
        real(dp), intent(in) :: loading(n,d,k), xi(n,k), sigma(d,d,k), pi(k)
        real(dp), intent(out) :: logdensity, resp(k), zmean(n), zsecond(n,n)
        real(dp), intent(out) :: component_mean(n,k), component_second(n,n,k)
        integer, intent(out) :: status
        real(dp) :: sinv(d,d), precision(d,d), cov(d,d), rhs(d), amean(d)
        real(dp) :: gc(n,d), active_cov(n,n), lw(k), lds, lda, peak, total
        logical :: ok
        integer :: b, i, j
        status = 1
        logdensity = 0.d0; resp = 0.d0; zmean = 0.d0; zsecond = 0.d0
        component_mean = 0.d0; component_second = 0.d0
        if( n < 1 .or. d < 1 .or. k < 1 .or. d > n ) return
        if( .not. all(ieee_is_finite(g)) .or. .not. all(ieee_is_finite(h)) ) return
        if( .not. all(ieee_is_finite(loading)) .or. .not. all(ieee_is_finite(xi)) ) return
        if( .not. all(ieee_is_finite(sigma)) .or. .not. all(ieee_is_finite(pi)) ) return
        if( .not. ieee_is_finite(e) .or. .not. ieee_is_finite(norm0) ) return
        if( any(pi <= 0.d0) .or. abs(sum(pi)-1.d0) > 1.d-8 ) return
        do b = 1, k
            call inverse_spd(sigma(:,:,b), d, sinv, lds, ok)
            if( .not. ok )then
                status = 2
                return
            endif
            gc = matmul(g, loading(:,:,b))
            precision = sinv + matmul(transpose(loading(:,:,b)), gc)
            call inverse_spd(precision, d, cov, lda, ok)
            if( .not. ok )then
                status = 3
                return
            endif
            rhs = matmul(transpose(loading(:,:,b)), h-matmul(g,xi(:,b)))
            amean = matmul(cov, rhs)
            lw(b) = log(pi(b)) - 0.5d0*(e - 2.d0*dot_product(h,xi(:,b)) + &
                &dot_product(xi(:,b),matmul(g,xi(:,b))) + lds + lda - dot_product(rhs,amean))
            component_mean(:,b) = xi(:,b) + matmul(loading(:,:,b),amean)
            active_cov = matmul(matmul(loading(:,:,b),cov),transpose(loading(:,:,b)))
            do j = 1, n
                do i = 1, n
                    component_second(i,j,b) = active_cov(i,j) + component_mean(i,b)*component_mean(j,b)
                enddo
            enddo
        enddo
        if( .not. all(ieee_is_finite(lw)) ) return
        peak = maxval(lw)
        resp = exp(lw-peak)
        total = sum(resp)
        if( total <= 0.d0 ) return
        resp = resp/total
        logdensity = peak + log(total) - 0.5d0*norm0
        zmean = matmul(component_mean,resp)
        do b = 1, k
            zsecond = zsecond + resp(b)*component_second(:,:,b)
        enddo
        status = 0
    end subroutine sp_posterior

    !> Reduced component second moments already include conditional covariance. Unique images
    !! contribute once, with responsibilities summing to one. No pooled-ancestor pseudo-counts.
    subroutine sp_update( n, d, k, count, first, second, loading, xi, sigma, pi, status )
        integer, intent(in) :: n, d, k
        real(dp), intent(in) :: count(k), first(n,k), second(n,n,k), loading(n,d,k)
        real(dp), intent(out) :: xi(n,k), sigma(d,d,k), pi(k)
        integer, intent(out) :: status
        real(dp) :: gram(d,d), inverse(d,d), left(d,n), centered(n,n), scratch(d,d), ld
        logical :: ok
        integer :: b, i, j
        status = 1; xi = 0.d0; sigma = 0.d0; pi = 0.d0
        if( n < 1 .or. d < 1 .or. d > n .or. k < 1 ) return
        if( .not. all(ieee_is_finite(count)) .or. any(count <= 0.d0) ) return
        if( .not. all(ieee_is_finite(first)) .or. .not. all(ieee_is_finite(second)) ) return
        if( .not. all(ieee_is_finite(loading)) ) return
        pi = count/sum(count)
        do b = 1, k
            xi(:,b) = first(:,b)/count(b)
            gram = matmul(transpose(loading(:,:,b)),loading(:,:,b))
            call inverse_spd(gram, d, inverse, ld, ok)
            if( .not. ok )then
                status = 2
                return
            endif
            left = matmul(inverse,transpose(loading(:,:,b)))
            centered = second(:,:,b)/count(b)
            do j = 1, n
                do i = 1, n
                    centered(i,j) = centered(i,j) - xi(i,b)*xi(j,b)
                enddo
            enddo
            sigma(:,:,b) = matmul(matmul(left,centered),transpose(left))
            sigma(:,:,b) = 0.5d0*(sigma(:,:,b)+transpose(sigma(:,:,b)))
            call inverse_spd(sigma(:,:,b), d, scratch, ld, ok)
            if( .not. ok )then
                status = 3
                return
            endif
        enddo
        status = 0
    end subroutine sp_update

    !> Exact same-rank coordinate change ONLY: caller establishes Wraw=Wnew*transform.
    !! Active Sigma is unchanged; union means/loadings both transform. A physical subspace
    !! change is an M-step followed by a new E-step, never reuse of old image posteriors.
    subroutine sp_transport( n, d, k, transform, loading, xi, status )
        integer, intent(in) :: n, d, k
        real(dp), intent(in) :: transform(n,n)
        real(dp), intent(inout) :: loading(n,d,k), xi(n,k)
        integer, intent(out) :: status
        real(dp) :: check(n,n), inverse(n,n), ld, moved(n,d,k), shifted(n,k)
        logical :: ok
        integer :: b
        status = 1
        if( .not. all(ieee_is_finite(transform)) ) return
        check = matmul(transpose(transform),transform)
        call inverse_spd(check, n, inverse, ld, ok)
        if( .not. ok ) return
        do b = 1, k
            moved(:,:,b) = matmul(transform,loading(:,:,b))
        enddo
        shifted = matmul(transform,xi)
        if( .not. all(ieee_is_finite(moved)) .or. .not. all(ieee_is_finite(shifted)) ) return
        loading = moved; xi = shifted
        status = 0
    end subroutine sp_transport

    !> Warm-start from point latents only; subsequent EM always uses full posterior moments.
    !! q shared directions come from pooled within-component scatter; r private directions
    !! come from each residual scatter. These are physical union-frame loadings, not PC IDs.
    subroutine sp_initialize( nobs, n, k, q, r, z, loading, xi, sigma, pi, status )
        integer, intent(in) :: nobs, n, k, q, r
        real(dp), intent(in) :: z(nobs,n)
        real(dp), intent(out) :: loading(n,q+r,k), xi(n,k), sigma(q+r,q+r,k), pi(k)
        integer, intent(out) :: status
        type(kmeans) :: km
        real(dp), allocatable :: rows(:,:), scatter(:,:,:), pooled(:,:), shared(:,:), residual(:,:), private(:,:)
        integer, allocatable :: labels(:), counts(:)
        real(dp) :: delta(n), projector(n,n), inverse(q+r,q+r), ld
        logical :: ok
        integer :: i, j, b, nr
        status = 1; loading = 0.d0; xi = 0.d0; sigma = 0.d0; pi = 0.d0
        if( k < 2 .or. q < 0 .or. r < 1 .or. q+r > n .or. nobs < k ) return
        if( .not. all(ieee_is_finite(z)) ) return
        nr = 0
        do i = 1, nobs
            if( any(z(i,:) /= 0.d0) ) nr = nr + 1
        enddo
        if( nr < k ) return
        allocate(rows(nr,n), labels(nr), counts(k), scatter(n,n,k), pooled(n,n))
        allocate(shared(n,q), residual(n,n), private(n,r))
        nr = 0
        do i = 1, nobs
            if( any(z(i,:) /= 0.d0) )then
                nr = nr + 1
                rows(nr,:) = z(i,:)
            endif
        enddo
        call km%new(rows,k)
        call km%cluster(labels,xi)
        call km%kill
        counts = 0; scatter = 0.d0
        do i = 1, nr
            b = labels(i)
            if( b < 1 .or. b > k ) return
            counts(b) = counts(b)+1
            delta = rows(i,:)-xi(:,b)
            do j = 1, n
                scatter(:,j,b) = scatter(:,j,b)+delta*delta(j)
            enddo
        enddo
        if( any(counts <= q+r) )then
            status = 2
            return
        endif
        pooled = sum(scatter,dim=3)/real(nr,dp)
        if( q > 0 )then
            call leading_vectors(pooled,n,q,shared,ok)
            if( .not. ok ) return
        endif
        projector = 0.d0
        do i = 1, n
            projector(i,i) = 1.d0
        enddo
        if( q > 0 ) projector = projector-matmul(shared,transpose(shared))
        do b = 1, k
            scatter(:,:,b) = scatter(:,:,b)/real(counts(b),dp)
            residual = matmul(matmul(projector,scatter(:,:,b)),projector)
            call leading_vectors(residual,n,r,private,ok)
            if( .not. ok )then
                status = 3
                return
            endif
            if( q > 0 ) loading(:,1:q,b) = shared
            loading(:,q+1:q+r,b) = private
            sigma(:,:,b) = matmul(matmul(transpose(loading(:,:,b)),scatter(:,:,b)),loading(:,:,b))
            call inverse_spd(sigma(:,:,b),q+r,inverse,ld,ok)
            if( .not. ok )then
                status = 4
                return
            endif
            pi(b) = real(counts(b),dp)/real(nr,dp)
        enddo
        status = 0
    end subroutine sp_initialize

    subroutine inverse_spd( matrix, n, inverse, logdet, ok )
        integer, intent(in) :: n
        real(dp), intent(in) :: matrix(n,n)
        real(dp), intent(out) :: inverse(n,n), logdet
        logical, intent(out) :: ok
        real(dp) :: chol(n,n), rhs(n), tmp(n), sol(n)
        integer :: j
        inverse = 0.d0; logdet = 0.d0; ok = .false.
        if( .not. all(ieee_is_finite(matrix)) ) return
        if( maxval(abs(matrix-transpose(matrix))) > 1.d-9*max(1.d0,maxval(abs(matrix))) ) return
        call cholesky(matrix,chol,n,ok)
        if( .not. ok ) return
        do j = 1, n
            if( chol(j,j) <= 0.d0 )then
                ok = .false.
                return
            endif
            logdet = logdet + 2.d0*log(chol(j,j))
            rhs = 0.d0; rhs(j) = 1.d0
            call chol_forward(chol,rhs,tmp,n)
            call chol_backward(chol,tmp,sol,n)
            inverse(:,j) = sol
        enddo
        inverse = 0.5d0*(inverse+transpose(inverse))
        ok = all(ieee_is_finite(inverse)) .and. ieee_is_finite(logdet)
    end subroutine inverse_spd

    subroutine leading_vectors( matrix, n, rank, vectors, ok )
        integer, intent(in) :: n, rank
        real(dp), intent(in) :: matrix(n,n)
        real(dp), intent(out) :: vectors(n,rank)
        logical, intent(out) :: ok
        real(dp) :: work(n,n), eigen(n), evector(n,n), maximum
        integer :: j, index, nrot
        ok = .false.; vectors = 0.d0
        work = 0.5d0*(matrix+transpose(matrix))
        call jacobi(work,n,n,eigen,evector,nrot)
        maximum = maxval(eigen)
        if( maximum <= 0.d0 ) return
        do j = 1, rank
            index = maxloc(eigen,dim=1)
            if( eigen(index) <= 1.d-12*maximum ) return
            vectors(:,j) = evector(:,index)
            eigen(index) = -huge(1.d0)
        enddo
        ok = .true.
    end subroutine leading_vectors
end module simple_flex_shared_private
