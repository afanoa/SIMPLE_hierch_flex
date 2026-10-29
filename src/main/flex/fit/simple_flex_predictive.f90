!@descr: normalized Gaussian latent predictive scores from fixed real observation statistics
module simple_flex_predictive
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
use simple_defs,   only: dp, DPI
use simple_linalg, only: cholesky, chol_forward
implicit none
private
public :: flex_predictive_score, flex_predictive_mixture

contains

    ! G=D' N^-1 D, b=D' N^-1(y-m), s=(y-m)' N^-1(y-m); latent covariance is L L'.
    ! ndof counts independent REAL coordinates and logdet_noise is log|N| in that same measure.
    ! Failure returns NaN/ok=false; no ridge, covariance repair or cancellation clipping is applied.
    pure subroutine flex_predictive_score( gram, b, s, logdet_noise, ndof, latent_chol, logp, ok )
        real(dp), intent(in)  :: gram(:,:), b(:), s, logdet_noise, latent_chol(:,:)
        integer,  intent(in)  :: ndof
        real(dp), intent(out) :: logp
        logical,  intent(out) :: ok
        real(dp) :: h(size(b),size(b)), factor(size(b),size(b)), v(size(b)), w(size(b))
        real(dp) :: logdet_h, removed, quadratic, bound, value, scale
        integer  :: q, i
        logical  :: factored
        logp = ieee_value(0._dp, ieee_quiet_nan)
        ok = .false.
        q = size(b)
        if( any(shape(gram) /= [q,q]) ) return
        if( any(shape(latent_chol) /= [q,q]) ) return
        if( ndof < 1 ) return
        if( .not. ieee_is_finite(s) .or. .not. ieee_is_finite(logdet_noise) ) return
        if( s < 0._dp ) return
        if( .not. all(ieee_is_finite(gram)) .or. .not. all(ieee_is_finite(b)) ) return
        if( .not. all(ieee_is_finite(latent_chol)) ) return
        logdet_h = 0._dp
        removed = 0._dp
        if( q > 0 )then
            scale = maxval(abs(gram))
            bound = 64._dp * epsilon(1._dp) * real(q,dp) * scale
            if( maxval(abs(gram - transpose(gram))) > bound ) return
            do i = 1, q
                if( latent_chol(i,i) <= 0._dp ) return
                if( any(latent_chol(i,i+1:q) /= 0._dp) ) return
            enddo
            if( .not. valid_joint_gram(gram, b, s) ) return
            h = matmul(transpose(latent_chol), matmul(gram, latent_chol))
            do i = 1, q
                h(i,i) = h(i,i) + 1._dp
            enddo
            if( .not. all(ieee_is_finite(h)) ) return
            call cholesky(h, factor, q, factored)
            if( .not. factored ) return
            if( .not. all(ieee_is_finite(factor)) ) return
            v = matmul(transpose(latent_chol), b)
            if( .not. all(ieee_is_finite(v)) ) return
            call chol_forward(factor, v, w, q)
            if( .not. all(ieee_is_finite(w)) ) return
            removed = dot_product(w,w)
            do i = 1, q
                logdet_h = logdet_h + 2._dp * log(factor(i,i))
            enddo
        endif
        if( .not. ieee_is_finite(removed) ) return
        quadratic = s - removed
        bound = 64._dp * epsilon(1._dp) * (real(ndof,dp) + 2._dp*real(q,dp) + 1._dp) * max(s,removed)
        if( quadratic < 0._dp ) return
        ! Near-total subtraction needs observation-space scoring instead of a fabricated residual.
        if( removed > 0._dp )then
            if( quadratic <= bound ) return
        endif
        value = -.5_dp * (real(ndof,dp)*log(2._dp*DPI) + logdet_noise + logdet_h + quadratic)
        if( .not. ieee_is_finite(value) ) return
        logp = value
        ok = .true.
    end subroutine flex_predictive_score

    ! Component densities already include their own Gaussian constants on the same observation support.
    ! Nonnegative weights must sum to one within roundoff; inactive component scores are ignored.
    pure subroutine flex_predictive_mixture( component_logp, weights, logp, responsibilities, ok )
        real(dp), intent(in)  :: component_logp(:), weights(:)
        real(dp), intent(out) :: logp, responsibilities(:)
        logical,  intent(out) :: ok
        real(dp) :: joint(size(weights)), largest, total, weight_sum, value
        integer  :: k, j
        logp = ieee_value(0._dp, ieee_quiet_nan)
        responsibilities = 0._dp
        ok = .false.
        k = size(weights)
        if( k < 1 ) return
        if( size(component_logp) /= k .or. size(responsibilities) /= k ) return
        if( .not. all(ieee_is_finite(weights)) ) return
        if( any(weights < 0._dp) ) return
        weight_sum = sum(weights)
        if( .not. ieee_is_finite(weight_sum) ) return
        if( abs(weight_sum - 1._dp) > 64._dp*epsilon(1._dp)*real(k,dp) ) return
        if( weight_sum <= 0._dp ) return
        largest = -huge(1._dp)
        joint = 0._dp
        do j = 1, k
            if( weights(j) <= 0._dp ) cycle
            if( .not. ieee_is_finite(component_logp(j)) ) return
            joint(j) = component_logp(j) + log(weights(j)) - log(weight_sum)
            if( .not. ieee_is_finite(joint(j)) ) return
            largest = max(largest, joint(j))
        enddo
        total = 0._dp
        do j = 1, k
            if( weights(j) <= 0._dp ) cycle
            joint(j) = exp(joint(j) - largest)
            total = total + joint(j)
        enddo
        value = largest + log(total)
        if( .not. ieee_is_finite(value) ) return
        logp = value
        responsibilities = joint / total
        ok = .true.
    end subroutine flex_predictive_mixture

    ! Pivoted semidefinite elimination validates the joint Gram without changing the scored statistics.
    pure logical function valid_joint_gram( gram, b, s ) result( valid )
        real(dp), intent(in) :: gram(:,:), b(:), s
        real(dp) :: a(size(b)+1,size(b)+1), row(size(b)+1), scale, tolerance, pivot
        integer  :: n, q, i, j, k, p
        valid = .false.
        q = size(b)
        n = q + 1
        a(1:q,1:q) = .5_dp*gram + .5_dp*transpose(gram)
        a(1:q,n) = b
        a(n,1:q) = b
        a(n,n) = s
        scale = maxval(abs(a))
        if( scale == 0._dp )then
            valid = .true.
            return
        endif
        a = a / scale
        tolerance = 64._dp * epsilon(1._dp) * real(n,dp)
        do i = 1, n
            p = i
            do j = i+1, n
                if( a(j,j) > a(p,p) ) p = j
            enddo
            row = a(i,:)
            a(i,:) = a(p,:)
            a(p,:) = row
            row = a(:,i)
            a(:,i) = a(:,p)
            a(:,p) = row
            pivot = a(i,i)
            if( pivot < -tolerance ) return
            if( pivot <= tolerance )then
                valid = maxval(abs(a(i:n,i:n))) <= tolerance
                return
            endif
            do j = i+1, n
                do k = j, n
                    a(k,j) = a(k,j) - a(k,i)*a(j,i)/pivot
                    a(j,k) = a(k,j)
                enddo
            enddo
        enddo
        valid = .true.
    end function valid_joint_gram

end module simple_flex_predictive
