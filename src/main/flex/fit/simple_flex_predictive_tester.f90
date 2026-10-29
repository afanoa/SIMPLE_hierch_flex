!@descr: independent dense Gaussian and failure checks for normalized predictive scores
module simple_flex_predictive_tester
use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
use simple_defs,            only: dp, DPI
use simple_flex_predictive, only: flex_predictive_score, flex_predictive_mixture
use simple_test_utils,      only: assert_true, assert_false, assert_double
implicit none
private
public :: test_flex_predictive

contains

    subroutine test_flex_predictive()
        real(dp) :: empty(0,0), empty_b(0), g(1,1), b(1), prior(1,1), logp, expected, rank0, two_real
        real(dp) :: component(3), weights(3), resp(3), nan
        logical  :: ok
        nan = ieee_value(0._dp, ieee_quiet_nan)
        call flex_predictive_score(empty, empty_b, 0._dp, 0._dp, 1, empty, logp, ok)
        call assert_true(ok, 'rank-zero standard normal is valid')
        call assert_double(-0.91893853320467274178_dp, logp, 'known one-dimensional normalization')
        call flex_predictive_score(empty, empty_b, 2._dp, log(6._dp), 3, empty, rank0, ok)
        expected = -.5_dp*(3._dp*log(2._dp*DPI) + log(6._dp) + 2._dp)
        call assert_true(ok, 'rank-zero nonunit noise is valid')
        call assert_double(expected, rank0, 'rank zero retains noise determinant and real dimension')
        call flex_predictive_score(empty, empty_b, 2._dp, log(6._dp), 2, empty, two_real, ok)
        call assert_double(-.5_dp*log(2._dp*DPI), rank0-two_real, 'one extra real degree adds its constant')
        call flex_predictive_score(empty, empty_b, 5._dp, 0._dp, 2, empty, two_real, ok)
        expected = -log(2._dp*DPI) - 2.5_dp
        call assert_double(expected, two_real, 'one complex value with variance two has two unit-variance real parts')
        g = 4._dp
        b = 2._dp
        prior = 3._dp
        call flex_predictive_score(g, b, 1._dp, 0._dp, 1, prior, logp, ok)
        expected = -.5_dp*(log(2._dp*DPI) + log(37._dp) + 1._dp/37._dp)
        call assert_true(ok, 'nonunit prior is valid')
        call assert_double(expected, logp, 'one-dimensional integrated variance is thirty-seven', 64._dp)
        call test_dense_predictive()
        component = [logp, logp, nan]
        weights = [.3_dp,.7_dp,0._dp]
        call flex_predictive_mixture(component, weights, rank0, resp, ok)
        call assert_true(ok, 'duplicate components and zero weight are valid')
        call assert_double(logp, rank0, 'duplicating an identical component preserves density', 32._dp)
        call assert_true(maxval(abs(resp-weights)) < 1.e-14_dp, 'duplicate responsibilities equal their weights')
        component = [-10000._dp,-10002._dp, nan]
        weights = [.25_dp,.75_dp,0._dp]
        call flex_predictive_mixture(component, weights, rank0, resp, ok)
        expected = -10000._dp + log(.25_dp + .75_dp*exp(-2._dp))
        call assert_true(ok, 'log-sum-exp remains valid for tiny densities')
        call assert_double(expected, rank0, 'log-sum-exp uses normalized mixture probabilities', 8192._dp)
        call assert_double(1._dp, sum(resp), 'mixture responsibilities sum to one', 32._dp)
        weights = [1._dp,-.1_dp,.1_dp]
        call flex_predictive_mixture(component, weights, rank0, resp, ok)
        call assert_false(ok, 'negative mixture weights fail')
        call assert_true(ieee_is_nan(rank0), 'failed mixture returns no usable density')
        weights = [.5_dp,.6_dp,0._dp]
        call flex_predictive_mixture(component, weights, rank0, resp, ok)
        call assert_false(ok, 'unnormalized mixture weights fail')
        weights = 0._dp
        call flex_predictive_mixture(component, weights, rank0, resp, ok)
        call assert_false(ok, 'all-zero mixture weights fail')
        weights = [1._dp,0._dp,0._dp]
        component(1) = nan
        call flex_predictive_mixture(component, weights, rank0, resp, ok)
        call assert_false(ok, 'nonfinite active component fails')
        weights(1) = nan
        call flex_predictive_mixture(component, weights, rank0, resp, ok)
        call assert_false(ok, 'nonfinite weight fails')
        g = 0._dp
        b = 0._dp
        prior = 1._dp
        call flex_predictive_score(g, b, 1._dp, 0._dp, 1, prior, logp, ok)
        call assert_true(ok, 'unobserved latent direction with singular Gram is valid')
        call assert_double(-.5_dp*(log(2._dp*DPI)+1._dp), logp, 'unobserved latent direction integrates to one')
        call flex_predictive_score(g, empty_b, 1._dp, 0._dp, 1, prior, logp, ok)
        call assert_false(ok, 'mismatched statistic shapes fail before array access')
        g = -1._dp
        b = 0._dp
        prior = 1._dp
        call flex_predictive_score(g, b, 1._dp, 0._dp, 1, prior, logp, ok)
        call assert_false(ok, 'nonpositive posterior system fails without a ridge')
        g = -.5_dp
        call flex_predictive_score(g, b, 1._dp, 0._dp, 1, prior, logp, ok)
        call assert_false(ok, 'indefinite Gram fails even if posterior precision would be positive')
        g = 1._dp
        b = 2._dp
        call flex_predictive_score(g, b, 1._dp, 0._dp, 1, prior, logp, ok)
        call assert_false(ok, 'inconsistent joint statistics fail')
        b = 0._dp
        prior = 0._dp
        call flex_predictive_score(g, b, 1._dp, 0._dp, 1, prior, logp, ok)
        call assert_false(ok, 'singular prior must use an explicit reduced rank')
        prior = -1._dp
        call flex_predictive_score(g, b, 1._dp, 0._dp, 1, prior, logp, ok)
        call assert_false(ok, 'negative Cholesky diagonal fails')
        prior = 1._dp
        call flex_predictive_score(g, b, 1._dp, 0._dp, 0, prior, logp, ok)
        call assert_false(ok, 'zero observation dimension fails')
        call flex_predictive_score(g, b, -1._dp, 0._dp, 1, prior, logp, ok)
        call assert_false(ok, 'negative noise quadratic fails')
        call flex_predictive_score(g, b, 1._dp, nan, 1, prior, logp, ok)
        call assert_false(ok, 'nonfinite noise determinant fails')
        b = 1.e15_dp
        prior = 1.e15_dp
        call flex_predictive_score(g, b, 1.e30_dp, 0._dp, 1, prior, logp, ok)
        call assert_false(ok, 'catastrophic sufficient-statistic cancellation fails closed')
        call assert_true(ieee_is_nan(logp), 'failed component returns no usable density')
    end subroutine test_flex_predictive

    subroutine test_dense_predictive()
        real(dp) :: d(3,2), l(2,2), noise(3), residual(3), g(2,2), b(2), projected(3,2), c(3,3)
        real(dp) :: inverse(3,3), determinant, s, logp, expected
        integer  :: i, j
        logical  :: ok
        d(:,1) = [1._dp,-.5_dp,.75_dp]
        d(:,2) = [.2_dp,1.2_dp,-.3_dp]
        l(:,1) = [1.4_dp,.3_dp]
        l(:,2) = [0._dp,.7_dp]
        noise = [.5_dp,2._dp,1.5_dp]
        residual = [.4_dp,-1.1_dp,.7_dp]
        do i = 1, 2
            b(i) = sum(d(:,i)*residual/noise)
            do j = 1, 2
                g(i,j) = sum(d(:,i)*d(:,j)/noise)
            enddo
        enddo
        s = sum(residual**2/noise)
        projected = matmul(d,l)
        c = matmul(projected,transpose(projected))
        do i = 1, 3
            c(i,i) = c(i,i) + noise(i)
        enddo
        ! Three-dimensional cofactors are independent of the latent Cholesky/Woodbury implementation.
        inverse(1,1) = c(2,2)*c(3,3)-c(2,3)*c(3,2)
        inverse(1,2) = c(1,3)*c(3,2)-c(1,2)*c(3,3)
        inverse(1,3) = c(1,2)*c(2,3)-c(1,3)*c(2,2)
        inverse(2,1) = c(2,3)*c(3,1)-c(2,1)*c(3,3)
        inverse(2,2) = c(1,1)*c(3,3)-c(1,3)*c(3,1)
        inverse(2,3) = c(1,3)*c(2,1)-c(1,1)*c(2,3)
        inverse(3,1) = c(2,1)*c(3,2)-c(2,2)*c(3,1)
        inverse(3,2) = c(1,2)*c(3,1)-c(1,1)*c(3,2)
        inverse(3,3) = c(1,1)*c(2,2)-c(1,2)*c(2,1)
        determinant = dot_product(c(1,:),inverse(:,1))
        inverse = inverse/determinant
        expected = -.5_dp*(3._dp*log(2._dp*DPI) + log(determinant) + dot_product(residual,matmul(inverse,residual)))
        call flex_predictive_score(g, b, s, sum(log(noise)), 3, l, logp, ok)
        call assert_true(ok, 'dense reference statistics are valid')
        call assert_double(expected, logp, 'dense covariance equals latent sufficient-statistic score', 256._dp)
        l(1,2) = .1_dp
        call flex_predictive_score(g, b, s, sum(log(noise)), 3, l, logp, ok)
        call assert_false(ok, 'nontriangular latent Cholesky factor fails')
        l(1,2) = 0._dp
        g(1,2) = g(1,2) + .1_dp
        call flex_predictive_score(g, b, s, sum(log(noise)), 3, l, logp, ok)
        call assert_false(ok, 'asymmetric Gram fails')
    end subroutine test_dense_predictive

end module simple_flex_predictive_tester
