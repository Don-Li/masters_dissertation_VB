
##### Fisher divergence for CAVI conjugate Gaussian #####
FD_CAVI_conjugate_guassian = function( x, a0, b0, mu0, lambda0, r_mu, r_tau, p_seq ){
    # r_mu is a sample of mu from the CAVI approximation over which we compute the 
        # Monte Carlo estimate of the Fisher divergence.
    # Same for r_tau.
    # p_seq is a sequence of powers for the power posterior.

    N = length(x)
    xbar = mean(x)
    sigma2 = sum( (x - xbar)^2 )/N
    
    deriv_p_mu_lik = r_tau * N * ( xbar - r_mu )
    deriv_p_mu_prior = - lambda0 * r_tau * ( r_mu - mu0 )
    deriv_p_tau_lik = N/2 * 1/r_tau - 1/2 * N * ( sigma2 + ( xbar - r_mu )^2 )
    deriv_p_tau_prior = 1/r_tau * ( 1/2 + a0 - 1 ) -
        1/2 * ( 
            lambda0 * ( r_mu - mu0 )^2 + 2 * b0
        )
    
    weight_mu = ( sum(x) + lambda0 * mu0)/(N + lambda0)
    
    d_dmu_q = -datalist$aN/datalist$bN * ( N + lambda0 ) * ( r_mu - datalist$muN )
    d_dtau_q = 1/r_tau * ( ( N + 1 )/2 + a0 - 1) -
        b0 -
        1/2 * (
            sum( x^2 ) - 2 * sum( x ) * datalist$muN + 
                N * (1/datalist$lambdaN + datalist$muN^2) +
            lambda0 * ( ( datalist$muN - mu0 )^2 + 1/datalist$lambdaN )
        )
    
    powers = sapply( p_seq, function(p_){
        deriv_p_mu = p_ * deriv_p_mu_lik + deriv_p_mu_prior
        deriv_p_tau = p_ * deriv_p_tau_lik + deriv_p_tau_prior
        
        mu_diff = mean( (d_dmu_q - deriv_p_mu)^2 )
        tau_diff = mean( (d_dtau_q - deriv_p_tau)^2 )
        
        c( mu = mu_diff, tau = tau_diff, FD = mu_diff + tau_diff )
    } )

    cbind( t( powers ), p = p_seq )
}

##### Fisher divergence for ADVI meanfield conjugate Gaussian #####
FD_ADVI_meanfield = function( x, a0, b0, mu0, lambda0, 
    qmu_mu, qmu_var, qlogtau_mu, qlogtau_var,
    r_mu, r_logtau, p_seq ){
    # r_mu is a sample of mu from the CAVI approximation over which we compute the 
    # Monte Carlo estimate of the Fisher divergence.
    # Same for r_tau.
    # p_seq is a sequence of powers for the power posterior.
    
    # qmu_mu is the variational mean for mu from the ADVI approximation
        # qmu_var is the same but for variational variance for mu
    # Corresponding defintiions for qlogtau_mu
    
    r_tau = exp( r_logtau )

    deriv_q_mu = d_dmu_log_norm( qmu_mu, qmu_var, r_mu )
    deriv_q_tau = d_dmu_log_lnorm( qlogtau_mu, qlogtau_var, r_tau )
    
    N = length(x)
    xbar = mean(x)
    sigma2 = sum( (x - xbar)^2 )/N
    
    deriv_p_mu_lik = r_tau * N * ( xbar - r_mu )
    deriv_p_mu_prior = - lambda0 * r_tau * ( r_mu - mu0 )
    
    deriv_p_tau_lik = N/2 * 1/r_tau - 1/2 * N * ( sigma2 + ( xbar - r_mu )^2 )
    deriv_p_tau_prior = 1/r_tau * ( 1/2 + a0 - 1 ) -
        1/2 * ( 
            lambda0 * ( r_mu - mu0 )^2 + 2 * b0
        )

    powers = sapply( p_seq, function(p_){
        deriv_p_mu = p_ * deriv_p_mu_lik + deriv_p_mu_prior
        deriv_p_tau = p_ * deriv_p_tau_lik + deriv_p_tau_prior
        mu_diff = mean( (deriv_p_mu - deriv_q_mu)^2 )
        tau_diff = mean( (deriv_p_tau - deriv_q_tau)^2 )
        c( mu = mu_diff, tau = tau_diff, FD = mu_diff + tau_diff )
    } )

    cbind( t(powers), p = p_seq )
}

##### Fisher divergence for ADVI fullrank conjugate Gaussian #####
FD_ADVI_fullrank = function( x, a0, b0, mu0, lambda0, 
    r_mu, r_logtau, q_mu, q_vcov,
    p_seq ){
    
    # q_mu is the vector of variational means, mu and logtau
    # q_vcov is the vector of variational variances with the same ordering as q_mu

    r_tau = exp( r_logtau )

    q_vcov = fullrank_results$parmm_vcov
    
    # Gradient for the approximation, factor in the Jacobian to transform
        # logtau to tau
    grad_q = {
        x_mu_logtau = cbind( r_mu, r_logtau )
        inv_q_vcov = solve( q_vcov )
        diff_matrix = x_mu_logtau - t(mu_)[ rep(1,nrow(x_mu_logtau)), ]

        mvn_grad = ( diff_matrix %*% inv_q_vcov )
        mvn_grad[,2] = (mvn_grad[,2] + 1)*(1/r_tau)
        -mvn_grad
    }
    
    N = length(x)
    xbar = mean(x)
    sigma2 = sum( (x - xbar)^2 )/N
    
    # Derivatives for the target posterior
    deriv_p_mu_lik = r_tau * N * ( xbar - r_mu )
    deriv_p_mu_prior = - lambda0 * r_tau * ( r_mu - mu0 )
    
    deriv_p_tau_lik = N/2 * 1/r_tau - 1/2 * N * ( sigma2 + ( xbar - r_mu )^2 )
    deriv_p_tau_prior = 1/r_tau * ( 1/2 + a0 - 1 ) -
        1/2 * ( 
            lambda0 * ( r_mu - mu0 )^2 + 2 * b0
        )

    powers = sapply( p_seq, function(p_){
        deriv_p_mu = p_ * deriv_p_mu_lik + deriv_p_mu_prior
        deriv_p_tau = p_ * deriv_p_tau_lik + deriv_p_tau_prior
        
        mu_diff = mean( (deriv_p_mu - grad_q[,1])^2 )
        tau_diff = mean( (deriv_p_tau - grad_q[,2])^2 )
        
        c( mu = mu_diff, tau = tau_diff, FD = mu_diff + tau_diff )
    } )
    
    cbind( t(powers), p = p_seq )
}

##### Fisher divergence for Laplace approx of conjugate Gaussian {mu, tau} #####
FD_Laplace_mu_tau = function( X, a0, b0, mu0, lambda0,
    r_mu, r_tau, qmu_mu, qmu_var, qtau_mu, qtau_var,
    p_seq ){
    N = length(X)
    
    deriv_q_mu = d_dmu_log_norm( qmu_mu, qmu_var, r_mu )
    deriv_q_tau = d_dmu_log_norm( qtau_mu, qtau_var, r_tau )
    
    x = X
    N = length(x)
    xbar = mean(x)
    sigma2 = sum( (x - xbar)^2 )/N
    
    deriv_p_mu_lik = r_tau * N * ( xbar - r_mu )
    deriv_p_mu_prior = - lambda0 * r_tau * ( r_mu - mu0 )
    
    deriv_p_tau_lik = N/2 * 1/r_tau - 1/2 * N * ( sigma2 + ( xbar - r_mu )^2 )
    deriv_p_tau_prior = 1/r_tau * ( 1/2 + a0 - 1 ) -
        1/2 * ( 
            lambda0 * ( r_mu - mu0 )^2 + 2 * b0
        )
    
    powers = sapply( p_seq, function(p_){
        deriv_p_mu = p_ * deriv_p_mu_lik + deriv_p_mu_prior
        deriv_p_tau = p_ * deriv_p_tau_lik + deriv_p_tau_prior
        
        mu_diff = mean( (deriv_p_mu - deriv_q_mu)^2 )
        tau_diff = mean( (deriv_p_tau - deriv_q_tau)^2 )
        
        c( mu = mu_diff, tau = tau_diff, FD = mu_diff + tau_diff )
    } )
    cbind( t(powers), p = p_seq )
}

##### Fisher divergence for Laplace approx of conjugate Gaussian {mu, logtau} #####
FD_Laplace_mu_logtau = function( X, a0, b0, mu0, lambda0, 
    r_mu, r_logtau, qmu_mu, qmu_var, qlogtau_mu, qlogtau_var,
    p_seq ){
    N = length(X)
    
    r_tau = exp( r_logtau )
    
    deriv_q_mu = d_dmu_log_norm( qmu_mu, qmu_var, r_mu )
    deriv_q_tau = d_dmu_log_lnorm( qlogtau_mu, qlogtau_var, r_tau )
    
    x = X
    N = length(x)
    xbar = mean(x)
    sigma2 = sum( (x - xbar)^2 )/N
    
    deriv_p_mu_lik = r_tau * N * ( xbar - r_mu )
    deriv_p_mu_prior = - lambda0 * r_tau * ( r_mu - mu0 )
    
    deriv_p_tau_lik = N/2 * 1/r_tau - 1/2 * N * ( sigma2 + ( xbar - r_mu )^2 )
    deriv_p_tau_prior = 1/r_tau * ( 1/2 + a0 - 1 ) -
        1/2 * ( 
            lambda0 * ( r_mu - mu0 )^2 + 2 * b0
        )
    
    powers = sapply( p_seq, function(p_){
        deriv_p_mu = p_ * deriv_p_mu_lik + deriv_p_mu_prior
        deriv_p_tau = p_ * deriv_p_tau_lik + deriv_p_tau_prior
        
        mu_diff = mean( (deriv_p_mu - deriv_q_mu)^2 )
        tau_diff = mean( (deriv_p_tau - deriv_q_tau)^2 )
        
        c( mu = mu_diff, tau = tau_diff, FD = mu_diff + tau_diff )
    } )
    cbind( t(powers), p = p_seq )
}