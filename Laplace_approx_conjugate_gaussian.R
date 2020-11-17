conjugate_gaussian_laplace_mu_tau = function( a0, b0, mu0, lambda0, X ){
    N = length(X)
    xbar = mean(X)
    sigma2 = mean( (X - xbar)^2 )
    
    # Theta vector
    mu = (N * xbar + lambda0 * mu0)/
        ( N + lambda0 )
    tau = 2 * ( (N+1)/2 + a0 - 1 ) /
        ( N * sigma2 + N * (xbar - mu)^2 + lambda0 * (mu - mu0)^2 + 2*b0 )
    theta = c(mu = mu, tau = tau)
    # Hessian
    d_dmu2 = -tau * ( N + lambda0 )
    d_dtau2 = -(1/tau^2) * ( (N+1)/2 + a0 - 1 )
    ## cross derivatives are zero
    H = -diag( c(d_dmu2, d_dtau2 ) )
    vcov_ = solve( H )
    
    list(
        theta = theta, H = H, vcov_ = vcov_
    )
}

conjugate_gaussian_laplace_mu_logtau = function( a0, b0, mu0, lambda0, X ){
    N = length(X)
    xbar = mean(X)
    sigma2 = mean( (X - xbar)^2 )
    
    mu_mu02 = (mu - mu0)^2
    
    # Theta vector
    mu = (N * xbar + lambda0 * mu0)/
        ( N + lambda0 )
    log_tau = log(
        2 * ( (N+1)/2 + a0 - 1 ) /
        ( N * sigma2 + N * (xbar - mu)^2 + lambda0 * mu_mu02 + 2*b0 )
    )
    tau = exp(log_tau)
    theta = c(mu = mu, log_tau = log_tau)
    
    # Hessian
    d_dmu2 = -tau * ( N + lambda0 )
    d_dlogtau2 = -tau/2 * (
        N * sigma2 + N * (xbar - mu)^2 + lambda0 * mu_mu02 + 2*b0 
    )
    # ## Cross derivatives are zero
    H = -diag( c(d_dmu2, d_dlogtau2) )
    vcov_ = solve( H )

    list(
        theta = theta, H = H, vcov_ = vcov_
    )
}