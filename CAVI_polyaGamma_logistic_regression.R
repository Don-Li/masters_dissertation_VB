library( data.table )
library( mvtnorm )

source( "logistic_regression_utilities.R" )

logistic_regression_CAVI = function( X, Y, mu0, sigma0, maxiter = 100, beta_init = NULL,
    z_init = NULL ){
    
    N = length(Y)
    D = ncol(X)
    lambda0 = solve( sigma0 )
    logdet_sigma0 = log(det( lambda0 ))
    
    storage = data.table( mu = rep( list(list()), maxiter+1) )
    
    # The covariance matrix storage
    storage[ , sigmaN := list(list()) ]
    # The latent zeta for data augmentation
    storage[ , Z := list(list()) ]
    storage[ , ELBO := NaN ]
    
    # Initialisations
    if ( is.null(beta_init) ){
        beta = rep( 0, D )
    }
    if ( is.null(z_init) ){
        z_init = 0.25
    }
    Z = rep( z_init, N )
    Y_05 = Y - 0.5
    Xt = t(X)
    
    for ( iter in 1:(maxiter+1) ){
        # Update global params
        lambdaN = t(X * c(Z)) %*% X + lambda0
        logdet_sigmaN = log(det( lambdaN ))
        sigmaN = solve( lambdaN )
        
        muN = sigmaN %*% ( Xt %*% Y_05 + lambda0 %*% mu0 )
        
        # Update xi
        eta = X %*% muN
        xi = sqrt( eta^2 + rowSums( X %*% sigmaN * X ) )
        Z = tanh(xi/2)/(2*xi)
        # Durante and Rigon (2019) recommend re-initialising Z if something goes wrong
            # with it. Not sure when this happens, but Z is bounded 0 to 0.25.
        Z[ is.na(Z) ] = z_init
        
        # Update the ELBO if you want
        if ( iter %% 1000 == 0 ){
            muN_mu0 = muN - mu0
            ELBO = 1/2 *(
                - logdet_sigma0 +
                logdet_sigmaN -
                t(muN_mu0) %*% lambda0 %*% (muN_mu0) +
                2 * sum( Y_05 * eta + log(sigmoid(xi)) - 0.5 * xi ) -
                sum( diag(lambda0 %*% sigmaN) )
            )
            ELBO = as.numeric(ELBO)
        } else{
            ELBO = NaN
        }
        v = list( list(muN), list( sigmaN ), list(Z), ELBO )
        set( storage, i = iter, names(storage), v )
    }
    
    list(
        ELBO = ELBO,
        muN = muN, sigmaN = sigmaN, lambdaN = lambdaN, Z = Z,
        mu0 = mu0, sigma0 = sigma0,
        storage = storage
    )
}
























