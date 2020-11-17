source( "logistic_regression_utilities.R" )

#### ELBO gradients #####
elbo_gradient_meanfield = function(eta, mu, w, x, y, tx, mu0, sigma0_inv){
    dim_ = ncol(x)
    sigma = exp(w)
    zeta = sigma * eta + mu
    p = sigmoid(x %*% zeta)
    grad_mu = tx %*% (y - p) - sigma0_inv %*% ( zeta - mu0 )
    grad_w = sigma * eta * grad_mu + 1
    c(mu = grad_mu, w = grad_w)
}

elbo_gradient_fullrank = function( eta, mu, L, x, y, tx, mu0, sigma0_inv ){
    dim_ = ncol(x)
    dim_L = length(L)
    
    L2 = L
    diag(L2) = exp(diag(L2))
    zeta = L2 %*% eta + mu
    p = sigmoid(x %*% zeta)
    grad_mu = tx %*% ( y - p ) - sigma0_inv %*% ( zeta - mu0 )

    grad_L = outer( c(grad_mu), eta )
    grad_L = grad_L + t(solve(L2))
    list( mu = grad_mu, L = grad_L )
}

##### ELBO #####
elbo_meanfield = function( n_ELBO, mu, w, x, y, mu0, sigma0_inv ){
    dim_ = ncol(x)
    
    elbo_sample = sapply( 1:n_ELBO, function(i){
        r_eta = rnorm( dim_ )
        sigma = exp(w)
        zeta = sigma * r_eta + mu
        p = sigmoid( x %*% zeta )
        sum( dbinom( y, 1, p, log = T ) ) - 
            1/2 * ( zeta - mu0 )^T %*% sigma0_inv %*% ( zeta - mu0 )
    } )
    
    gaussian_entropy = 1 + log( 2 * pi ) + sum(w)
    mean(elbo_sample) + gaussian_entropy
}

elbo_fullrank = function( n_ELBO, mu, L, x, y, mu0, sigma0_inv  ){
    dim_ = ncol(x)
    dim_L = length(L)
    L2 = L
    diag(L2) = exp(diag(L2))

    elbo_sample = sapply( 1:n_ELBO, function(i){
        r_eta = rnorm( dim_ )
        zeta = L2 %*%  r_eta + mu
        p = sigmoid( x %*% zeta )
        sum( dbinom( y, 1, p, log = T ) ) - 
            1/2 * ( c(zeta) - mu0 )^T %*% sigma0_inv %*% ( c(zeta) - mu0 )
    } )
    gaussian_entropy = 0.5 * dim_ * ( 1 + log(2*pi) ) +
        0.5 * log( det( L2 %*% t(L2) ) )
    mean(elbo_sample) + gaussian_entropy
}

##### Tune the scale for gradient descent #####
tune_scale = function( x, y, n_tune_iter = 1000, seed_mu = 0, seed_w = 1, scale_sequence,
    n_ELBO = 1000, algorithm = "meanfield", mu0, sigma0, seed_L = NULL ){
    
    test_elbo_seq = sapply( scale_sequence, function( scale_factor ){
        cat( "-------------------\n" )
        cat( "Scale tuning:", scale_factor, "\n" )
        test_returns = try({
            logistic_advi( x = x, y = y, max_iter = n_tune_iter, elbo_eval_iter = n_tune_iter, 
            seed_mu = seed_mu, seed_w = seed_w, scale_factor = scale_factor, n_ELBO = n_ELBO,
            algorithm = algorithm, mu0 = mu0, sigma0 = sigma0, seed_L = seed_L )
        })
        if ( "try-error" %in% class(test_returns) ){
            test_ELBO = -Inf
        } else{
            test_ELBO = test_returns$elbo_storage[n_tune_iter]
        }
        cat( "ELBO:", test_ELBO, "\n" )
        test_ELBO
    } )
    best_scale = scale_sequence[ which.max(test_elbo_seq ) ]

    cat( "#####\n" )
    cat( "Using best scale", best_scale, "\n" )
    cat( "--------------------\n" )
    best_scale
}

# Convert a vector to a lower-triangle matrix
vect_to_lower_tri = function( X_vect ){
    cols = 1/2 * ( sqrt( 8 *length(X_vect) + 1 ) - 1 )
    zeros = matrix( 0, ncol = cols, nrow = cols )
    zeros[ lower.tri(zeros, diag =T ) ] = X_vect
    zeros
}

##### ADVI #####
logistic_advi = function( x, y, max_iter = 1e5, elbo_eval_iter = 100, seed_mu = 0, seed_w = 0, 
    seed_L = NULL, n_ELBO = 1, scale_factor, n_grad_samples = 50,
    print_iter = 500, algorithm = "meanfield", mu0, sigma0 ){
    
    dim_ = ncol( x )
    tx = t(x)
    sigma0_inv = solve( sigma0 )
    
    # means
    mu_storage = matrix(NaN, max_iter + 1, ncol = dim_ )
    colnames(mu_storage) = paste0( "mu", 1:dim_ )
    mu_storage[1,] = seed_mu
    
    if ( algorithm == "meanfield" ){
        # log-variances, mean-field
        w_storage = matrix(NaN, max_iter + 1, ncol = dim_ )
        colnames(w_storage) = paste0( "w", 1:dim_ )
        w_storage[1,] = seed_w
        
        # rho is the stochastic gradient ascent step size
        rho = matrix(NaN, max_iter, ncol = 2 * dim_ )
        # s is used for the step size
        s = matrix(NaN, max_iter, ncol = 2 * dim_ )
        # w is the log-variance
        w = w_storage[1,]
    } else{
        if ( is.null( seed_L) ) {
            L = diag(dim_)
        } else{
            L = seed_L
        }
        L_storage = matrix(NaN, max_iter + 1, dim_ * (dim_+1)/2 )
        L_storage[1,] = L[ lower.tri( L, diag = T ) ]
        
        # rho is the stochastic gradient ascent step size
        rho = matrix(NaN, max_iter, ncol = dim_ * (dim_+1)/2 + dim_ )
        # s is used for the step size
        s = matrix(NaN, max_iter, ncol = dim_ * (dim_+1)/2 + dim_ )
    }
    
    elbo_storage = rep(NaN, max_iter)
    
    # Do some things to make the containers
    if ( algorithm == "meanfield" ){
        test_w = rep(0, dim_)
        test_grad = elbo_gradient( rep(0,dim_), rep(0,dim_), test_w, 
            x[1,, drop = F], y[1], tx[,1,drop = F], mu0 = mu0, sigma0_inv = sigma0_inv )
    } else{
        test_L = diag( rep(0.9, dim_ ) )
        test_grad = elbo_gradient_fullrank( rep(0,dim_), rep(0,dim_), test_L, 
            x[1,, drop = F], y[1], tx[,1,drop = F], mu0, sigma0_inv )
        test_grad = c( 
            mu = test_grad$mu,
            L = test_grad$L[ lower.tri( test_grad$L, diag = T ) ]
        )
    }
    
    mu_names = grepl( "mu", names(test_grad) )
    w_names = grepl( "w", names(test_grad) )
    L_names = grepl( "L", names(test_grad) )
    
    i = 1
    for ( i in 1:max_iter ){
        if ( i %% print_iter == 0 ){
            cat( "ADVI iteration:", i, "\n" )
        }
        eta = rnorm( dim_ )
        mu = mu_storage[i, ]
        
        if ( algorithm == "meanfield" ){
            eta_ = matrix( rnorm( dim_ * n_grad_samples ), nrow = dim_ )
            grad_vector = c( rowMeans( apply( eta_, 2, function(eta ){
                grad_vector = elbo_gradient(eta, mu, w, x, y, tx, mu0 = mu0, sigma0_inv = sigma0_inv)
            } ) ) )
        } else{
            eta_ = matrix( rnorm( dim_ * n_grad_samples ), nrow = dim_ )
            grad_vector = c(rowMeans( apply( eta_, 2, function( eta ){
                grad_vector = elbo_gradient_fullrank(eta, mu, L, x, y, tx, mu0, sigma0_inv)
                grad_vector = c(
                    mu = grad_vector$mu,
                    L = grad_vector$L[ lower.tri( grad_vector$L, diag = T ) ]
                )
            } ) ))
        }
        
        # 0.1 from Robbins Monroe sequence, weighting factor
        s[i,] = 0.1 * grad_vector^2 + 
            ifelse( i > 1, 0.9 * s[i-1,], 0 ) 
        
        # Adjust SGD step size
        # See equation 10 in Kucukelbir et al
        # epsilon = 1e-16 to satisfy Robbins Monroe sequence
        rho[i,] = scale_factor * i^(-0.5 + 1e-16)/(1 + sqrt(s[i,]))
        
        # Update the vectors
        mu = mu + rho[ i, mu_names] * grad_vector[ mu_names ]
        mu_storage[ i+1, ] = mu
        
        # Update the variances
        if ( algorithm == "meanfield" ){
            w = w + rho[ i, w_names ] * grad_vector[ w_names ]
            w_storage[ i+1, ] = w
        } else{
            delta_L = vect_to_lower_tri( rho[ i, L_names ] * grad_vector[ L_names ] )
            L = L + delta_L
            L_storage[ i+1, ] = L[ lower.tri(L, diag = T) ]
        }
        
        # Estimate the ELBO
        if ( (i %% elbo_eval_iter) == 0 ){
            if ( algorithm == "meanfield" ){
                elbo_storage[i] = elbo(n_ELBO = n_ELBO, mu = mu, w = w, x = x, y = y,
                    mu0 = mu0, sigma0_inv = sigma0_inv )
            } else{
                elbo_storage[i] = elbo_fullrank( n_ELBO = n_ELBO, mu = mu, L = L, 
                    x = x, y = y, mu0 = mu0, sigma0_inv = sigma0_inv )
            }
        } 
    }
    
    if ( algorithm == "meanfield" ){
        return_list = list( mu = mu, w = w,
            mu_storage = mu_storage, w_storage = w_storage,
            elbo_storage = elbo_storage )
    } else{
        return_list = list( mu = mu, L = L,
            mu_storage = mu_storage, L_storage = L_storage,
            elbo_storage = elbo_storage )
    }
    
    return_list
}