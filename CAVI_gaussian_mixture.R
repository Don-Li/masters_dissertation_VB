# Code for fitting VB Gaussian mixture
    # Function definitions after the loop
# Load data
old_faithful = data.table( faithful )
# Make environment
VB_gmm = new.env()
add_data_to_VB( VB_gmm, old_faithful )

# Set priors
priors = list(
    a0 = 1e-3,
    b0 = 1,
    m0 = c(3.5, 71),
    nu0 = VB_gmm$D + 1,
    w0 = diag( 1, VB_gmm$D )
)

settings = list(
    max_iter = 100,
    ELBO_tol = 1e-15,
    n_components = 10
)

add_list_to_VB( VB_gmm, priors )
add_list_to_VB( VB_gmm, settings )
setup_VB( VB_gmm )

# Initialise
    #"square" for initialise on a grid
    # There is also "kmeans" for initialising with k-meeans
GMM_VB_init( VB_gmm, init_option = "square" )

iter = 1

for ( iter in 1:VB_gmm$max_iter ){
    update_r_nk( VB_gmm )

    # Variational M step
    # Compute the soft statistics
    ## N_k, S_k, xbar_k
    update_soft_statistics( VB_gmm )
    # print( VB_gmm$m_k )

    # Variational E step
    ## Update Dirichlet;  a_k
    ## Update Gaussian-Wishart
    ### b_k, nu_k, W_k
    update_variational_params( VB_gmm )

    update_ELBO( VB_gmm, iter )
}


add_data_to_VB = function( VB_gmm, dataset ){
    VB_gmm$X = as.matrix( dataset )
    VB_gmm$D = ncol( dataset )
    VB_gmm$N = nrow( dataset )
}

add_list_to_VB = function( VB_gmm, list_ ){
    list2env( list_, envir = VB_gmm )
}

setup_VB = function( VB_gmm ){
    D = VB_gmm$D
    n_components = VB_gmm$n_components
    N = VB_gmm$N
    max_iter = VB_gmm$max_iter
    
    VB_gmm$W0_inv = solve( VB_gmm$w0 )
    
    # Store responsibilities
    r_nk = matrix( 0, N, n_components )
    log_r_nk = r_nk 
    log_rho_nk = r_nk
    
    # Soft summaries
    xbar_k = matrix( 0, n_components, D )
    S_k = array( 0, c(D, D, n_components) )
    N_k = rep( 0, n_components )
    
    # Probs and stuff
    expectation_log_pi = rep( 0, n_components )
    expectation_log_det_lambda = rep( 0, n_components )
    
    ELBO = rep( Inf, max_iter )
    
    VB_gmm$r_nk = r_nk
    VB_gmm$log_r_nk = log_r_nk
    VB_gmm$log_rho_nk = log_rho_nk
    VB_gmm$xbar_k = xbar_k
    VB_gmm$S_k = S_k
    VB_gmm$N_k = N_k
    VB_gmm$expectation_log_pi = expectation_log_pi
    VB_gmm$expectation_log_det_lambda = expectation_log_det_lambda
    VB_gmm$ELBO = ELBO
    VB_gmm
}

GMM_VB_init = function( VB_gmm, init_option = "kmeans", n_restart = 10 ){
    D = VB_gmm$D
    n_components = VB_gmm$n_components
    
    if ( init_option == "kmeans" ){
        init_m_k = kmeans( VB_gmm$X, n_components, nstart = n_restart, 
            iter.max = 1e5 )$centers
    }
    if ( init_option == "square" ){
        init_m_k = apply( VB_gmm$X, 2, function(Y){
            runif( n_components, min(Y), max(Y) )
        } )
    }
    if ( init_option == "point" ){
        NN = 1:nrow( VB_gmm$X )
        id = rep( sample( NN, size = 1 ), n_components )
        init_m_k = VB_gmm$X[ id, ]
        # init_m_k[] = 100
        init_m_k = init_m_k + runif( n_components * D, -1e-5, 1e-5 )
        
    }
    init_m_k = init_m_k[ do.call( order, as.data.table(init_m_k) ), ]
    VB_gmm$m_k = init_m_k
        
    VB_gmm$b_k = rep( VB_gmm$b0, n_components )
    nu_k = rep( VB_gmm$nu0, n_components )
    VB_gmm$nu_k = nu_k
    a_k = rep( VB_gmm$a0, n_components )
    VB_gmm$a_k = a_k
    
    W_k = array( 0, c(D,D,n_components) )
    W_k[] = VB_gmm$w0
    VB_gmm$W_k = W_k
    
    VB_gmm$expectation_log_pi[] = digamma(a_k) - digamma(sum(a_k))
    
    for ( comp in 1:n_components ){
        digamma_ = sum( digamma( (nu_k[comp] + 1 - 1:D)/2 ) )
        VB_gmm$expectation_log_det_lambda[comp] = digamma_ +
            D*log(2) + log(det(W_k[,,comp]))
    }
    VB_gmm
}

update_r_nk = function( VB_gmm ){
    D = VB_gmm$D
    n_components = VB_gmm$n_components
    N = VB_gmm$N
    max_iter = VB_gmm$max_iter
    
    comp = 1
    for ( comp in 1:n_components ){
        xn_diff_mk = VB_gmm$X - VB_gmm$m_k[rep(comp,N),]
        
        # E_{mu, lambda}[ (x_n - mu_k)^T %*% lambda_k %*% (x_n - mu_k) ]
        expectation_quadratic_form = 
            D / VB_gmm$b_k[comp] + 
                VB_gmm$nu_k[comp] * rowSums(  xn_diff_mk %*% VB_gmm$W_k[,,comp] *  xn_diff_mk )
        
        # E[ log det lambda_k ]
        digamma_term = sum( digamma( (VB_gmm$nu_k[comp] + 1 - 1:D)/2 ) )
        VB_gmm$expectation_log_det_lambda[comp] = digamma_term +
            D * log(2) + log( det( VB_gmm$W_k[,,comp] ) )
        
        # E[ log pi_k ]
        VB_gmm$expectation_log_pi[comp] = digamma( VB_gmm$a_k[comp] ) - digamma( sum( VB_gmm$a_k ) )
            
        VB_gmm$log_rho_nk[ , comp ] = VB_gmm$expectation_log_pi[comp] + 
            VB_gmm$expectation_log_det_lambda[comp]/2 -
            expectation_quadratic_form/2
        
    }
    
    VB_gmm$log_r_nk[] = VB_gmm$log_rho_nk - matrixStats::rowLogSumExps( VB_gmm$log_rho_nk )
    VB_gmm$r_nk[] = exp( VB_gmm$log_r_nk )

    VB_gmm
}

update_soft_statistics = function( VB_gmm ){
    D = VB_gmm$D
    n_components = VB_gmm$n_components
    N = VB_gmm$N
    max_iter = VB_gmm$max_iter
    zero_tol = 1e-15
    
    xbar_k_temp = matrix( 0, n_components, D )
    S_k_temp = array( 0, c(D,D,n_components) )
    N_k_temp = rep( 0, n_components )
    
    N_k_temp[] = colSums( VB_gmm$r_nk ) + zero_tol
    
    for ( comp in 1:n_components ){
        if ( N_k_temp[comp] >= zero_tol ){
            xbar_k_temp[comp,] = colSums( VB_gmm$r_nk[,comp] * VB_gmm$X )/N_k_temp[comp]
            x_soft_diff = VB_gmm$X - xbar_k_temp[ rep(comp,N), ]
            S_k_temp[,,comp] = t(x_soft_diff) %*% (VB_gmm$r_nk[,comp] * x_soft_diff)/N_k_temp[comp]
        }
    }
    
    VB_gmm$N_k[] = N_k_temp
    VB_gmm$S_k[] = S_k_temp
    VB_gmm$xbar_k[] = xbar_k_temp
    
    VB_gmm
}

update_variational_params = function( VB_gmm ){
    D = VB_gmm$D
    n_components = VB_gmm$n_components
    N = VB_gmm$N
    max_iter = VB_gmm$max_iter
    
    N_k = VB_gmm$N_k
    m0_matrix = matrix( VB_gmm$m0, nrow = n_components, ncol = D, byrow = T )
    
    a_k_temp = VB_gmm$a0 + N_k
    b_k_temp = VB_gmm$b0 + N_k
    m_k_temp = 1/b_k_temp * ( VB_gmm$b0 * m0_matrix + N_k * VB_gmm$xbar_k )
    nu_k_temp = VB_gmm$nu0 + N_k
    
    W_k_temp = array( 0, c(D,D,n_components) )
    for ( comp in 1:n_components ){
        xbar_k_diff = VB_gmm$xbar_k[comp,,drop = F] - VB_gmm$m0
        W_k_inv = VB_gmm$W0_inv +
            N_k[comp] * VB_gmm$S_k[,,comp] + 
            (VB_gmm$b0*N_k[comp])/(VB_gmm$b0 + N_k[comp]) * t(xbar_k_diff) %*% xbar_k_diff
        W_k_temp[,,comp] = solve( W_k_inv )
    }
    
    VB_gmm$a_k[] = a_k_temp
    VB_gmm$b_k[] = b_k_temp
    VB_gmm$nu_k[] = nu_k_temp
    VB_gmm$W_k[] = W_k_temp
    VB_gmm$m_k[] = m_k_temp
    
    VB_gmm
}

log_B_wishart = function( W, nu, D ){
    -nu/2 * log( det( W ) ) -
        ( nu*D/2*log(2) +
            D*(D-1)/4 * log(pi) +
            sum( lgamma( (nu + 1 - 1:D)/2 ) )
        )
}

update_ELBO = function( VB_gmm, iter ){
    D = VB_gmm$D
    n_components = VB_gmm$n_components
    N = VB_gmm$N
    max_iter = VB_gmm$max_iter
    
    N_k = VB_gmm$N_k
    m_k = VB_gmm$m_k
    expectation_log_det_lambda = VB_gmm$expectation_log_det_lambda
    W_k = VB_gmm$W_k
    nu_k = VB_gmm$nu_k
    b_k = VB_gmm$b_k
    S_k = VB_gmm$S_k
    a_k = VB_gmm$a_k
    
    expectation_log_pi = VB_gmm$expectation_log_pi
    
    a0 = VB_gmm$a0
    b0 = VB_gmm$b0
    m0 = VB_gmm$m0
    W0_inv = VB_gmm$W0_inv
    W0 = VB_gmm$w0
    nu0 = VB_gmm$nu0
    
    # Expectation of log likelihood X|Z,mu,lambda
    expectation_true_log_gaussian = 0
    for ( k in 1:n_components ){
        x_diff = VB_gmm$xbar_k[k,] - m_k[k,]
        expectation_true_log_gaussian = expectation_true_log_gaussian +
            N_k[k] * (
            # From the computation of r_{n,k}
            expectation_log_det_lambda[k] -
                D* 1/b_k[k] -
                nu_k[k] * sum(diag( S_k[,,k] %*% W_k[,,k] )) -
                nu_k[k] * t(x_diff) %*% W_k[,,k] %*% x_diff -
                D * log(2*pi)
        )
    }
    
    expectation_true_log_gaussian = expectation_true_log_gaussian/2
    
    # Expectation of log p(z|pi)
    expectation_true_log_z_given_pi = 
        sum( VB_gmm$r_nk * matrix( expectation_log_pi, nrow = N, ncol = n_components, byrow = T ) )
    
    # Expectation of log p(pi)
    log_dirichlet_constant = lgamma( n_components*a0 ) - n_components*lgamma(a0)
    expectation_true_log_pi = log_dirichlet_constant +
        (a0 - 1)*sum(expectation_log_pi)
    
    # Expectation of log p(mu, lambda)
    expectation_true_log_mu_lambda = 0
    expectation_true_log_mu_lambda_trace_sum = 0

    for ( k in 1:n_components ){
        m_diff = m_k[k,] - m0
        expectation_true_log_mu_lambda = expectation_true_log_mu_lambda +
            D * log( b0/(2*pi) ) +
                expectation_log_det_lambda[k] -
                (D * b0)/(b_k[k]) -
                b0 * nu_k[k] * m_diff %*% W_k[,,k] %*% m_diff
        
        expectation_true_log_mu_lambda_trace_sum = 
            expectation_true_log_mu_lambda_trace_sum +
            nu_k[k] * sum(diag( W0_inv %*% W_k[,,k] ) )
    }
    
    expectation_true_log_mu_lambda = expectation_true_log_mu_lambda/2
    expectation_true_log_mu_lambda_trace_sum = expectation_true_log_mu_lambda_trace_sum/2

    expectation_true_log_mu_lambda = expectation_true_log_mu_lambda +
        n_components * log_B_wishart( W0, nu0, D ) +
        (nu0 - D - 1)/2 * sum(expectation_log_det_lambda) -
        expectation_true_log_mu_lambda_trace_sum
    
    # Expectation of log q(z)
    expectation_VB_log_z = sum( VB_gmm$r_nk * VB_gmm$log_r_nk )
    
    # Expectation of log q(pi)
    log_dririchlet_constant = lgamma( sum(a_k) ) - sum( lgamma( a_k ) )
    expectation_VB_log_pi =
        log_dririchlet_constant +
        sum( (a_k - 1)*expectation_log_pi )
    
    # Expectation log q(mu,lambda)
    expectation_VB_log_mu_lambda  = 0
    for ( k in 1:n_components ){
        expectation_VB_log_mu_lambda = expectation_VB_log_mu_lambda +
            expectation_log_det_lambda[k]/2 +
                D/2 * log( b_k[k] / (2*pi) ) -
                D/2 -
                (
                    # this term is the entropy of the inverse Wishart
                    -log_B_wishart( W_k[,,k], nu_k[k], D ) -
                    (nu_k[k] - D - 1)/2 * expectation_log_det_lambda[k] +
                    nu_k[k] * D / 2
                )
    }
    
    ELBO = expectation_true_log_gaussian + expectation_true_log_z_given_pi +
        expectation_true_log_pi + expectation_true_log_mu_lambda -
        expectation_VB_log_z  - expectation_VB_log_pi - 
        expectation_VB_log_mu_lambda
    
    VB_gmm$ELBO[iter] = ELBO
    VB_gmm
}