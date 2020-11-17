# Main function for fitting the conjugate Gaussian VB approximation
    # Function definitions below
    # Methods and class definitions below that
fit_CAVI = function( prior, x, verbose = TRUE ){
    N = length(x)
    mu0 = prior["mu0"]
    lambda0 = prior["lambda0"]
    a0 = prior["a0"]
    b0 = prior["b0"]
    
    VB_model = VB_uniGauss( x,
        mu0 = mu0, lambda0 = lambda0, a0 = a0, b0 = b0,
        muN = mu0, lambdaN = lambda0, aN = a0, bN = b0,
    maxiter = 1000 )
    param_names = VB_model@container$VB_param_list
    # VB_env_checker( VB_model )
    VB_update_ELBO( VB_model, verbose = verbose )

    for ( i in 1:VB_model@container$maxiter ){
        conv = VB_check_convergence( VB_model )
        if ( conv ) break
        
        VB_update_iter( VB_model, verbose )
        VB_iter = VB_model@container$iter
    
        for ( p in VB_model@container$VB_param_list ){
            VB_update( VB_model, p, verbose )
        }
        VB_update_ELBO( VB_model, verbose )
    }
    
    return_ = as.data.table( mget( param_names, VB_model@container ) )
    list( params = return_, model  = VB_model )
}

# Some functions for computing things
uniGaussian_mu_update = function( mu0, lambda0, sum_x, N ){
    (lambda0 * mu0 + sum_x)/(lambda0 + N)
}

uniGaussian_lambda_update = function( aN, bN, lambda0, N ){
    (lambda0 + N)*(aN/bN)
}

uniGaussian_a_update = function( a0, N ){
    (N + 1)/2 + a0
}

uniGaussian_b_update = function( muN, lambdaN, mu0, lambda0, a0, b0, sum_x, sum_x2, N ){
    inner_sum = sum_x2 + lambda0 * mu0^2 -
        2 * muN * (sum_x + lambda0 * mu0) +
        (lambda0 + N) * (1/lambdaN + muN^2)
        0
    b0 + 1/2 * inner_sum
}

uniGaussian_ELBO = function( muN, lambdaN, aN, bN, a0, b0, lambda0, mu0, x, N ){
    weight_mu = ( sum(x) + lambda0 * mu0 )/( N + lambda0 )
    c_constants = -( N + lambda0 ) * ( (sum(x) + lambda0*mu0)/(N+lambda0) )^2 +
        sum( x^2 ) + lambda0 * mu0^2 + 2 * b0
    
    ( a0 + (N-1)/2 ) * ( digamma( aN ) - log( bN ) ) -
        1/2 * (aN/bN) * (
            ( N + lambda0 ) * ( 1/lambdaN + (muN - weight_mu)^2 ) +
            c_constants
        ) -
        1/2 * log( lambdaN ) +
        lgamma( aN ) -
        ( aN - 1 ) * digamma( aN ) -
        log( bN ) +
        aN
}

uniGaussian_rng = function( muN, lambdaN, aN, bN, size = 1e6 ){
    mu = rnorm( size, muN, 1/sqrt(lambdaN ) )
    tau = rgamma( size, aN, bN )
    list( mu = mu, tau = tau )
}

uniGaussian_FD_mu = function(r_mu, r_tau, mu0, lambda0, x, p_ = 1){
    N = length(x)
    weight_mu_p = ( p_ * sum(x) + lambda0 * mu0)/(p_ * N + lambda0)
    weight_mu = ( sum(x) + lambda0 * mu0)/(N + lambda0)
    
    d_dmu_p = -r_tau * ( p_ * N + lambda0 ) * ( r_mu - weight_mu_p )
    d_dmu_q = -mean(r_tau) * ( N + lambda0 ) * ( r_mu - weight_mu )
    
    mean( (d_dmu_p - d_dmu_q)^2 )
}

uniGaussian_FD_tau = function(r_mu, r_tau, mu0, lambda0, a0, b0, x, p_ = 1){
    N = length(x)
    
    d_dtau_p = 1/r_tau * ( ( p_ * N + 1 )/2 + a0 - 1 ) -
        1/2 * ( 
            p_ * sum(x^2) - 2 * p_ * sum(x) * r_mu + p_ * N * r_mu^2 +
            lambda0 * (r_mu - mu0)^2 +
            2 * b0
            )
    d_dtau_q = 1/r_tau * ( ( N + 1 )/2 + a0 - 1) -
        b0 -
        1/2 * mean(
            sum( x^2 ) - 2 * sum( x ) * r_mu + N * r_mu^2 +
            lambda0 * ( r_mu - mu0)^2
        )
    
    mean( (d_dtau_p - d_dtau_q)^2 )
}

uniGauss_FD = function( r_mu, r_tau, mu0, lambda0, a0, b0, x, p_ = 1 ){
    N = length(x)
    FD_mu = uniGaussian_FD_mu( r_mu = r_mu, r_tau = r_tau, mu0 = mu0, 
        lambda0 = lambda0, x = x, p_ = p_ )
    FD_tau = uniGaussian_FD_tau( r_mu = r_mu, r_tau = r_tau, mu0 = mu0, 
        lambda0 = lambda0, a0 = a0, b0 = b0, x = x, p_ = p_)
    FD = FD_mu + FD_tau
    
    c( mu = FD_mu, tau = FD_tau, FD = FD )
}

uniGauss_FD_analytical = function( lambdaN, aN, bN, lambda0, N){
    var_tau = aN/bN^2
    FD_mu = ( N + lambda0 )^2 * ( 1/lambdaN ) * var_tau
    FD_tau = 1/2 * ( N + lambda0 )^2 * 1/lambdaN^2
    FD = FD_mu + FD_tau
    
    c( mu = FD_mu, tau = FD_tau, FD = FD )
}

# Make some methods and classes
setClass( "VB",
    representation = representation(
        container = "environment",
        slot_info = "character"
    )
)

setGeneric( "VB_update_iter", function( object, verbose ) standardGeneric( "VB_update_iter" ) )

setMethod( "VB_update_iter",
    signature( "VB", "logical" ),
    function( object, verbose = FALSE ){
        envir = object@container
        envir$iter = envir$iter + 1
        
        if ( verbose ){
            msg = paste0( "Iteration: ", envir$iter )
            cat( msg, "\n" )
        }
        invisible()
} )

setGeneric( "VB_check_convergence",
    function(object) standardGeneric("VB_check_convergence") )

setGeneric( "VB_rng",
    function(object) standardGeneric("VB_rng") )

setClass( "VB_conjugate_univariate_Guassian", contains = "VB" )

setMethod( "$", signature("VB"),
    function( x, name ){
        x@container[[name]]
} )

VB_uniGauss = function( X, mu0, lambda0, a0, b0, muN, lambdaN, aN, bN,
    maxiter, tolerance = 1e-5 ){
    
    slot_info = c(
        a0 = "numeric", b0 = "numeric", mu0 = "numeric", lambda0 = "numeric",
        sum_x = "numeric", sum_x2 = "numeric", N = "integer",
        muN = "numeric", aN = "numeric", bN = "numeric", lambdaN = "numeric",
        muN_history = "numeric", lambdaN_history = "numeric",
        aN_history = "numeric", bN_history = "numeric",
        X = "numeric",
        FD_mu = "numeric", FD_tau = "numeric", FD = "numeric",
        ELBO = "numeric",
        VB_param_list = "character",
        iter = "numeric", maxiter = "numeric",
        tolerance = "numeric", convergence = "logical"
    )
    
    empty_vect = rep( NaN, maxiter+1 )
    inf_vect = rep( Inf, maxiter+1 )
    env_list = list(
        X = X, sum_x = sum( X ), sum_x2 = sum( X^2 ), N = length(X),
        ELBO = inf_vect,
        FD_mu = inf_vect, FD_tau = inf_vect, FD  = inf_vect,
        mu0 = mu0, lambda0 = lambda0, a0 = a0, b0 = b0,
        muN = muN, lambdaN = lambdaN, aN = aN, bN = bN,
        muN_history = empty_vect, lambdaN_history = empty_vect, 
        aN_history = empty_vect, bN_history = empty_vect,
        VB_param_list = c("muN", "lambdaN", "aN", "bN"),
        iter = 1, maxiter = maxiter,
        tolerance = tolerance,
        convergence = FALSE
    )
    
    vb_env = list2env( env_list )
    
    vb_env$muN_history[1] = muN
    vb_env$lambdaN_history[1] = lambdaN
    vb_env$aN_history[1] = aN
    vb_env$bN_history[1] = bN
    vb_env$ELBO[1] = NaN
    
    vb_obj = new("VB_conjugate_univariate_Guassian",
        container = vb_env,
        slot_info = slot_info)
    
    vb_obj
}

setGeneric( "VB_env_checker",
    function(object) standardGeneric("VB_env_checker")
)

VB_env_checker_ = function( object ){
    env_stuff = sapply( object@container, class )
    check_eq = setequal( names(env_stuff), names(VB_model@slot_info) )
    if ( ! check_eq ){
        cat( "Object invalid, missing slots.\n" )
        stop()
    }
    e = env_stuff[ order( names( env_stuff ) ) ]
    s = VB_model@slot_info[ order( names( VB_model@slot_info ) ) ]
    
    if ( any( e != s ) ){
        cat( "Object slot type mismatch.\n" )
        cat( names(e)[ e != s ], "\n " )
        stop()
    }
}

setMethod( "VB_env_checker", signature("VB"),
    VB_env_checker_
)

model_param_checker = function( container, param_name ){
    if ( ! param_name %in% container$VB_param_list ){
        msg = paste0( param_name, " is not a valid parameter for this model" )
        stop( msg )
    }
}

setGeneric( "VB_update", function( object, param_name, verbose ) standardGeneric("VB_update") )

uniGaussian_VB_update_ = function( object, param_name, verbose = FALSE ){
        envir = object@container
        model_param_checker( envir, param_name )

        if ( param_name == "muN" ){
            new_param_val = uniGaussian_mu_update( mu0 = envir$mu0, lambda0 = envir$lambda0,
                sum_x = envir$sum_x, N = envir$N )
        }
        if ( param_name == "lambdaN" ){
            new_param_val = uniGaussian_lambda_update( lambda0 = envir$lambda0,
                N = envir$N, aN = envir$aN, bN = envir$bN )
        }
        if ( param_name == "aN" ){
            new_param_val = uniGaussian_a_update( a0 = envir$a0, N = envir$N )
        }
        if ( param_name == "bN" ){
            new_param_val = uniGaussian_b_update( muN = envir$muN, lambdaN = envir$lambdaN, 
                mu0 = envir$mu0, lambda0 = envir$lambda0, 
                a0 = envir$a0, b0 = envir$b0, sum_x =  envir$sum_x, 
                sum_x2 = envir$sum_x2, N = envir$N )
        }
        
        if ( verbose ){
            msg = paste0( "Iteration: ", envir$iter, ". New ", param_name, ": ", new_param_val )
            cat( msg, "\n" )
            }
        
        envir[[param_name]] = new_param_val
        history_var = paste0( param_name, "_history" )
        envir[[history_var]][envir$iter]  = new_param_val
        
        invisible()
}

setMethod( "VB_update", 
    signature("VB_conjugate_univariate_Guassian", "character", "logical"),
    uniGaussian_VB_update_
     )



uniGaussian_ELBO_update_ = function( object, verbose = FALSE ){
    envir = object@container
    
    new_ELBO = uniGaussian_ELBO( muN = envir$muN, lambdaN = envir$lambdaN, 
        aN = envir$aN, bN = envir$bN, a0 = envir$a0, b0 = envir$b0, 
        lambda0 = envir$lambda0, mu0 = envir$mu0, x = envir$X,
        N = envir$N )
    envir[["ELBO"]][ envir$iter ] = new_ELBO
    
    if ( verbose ){
        msg = paste0( "Iteration: ", envir$iter, ". ELBO ", new_ELBO )
        cat( msg, "\n" )
    }
    invisible()
}

setGeneric( "VB_update_ELBO", function( object, verbose ) standardGeneric( "VB_update_ELBO") )

setMethod( "VB_update_ELBO",
    signature( "VB_conjugate_univariate_Guassian", "logical" ),
    uniGaussian_ELBO_update_
)

setMethod( "VB_check_convergence", signature("VB"),
    function(object){
        envir = object@container
        
        iter = envir$iter
        ELBO = envir$ELBO
        tolerance = envir$tolerance
        
        if ( iter > 2 ){
            if ( iter >= envir$maxiter ){
                cat( "Max iterations reached without convergence\n" )
            } else if ( (ELBO[iter-1] - ELBO[iter]) > tolerance ){
                cat( "Warning, ELBO has computationally decreased. Check model or rounding error\n" )
            }
        }
        if ( iter > 11 ){
            E1 = ELBO[ (iter-11):iter ]
            E2 = E1[length(E1)]
            conv = all( (E2 - E1) < tolerance )
            
            if (conv){
                cat( "Convergence at", iter, "\n" )
                envir$convergence = TRUE
            }
        }
        
        envir$convergence
} )

setMethod( "show", signature("VB_conjugate_univariate_Guassian"),
    function(object){
        envir = object@container
        param_vect = unlist( mget( envir$VB_param_list, envir ) )
        print( param_vect )
        conv = envir$convergence
        cat( "Convergence", conv, "\n" )
        ELBO = envir$ELBO[ envir$iter ]
        cat( "ELBO", ELBO, "\n" )
} )

setGeneric( "VB_rng", function(object, size = 1e6) standardGeneric("VB_rng") )

setMethod( "VB_rng", signature("VB_conjugate_univariate_Guassian"),
    function(object, size = 1e6 ){
        envir = object@container
        uniGaussian_rng( muN = envir$muN, lambdaN = envir$lambdaN, 
            aN = envir$aN, bN = envir$bN, size = size )
} )

setGeneric( "VB_Fisher_divergence", function( object, p_ = 1, ...) standardGeneric("VB_Fisher_divergence") )

setMethod( "VB_Fisher_divergence", signature("VB_conjugate_univariate_Guassian"),
    function( object, p_, analytical = FALSE ){
        envir = object@container
        rng_stuff = VB_rng( object, size = 1e6 )

        if ( analytical ){
            fisher_D = uniGauss_FD_analytical( envir$lambdaN, envir$aN, envir$bN, envir$lambda0,
                envir$N )
        } else{
            fisher_D = uniGauss_FD( r_mu = rng_stuff$mu, r_tau = rng_stuff$tau, 
                mu0 = envir$mu0, lambda0 = envir$lambda0, a0 = envir$a0, 
                b0 = envir$b0, x = envir$X, p_ = p_ )
        }
        fisher_D
} )


