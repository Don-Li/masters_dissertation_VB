logistic_regression_grad = function(beta, X, tX, Y, mu0, lambda0, lambda0inv, det_lambda0){
    mu = sigmoid( X %*% beta )
    k = length(beta)

    # c() to make 1-column matrix as vector, since R won't recycle that anymore
    loglik = c(
        sum( dbinom( Y, 1, mu, log = T ) ) +
        -k/2 * log( 2 * pi ) + 1/2 * log( det_lambda0 ) -
        1/2 * t(beta - mu0) %*% lambda0inv %*% (beta - mu0)
    )
    grad = tX %*% ( Y - mu ) -
        lambda0 %*% (beta  - mu0)
    hess = - t(X * c(mu*(1-mu))) %*% X -
        lambda0

    list( loglik = loglik,
        grad = c(grad),
        hess = hess
    )
}

logistic_regression_laplace_gradient_ascent = function( beta_, mu0, lambda0, X, Y, step_size, 
    maxiter = 100, tol = 1e-10,
    param_tol = Inf, newton = FALSE ){
    
    lambda0inv = solve( lambda0 )
    det_lambda0 = det( lambda0 )
    tX = t(X)
    
    grad_info = logistic_regression_grad( beta_, X, tX, Y, mu0, lambda0, lambda0inv, det_lambda0 )
    
    old_loglik = grad_info$loglik
    monitor = matrix( c(beta_, grad_info$loglik), ncol = length(beta_)+1,
        nrow = maxiter+1, byrow = T )
    nbeta = length(beta_)
    
    for ( i in 2:(maxiter+1) ){
        cat( "Iteration:", i, "\n" )
        if ( !newton ){
            beta_ = beta_ + step_size * grad_info$grad
        } else{
            beta_ = beta_ - c(solve(grad_info$hess) %*% grad_info$grad)
        }
        grad_info = logistic_regression_grad( beta_, X, tX, Y, mu0, lambda0, lambda0inv, det_lambda0 )
        monitor[ i, ] = c(beta_, grad_info$loglik)
        
        if ( (grad_info$loglik - old_loglik) < tol &
                all( abs(monitor[ i-1, 1:nbeta] - beta_) < param_tol )
            ){
            break
        }
        old_loglik = grad_info$loglik
    }
    
    monitor = monitor[ 1:i,]
    
    list(
        betaN = beta_,
        loglik = grad_info$loglik,
        grad = grad_info$grad,
        hess = grad_info$hess,
        monitor = monitor,
        vcov_ = solve(-grad_info$hess)
    )
}