fit_pl <- function(
        manifest,        # manifest dataset (n x p)
        constraintMat,   # constraints for loadings matrix (p x q)
        corrFLAG = 0,    # estimate latent covariance matrix
        silent = F,      # silent mode
        ncores = 1,
        method = 'nlminb',
        maxiter = 1000,
        init = 'none',
        control =list()
){
    start_time <- Sys.time()
    RcppParallel::setThreadOptions(numThreads = ncores)

    #### PREPARING MODEL INPUT
    n <- nrow(manifest) ### number of subjects
    p <- ncol(manifest) ### number of items
    q <- ncol(constraintMat) ### number of latent variables

    categories <- apply(manifest, 2, max) + 1 ### number of categories in each item

    start_time_init <- Sys.time()
    lambda0_init <- c()
    s <- 0

    for (i in 1:length(categories)) {
        vec <- 1:(categories[i]-1)
        vec <- (vec -min(vec))/(max(vec)-min(vec))*(2)-1
        lambda0_init[(s + 1):(s + categories[i] - 1)] <- vec
        s <- s + categories[i] - 1
    }
    lambda_init = rep(0.5, sum(constraintMat))
    transformed_rhos_init = rep(0, q*(q-1)/2)
    if(init == 'cfa'){

    }

    end_time_init <- Sys.time()
    #### PREPARING FITTING FUNCTIONS
    if(silent == F) cat('Computing frequency table...')
    freq_tab <- pairs_freq(manifest, categories)
    if(silent == F) cat('Done!\n')

    # function for nll
    pair_nll <- function(par_vec){
        lambda0_ <- par_vec[1:length(lambda0_init)]
        lambda_ <- par_vec[(length(lambda0_init)+1):(length(lambda0_init)+length(lambda_init))]
        transformed_rhos_ <- par_vec[(length(lambda0_init)+1+length(lambda_init)):length(par_vec)]
        mod <-multiThread_completePairwise(
            Y = manifest,
            C_VEC = categories,
            A = constraintMat,
            FREQ = freq_tab,
            TAU = lambda0_,
            LAMBDA = lambda_,
            TRANSFORMED_RHOS = transformed_rhos_,
            CORRFLAG = corrFLAG,
            GRFLAG = 0,
            SILENTFLAG = 1
        )
        out <- mod$iter_nll
        return(out)
    }

    # function for gradient
    pair_gr <- function(par_vec){
        lambda0_ <- par_vec[1:length(lambda0_init)]
        lambda_ <- par_vec[(length(lambda0_init)+1):(length(lambda0_init)+length(lambda_init))]
        transformed_rhos_ <- par_vec[(length(lambda0_init)+1+length(lambda_init)):length(par_vec)]
        mod <-multiThread_completePairwise(
            Y = manifest,
            C_VEC = categories,
            A = constraintMat,
            FREQ = freq_tab,
            TAU = lambda0_,
            LAMBDA = lambda_,
            TRANSFORMED_RHOS = transformed_rhos_,
            CORRFLAG = corrFLAG,
            GRFLAG = 1,
            SILENTFLAG = 1
        )

        out <- mod$iter_ngradient
        return(out)
    }

    pair_H <- function(par_vec){
        lambda0_ <- par_vec[1:length(lambda0_init)]
        lambda_ <- par_vec[(length(lambda0_init)+1):(length(lambda0_init)+length(lambda_init))]
        transformed_rhos_ <- par_vec[(length(lambda0_init)+1+length(lambda_init)):length(par_vec)]
        mod <- multiThread_completePairwise(
            Y = manifest,
            C_VEC = categories,
            A = constraintMat,
            FREQ = freq_tab,
            TAU = lambda0_,
            LAMBDA = lambda_,
            TRANSFORMED_RHOS = transformed_rhos_,
            CORRFLAG = corrFLAG,
            GRFLAG = 2,
            SILENTFLAG = 1
        )

        out <- mod$H_approx
        return(out)
    }

    # parameter vector
    par <- c(lambda0_init, lambda_init, transformed_rhos_init)
    convergence <- ''
    #### FIT
    if(silent == F) cat('Estimating the model...')
    if(method == 'nlminb'){
        opt <- nlminb(start = par, objective = pair_nll, gradient = pair_gr, control =control)
        par <- opt$par
        convergence <- opt$convergence
        message <- opt$message
        obj <- opt$objective
    } else if(method == 'ucminf'){
        require(ucminf)
        opt <- ucminf::ucminf(par, pair_nll, pair_gr)
        par <- opt$par
        message <- opt$convergence
        obj <- opt$objective

    } else if(method == 'nls.lm'){
        require(minpack.lm)
        sol <- nls.lm(par = par, fn = pair_gr, control = list(nprint = 0, maxiter = maxiter, maxfev = 10000))
        par <- sol$par
        message <- sol$info
        obj <- NULL
    } else if(method == 'nleqslv'){
        require(nleqslv)
        sol <- nleqslv(par, pair_gr, method = "Broyden", global= "dbldog", xscalm = "auto", control = list(allowSingular = T, maxit = maxiter, trace = 0))
        par <- sol$x
        message <- sol$termcd
        obj <- NULL
    }else if(method == 'ucminf2'){
        require(ucminf)
        opt <- ucminf::ucminf(
            par, pair_nll, pair_gr,
            control = list(invhessian.lt=solve(pair_H(par))[lower.tri(pair_H(par),diag=TRUE)])
        )
        par <- opt$par
        message <- opt$message
        convergence <- opt$convergence
        obj <- opt$objective

    }
    if(silent == F) cat('Done!\n')
    #### REARRANGE RESULTS
    if(silent == F) cat('Storing results...')
    fit = list()

    # Fit message
    fit$convergence <- convergence
    fit$message <- message

    # intercept
    # intercepts
    fit$theta <- par
    intercepts <- par[1:length(lambda0_init)]
    fit$intercepts<- list()
    s <- 1
    for(i in 1:p){
        fit$intercepts[[i]] <- intercepts[s:(s + categories[i] - 2)]
        s <- s + categories[i] -1
    }

    # loadings
    lambda <- par[(length(lambda0_init)+1):(length(lambda0_init)+length(lambda_init))]
    fit$loadings <- constraintMat
    s = 1
    for(h in 1:q){
        for(j in 1:p){
            if(fit$loadings[j, h] != 0.0)
            {
                fit$loadings[j, h] = lambda[s]
                s = s+1
            }
        }
    }


    # correlations
    trhos <- par[(length(lambda0_init)+length(lambda_init)+1):length(par)]
    rhos <- (exp(2*trhos)-1)/(exp(2*trhos)+1)

    fit$factor_cov <- matrix(1, q, q)
    s = 1
    for( h in 1:q){
        for(j in 1:q){
            if(j > h)
            {
                fit$factor_cov[j, h] = rhos[s]
                fit$factor_cov[h, j] = rhos[s]
                s = s + 1
            }
        }
    }

    # nll
    fit$nll <- obj
    fit$last_iter <- opt$iterations


    # time
    fit$time_init <- difftime(end_time_init, start_time_init, units = ("secs"))[[1]]
    end_time <- Sys.time()
    fit$time <- difftime(end_time, start_time, units = ("secs"))[[1]]
    if(silent == F) cat('Done!\n')

    return(fit)


}
