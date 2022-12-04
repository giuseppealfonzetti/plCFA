#'@export
stepsize_tuning <- function(
        DATA_LIST = list('DATA', 'CONSTRMAT', 'CORRFLAG'),
        METHOD = 'st_hyper',
        CPP_CONTROL = list(),
        INIT = NULL,
        STEPSIZE_GRID = NULL,
        VERBOSEFLAG = 0
){

    out <- list()
    start_time <- Sys.time()
    # Identify model dimensions
    p <- ncol(DATA_LIST$DATA)
    n <- nrow(DATA_LIST$DATA)
    q <- ncol(DATA_LIST$CONSTRMAT)
    categories <- apply(DATA_LIST$DATA, 2, max) + 1 ### number of categories in each item
    d = sum(categories)-p + sum(constrMat) + q*(q-1)/2

    out$categories <- categories
    # Check Constraints
    if(is.null(DATA_LIST$CONSTRMAT)){
        stop('CONSTRMAT not declared')
    }else if(nrow(DATA_LIST$CONSTRMAT)!=ncol(DATA_LIST$DATA) || ncol(DATA_LIST$CONSTRMAT)>=nrow(DATA_LIST$CONSTRMAT)){

        stop('CONSTRMAT dimensions not acceptable. Check Items x Factors.')
    }
    out$constraints <- DATA_LIST$CONSTRMAT

    # Check Initialisation
    if(is.vector(INIT)){
        if(length(INIT)!=d)
            stop(paste0('init vector has length ', length(INIT), ' instead of ', d ,'.'))
        else
            message('1. Initialising at init vector.')
        out$theta_init <-  INIT
    }else{
        if(is.null(INIT))
            message('1. Initialising at default values')
        lambda0_init <- c()
        s <- 0

        for (i in 1:length(categories)) {
            vec <- 1:(categories[i]-1)
            vec <- (vec -min(vec))/(max(vec)-min(vec))*(2)-1
            lambda0_init[(s + 1):(s + categories[i] - 1)] <- vec
            s <- s + categories[i] - 1
        }
        lambda_init = rep(0.5, sum(DATA_LIST$CONSTRMAT))
        transformed_rhos_init = rep(0, q*(q-1)/2)
        out$theta_init <-  c(lambda0_init, lambda_init, transformed_rhos_init)
    }

    # Check if method entry is correct
    if(!(METHOD %in% c('st_bernoulli', 'st_hyper'))) stop('Method not available.')
    out$method <- METHOD

    # Check stochastic control parameters
    cpp_ctrl <- check_CSGD_args(CPP_CONTROL, N = n)

    # Collect and rearrange arguments to pass to cpp function
    args <- append(
        list( 'TAU' = out$theta_init[1:(sum(out$categories)-p)],
              'LAMBDA' = get_lambda(out$theta_init, out$categories, out$constraints),
              'TRANSFORMED_RHOS' = get_corr(out$theta_init, out$constraints),
              'C_VEC' = out$categories),
        c( DATA_LIST, cpp_ctrl) )

    args$METHODFLAG <- dplyr::if_else(METHOD == 'st_hyper', 0, 1)

    fun_grid <- sapply(STEPSIZE_GRID, function(x){
        args$ETA <- x
        fit <- do.call(plCFA, args)
        fun <- fit$path_nll[length(fit$path_nll)]
        return(fun)
        })

    out <- list(
        'tab' = cbind('stepsize' = STEPSIZE_GRID, 'ncl' = fun_grid),
        'chosen_eta' = STEPSIZE_GRID[which.min(fun_grid)]
    )

    return(out)

}
