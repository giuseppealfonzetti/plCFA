#'@export
fit_plCFA <- function(
        DATA_LIST = list('DATA', 'CONSTRMAT', 'CORRFLAG'),
        METHOD = 'ucminf',
        CPP_CONTROL = list(),
        UCMINF_CONTROL = list('ctrl' = list(), 'hessian' = 0 ),
        INIT = NULL,
        ITERATIONS_SUBSET = NULL,
        VERBOSEFLAG = 0,
        NCORES = 1
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
        if(!(METHOD %in% c('ucminf','st_bernoulli', 'st_hyper'))) stop('Method not available.')
        out$method <- METHOD
        RcppParallel::setThreadOptions(numThreads = NCORES)


        # Numerical optimisation
        if(METHOD == 'ucminf'){

                message('2. Optimising with ucminf...')

                # Compute frequency table bivariate patterns
                freq_tab <- pairs_freq(DATA_LIST$DATA, categories)

                Rwr_ncl <- function(par_vec){
                        lambda0_ <- par_vec[1:length(lambda0_init)]
                        lambda_ <- par_vec[(length(lambda0_init)+1):(length(lambda0_init)+length(lambda_init))]
                        transformed_rhos_ <- par_vec[(length(lambda0_init)+1+length(lambda_init)):length(par_vec)]
                        mod <-multiThread_completePairwise(
                                Y = DATA_LIST$DATA,
                                C_VEC = categories,
                                A = DATA_LIST$CONSTRMAT,
                                FREQ = freq_tab,
                                TAU = lambda0_,
                                LAMBDA = lambda_,
                                TRANSFORMED_RHOS = transformed_rhos_,
                                CORRFLAG = DATA_LIST$CORRFLAG,
                                GRFLAG = 0,
                                SILENTFLAG = 1
                        )
                        out <- mod$iter_nll
                        return(out)
                }

                # function for gradient
                Rwr_ngr <- function(par_vec){
                        lambda0_ <- par_vec[1:length(lambda0_init)]
                        lambda_ <- par_vec[(length(lambda0_init)+1):(length(lambda0_init)+length(lambda_init))]
                        transformed_rhos_ <- par_vec[(length(lambda0_init)+1+length(lambda_init)):length(par_vec)]
                        mod <-multiThread_completePairwise(
                                Y = DATA_LIST$DATA,
                                C_VEC = categories,
                                A = DATA_LIST$CONSTRMAT,
                                FREQ = freq_tab,
                                TAU = lambda0_,
                                LAMBDA = lambda_,
                                TRANSFORMED_RHOS = transformed_rhos_,
                                CORRFLAG = DATA_LIST$CORRFLAG,
                                GRFLAG = 1,
                                SILENTFLAG = 1
                        )

                        out <- mod$iter_ngradient
                        return(out)
                }

                # list of ucminf args
                args <- list(
                        'par' = out$theta_init,
                        'fn' = Rwr_ncl,
                        'gr' = Rwr_ngr,
                        'control' = UCMINF_CONTROL$ctrl,
                        'hessian' = UCMINF_CONTROL$hessian)

                # optimisation
                opt <- do.call(ucminf::ucminf, args)
                out$fit <- opt

                out$control <- UCMINF_CONTROL
                out$theta   <- opt$par

                end_time <- Sys.time()
                out$time <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])
                message('3. Done! (', round(out$time,2),' secs)')
                return(out)
        }

        # Stochastic approximation of numerical optimiser
        if(METHOD == 'st_bernoulli' | METHOD == 'st_hyper' ){

                message(paste0('2. Optimising with ', METHOD, '...'))

                # Check stochastic control parameters
                cpp_ctrl <- check_CSGD_args(CPP_CONTROL, N = n)

                # Check iterations selected
                if(!is.null(ITERATIONS_SUBSET)){
                        out$iterations_subset <- unique(c(0, ITERATIONS_SUBSET[ITERATIONS_SUBSET < cpp_ctrl$MAXT], cpp_ctrl$MAXT))
                }else{
                        out$iterations_subset <- 0:cpp_ctrl$MAXT
                }


                # Guarantee reproducibility stochastic optimisation
                # Note: R set.seed() is only needed by Bernoulli sampling.
                # For Hypergeometric sampling the seed is directly passed to cpp
                set.seed(cpp_ctrl$SEED)

                # Collect and rearrange arguments to pass to cpp function
                args <- append(
                        list( 'TAU' = out$theta_init[1:(sum(out$categories)-p)],
                              'LAMBDA' = get_lambda(out$theta_init, out$categories, out$constraints),
                              'TRANSFORMED_RHOS' = get_corr(out$theta_init, out$constraints),
                              'C_VEC' = out$categories),
                        c( DATA_LIST, cpp_ctrl) )

                args$METHODFLAG <- dplyr::if_else(METHOD == 'st_hyper', 0, 1)

                fit <- do.call(plCFA, args)
                end_time <- Sys.time()
                out$time <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])

                out$iterations_subset <- c(out$iterations_subset, fit$last_iter)#[c(out$iterations_subset, fit$las_iter-1)<=(fit$las_iter-1)]
                out$iterations_subset <- out$iterations_subset[out$iterations_subset<=fit$last_iter]
                fit$path_theta    <- fit$path_theta[out$iterations_subset + 1,]
                fit$path_av_theta <- fit$path_av_theta[out$iterations_subset + 1,]
                fit$path_grad     <- fit$path_grad[out$iterations_subset,]

                out$control <- cpp_ctrl
                out$fit <- fit
                out$theta <- fit$path_av_theta[nrow(fit$path_av_theta),]

                if('RcppClock'%in% (.packages())) out$clock <- summary(clock, units = 's')
                message('\n3. Done! (', round(out$time,2),' secs)')


                return(out)
        }

}
