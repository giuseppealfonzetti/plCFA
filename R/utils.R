#'@export
get_corr <- function(THETA, CONSTRMAT, FISHERFLAG = F){
    q <- ncol(CONSTRMAT)
    d <- length(THETA)
    ncorr <- q*(q-1)/2
    trho <- THETA[(d-ncorr+1):d]
    if(FISHERFLAG == F){
        out <- (exp(2*trho)-1)/(exp(2*trho)+1)
    }else if (FISHERFLAG==T){
        out <- trho
    }
    return(out)
}

#'@export
get_lambda <- function(THETA, CAT, CONSTRMAT){
    q <- ncol(CONSTRMAT)
    p <- nrow(CONSTRMAT)
    d <- length(THETA)
    c <- sum(CAT)

    ncorr <- q*(q-1)/2
    lambda <- THETA[(c-p+1):(d-ncorr)]
    lambda
}

#'@export
get_theta <- function(TAU, LOADINGS, LATENT_COV, CAT, A, TAUREPFLAG = 1){
    # thresholds params
    thr_vec <- c()
    if(TAUREPFLAG==1){
        item = 1; thr = 0
        for (i in 1:length(tau)) {
            maxitem = CAT[item]-1
            if(thr==0){
                thr_vec[i] = TAU[i]
            }else{
                thr_vec[i] = log(TAU[i] - TAU[i-1])
            }

            if(thr==maxitem-1){
                item = item + 1
                thr = 0
            }else{
                thr = thr + 1
            }

        }
    }else if(TAUREPFLAG==0){
        thr_vec <- TAU
    }


    load_vec <- c()
    s <- 1
    for (j in 1:ncol(LOADINGS)) {
        for (i in 1:nrow(LOADINGS)) {
            if(A[i,j]!=0){
                load_vec[s] <- LOADINGS[i,j]
                s = s+1
            }
        }
    }

    corr_vec <- c()
    s <- 1
    for(j in 1:ncol(LATENT_COV)){
        for (i in 1:nrow(LATENT_COV)){
            if(i>j){
                rho <- LATENT_COV[i,j]
                corr_vec[s] <- dplyr::if_else(abs(rho)<=1,.5*log((rho+1)/(1-rho)), 0)
                s = s+1
            }
        }
    }
    theta <- c(thr_vec, load_vec, corr_vec)
    return(theta)
}

#'@export
check_CSGD_args <- function(ARGS, N){

    out <- ARGS

    if(is.null(ARGS$MAXT)) out$MAXT <- round(N^1.1,0)
    if(is.null(ARGS$BURN)) out$BURN <- 0
    if(is.null(ARGS$ETA)) out$ETA <- 1e-3
    if(is.null(ARGS$SEED)) out$SEED <- 123
    if(is.null(ARGS$PAIRS_PER_ITERATION)) out$PAIRS_PER_ITERATION <- 1
    if(is.null(ARGS$PROB)) out$PROB <- .1
    return(out)
}

#'@export
get_tidy_path <- function(MOD_OBJ, PATH_LAB){
    iters <- MOD_OBJ$iterations_subset
    path  <- MOD_OBJ$fit[[PATH_LAB]]

    out <- dplyr::tibble(iter = iters) %>%
        dplyr::mutate(
            path_chosen = split(t(path), rep(1:nrow(path), each = ncol(path)))
        )

    colnames(out) <- c('iter', PATH_LAB)

    return(out)
}
