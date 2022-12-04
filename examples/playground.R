seed <- 1
set.seed(seed)
p <- 20; q <- 2; n <- 5000
constrMat <- build_constrMat(p,q,'simple')
true_load <- gen_loadings(FIXED = NULL, CONSTRAINT_MAT = constrMat, SEED = seed)
true_tau <- c(-1.2, 0, 1.2)
true_latent <- matrix(.7,q,q); diag(true_latent) <- 1
true_theta <- get_theta(rep(true_tau, p), true_load, true_latent, cat, constrMat, 0)
manifest <- gen_URV_data(n, true_load, true_tau, true_latent)

categories <- apply(manifest, 2, max) + 1 ### number of categories in each item
lambda0_init <- c()
s <- 0

CONSTRMAT <- constrMat
for (i in 1:length(categories)) {
    vec <- 1:(categories[i]-1)
    vec <- (vec -min(vec))/(max(vec)-min(vec))*(2)-1
    lambda0_init[(s + 1):(s + categories[i] - 1)] <- vec
    s <- s + categories[i] - 1
}
lambda_init = rep(0.5, sum(CONSTRMAT))
transformed_rhos_init = rep(0, q*(q-1)/2)


PL_BFGS <- fit_pl(
    manifest = manifest,
    constraintMat = constrMat,
    method = 'ucminf2',
    corrFLAG = 1,
    silent = F,
    ncores = 1
)
PL_BFGS
############################

par <- c(lambda0_init, lambda_init, transformed_rhos_init)
freq_tab <- pairs_freq(manifest, categories)
dim(freq_tab)
freq_tab[,1:10]

# function for nll
pair_nll <- function(par_vec){
    lambda0_ <- par_vec[1:length(lambda0_init)]
    lambda_ <- par_vec[(length(lambda0_init)+1):(length(lambda0_init)+length(lambda_init))]
    transformed_rhos_ <- par_vec[(length(lambda0_init)+1+length(lambda_init)):length(par_vec)]
    mod <-multiThread_completePairwise(
        Y = manifest,
        C_VEC = categories,
        A = CONSTRMAT,
        FREQ = freq_tab,
        TAU = lambda0_,
        LAMBDA = lambda_,
        TRANSFORMED_RHOS = transformed_rhos_,
        CORRFLAG = 1,
        GRFLAG = 0,
        SILENTFLAG = 1
    )
    out <- mod$iter_nll
    return(out)
}

pair_nll(par)

# function for gradient
pair_gr <- function(par_vec){
    lambda0_ <- par_vec[1:length(lambda0_init)]
    lambda_ <- par_vec[(length(lambda0_init)+1):(length(lambda0_init)+length(lambda_init))]
    transformed_rhos_ <- par_vec[(length(lambda0_init)+1+length(lambda_init)):length(par_vec)]
    mod <-multiThread_completePairwise(
        Y = manifest,
        C_VEC = categories,
        A = CONSTRMAT,
        FREQ = freq_tab,
        TAU = lambda0_,
        LAMBDA = lambda_,
        TRANSFORMED_RHOS = transformed_rhos_,
        CORRFLAG = 1,
        GRFLAG = 1,
        SILENTFLAG = 1
    )

    out <- mod$iter_ngradient
    return(out)
}
pair_gr(par)

pair_H <- function(par_vec){
    lambda0_ <- par_vec[1:length(lambda0_init)]
    lambda_ <- par_vec[(length(lambda0_init)+1):(length(lambda0_init)+length(lambda_init))]
    transformed_rhos_ <- par_vec[(length(lambda0_init)+1+length(lambda_init)):length(par_vec)]
    mod <- multiThread_completePairwise(
        Y = manifest,
        C_VEC = categories,
        A = CONSTRMAT,
        FREQ = freq_tab,
        TAU = lambda0_,
        LAMBDA = lambda_,
        TRANSFORMED_RHOS = transformed_rhos_,
        CORRFLAG = 1,
        GRFLAG = 2,
        SILENTFLAG = 1
    )

    out <- mod$H_approx
    return(out)
}
pair_H(par)

library(RcppClock)
###########
mod <- plCFA(
    DATA = manifest,
    C_VEC = categories,
    CONSTRMAT = constrMat,
    TAU = lambda0_init,
    LAMBDA = lambda_init,
    TRANSFORMED_RHOS = transformed_rhos_init,
    CORRFLAG = 1,
    SILENTFLAG = 1,
    METHODFLAG = 0,
    PAIRS_PER_ITERATION = 1,
    PROB = 10/100,
    ETA = 7,
    BURN = 200,
    MAXT = 1000,
    TOLGRAD = 5e-3,
    TOLPAR = 1e-5,
    TOLOBJ = 1e-3,
    TOLCOUNT = 10,
    SEED = 123
)
summary(clock, units = 's')
mod
mod$scale
1/(n*100)
n_pairs <- (p*(p-1)/2)
weights <- sampling_step(
    FULL_POOL = 1:n_pairs,
    SAMPLINGFLAG = 1,
    PROB = .1,
    PAIRS_PER_ITERATION = round(n_pairs/10,0),
    N_ITEMS = p,
    SEED = 123,
    SILENTFLAG = 1,
    ITER = 1
)
weights
length(weights)

stepsize_grid <- seq(3, 20, 1)
nll_grid <- purrr::map_dbl(stepsize_grid, ~{
    obj <- plCFA(
        Y = manifest,
        C_VEC = categories,
        A = constrMat,
        TAU = lambda0_init,
        LAMBDA = lambda_init,
        TRANSFORMED_RHOS = transformed_rhos_init,
        CORRFLAG = 1,
        SILENTFLAG = 1,
        METHODFLAG = 0,
        PAIRS_PER_ITERATION = 1,
        PROB = 10/100,
        ETA = .x,
        BURN = 200,
        MAXT = 100,
        TOLGRAD = 5e-3,
        TOLPAR = 1e-5,
        TOLOBJ = 1e-3,
        TOLCOUNT = 10,
        SEED = 123
    )

    obj$path_nll[obj$last_iter]
})
nll_grid
stepsize_grid[which.min(nll_grid)]

########
cpp_ctrl <- list(
    MAXT = 1000,
    BURN = 100,
    ETA = 7.9,
    PAIRS_PER_ITERATION = 10,
    PROB = .1
)
obj <- fit_plCFA(
    DATA_LIST = list('DATA' = manifest, 'CONSTRMAT' = constrMat, 'CORRFLAG'=1),
    METHOD = 'st_hyper',
    CPP_CONTROL = cpp_ctrl,
    INIT = NULL
)
nll8 <- obj$fit$path_nll[length(obj$fit$path_nll)]
mse8 <- mean((obj$theta-true_theta)^2)

nll1; nll4; nll7; nll8
mse1; mse4; mse7; mse8
obj$args
obj$fit$path_av_theta[nrow(obj$fit$path_av_theta),]
mse8 <- mean((obj$theta-true_theta)^2)
mse7;mse8
obj$clock
get_Lam(theta = obj$theta, c_vec=categories, A = constrMat)
get_Sigma_u(theta = obj$theta, A = constrMat)
PL_BFGS$loadings


cpp_ctrl <- list(
    MAXT = 200,
    BURN = 100,
    PAIRS_PER_ITERATION = 10,
    PROB = .1
)
stepsize_tuning(
    DATA_LIST = list('DATA' = manifest, 'CONSTRMAT' = constrMat, 'CORRFLAG'=1),
    METHOD = 'st_hyper',
    CPP_CONTROL = cpp_ctrl,
    STEPSIZE_GRID = seq(7,9,.1),
    INIT = NULL
)
