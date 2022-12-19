seed <- 1
set.seed(seed)
p <- 40; q <- 4; n <- 2000
constrMat <- build_constrMat(p,q,'simple')
true_load <- gen_loadings(FIXED = NULL, CONSTRAINT_MAT = constrMat, SEED = seed)
true_tau <- c(-1.2, 0, 1.2)

mat <- matrix(runif(q*q, 0,1),q,q)
mat[!lower.tri(mat)] <- 0
mat <- mat + t(mat)
diag(mat) <- 1
true_latent <- mat
#true_latent <- matrix(runif(q*q, 0,1),q,q); diag(true_latent) <- 1
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


PL_BFGS <- fit_plCFA(
    DATA_LIST = list('DATA' = manifest, 'CONSTRMAT' = constrMat, 'CORRFLAG'=1),
    METHOD = 'ucminf',
    CPP_CONTROL = cpp_ctrl,
    INIT = NULL
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
    METHODFLAG = 1,
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
# cpp_ctrl <- list(
#     MAXT = 3000,
#     BURN = 100,
#     PAIRS_PER_ITERATION = 1,
#     PROB = .1
# )
# stepsize_tuning(
#     DATA_LIST = list('DATA' = manifest, 'CONSTRMAT' = constrMat, 'CORRFLAG'=1),
#     METHOD = 'st_hyper',
#     CPP_CONTROL = cpp_ctrl,
#     STEPSIZE_GRID = 5:15,
#     INIT = NULL
# )

K = p*(p-1)/2; K
cpp_ctrl <- list(
    MAXT = K*20,
    BURN = 500,
    ETA = 7,
    PAIRS_PER_ITERATION = 1,
    PROB = 5/K
)
my_seq <- seq(0, cpp_ctrl$MAXT, 100)
st_1 <- fit_plCFA(
    DATA_LIST = list('DATA' = manifest, 'CONSTRMAT' = constrMat, 'CORRFLAG'=1),
    METHOD = 'st_hyper',
    CPP_CONTROL = cpp_ctrl,
    INIT = NULL,
    ITERATIONS_SUBSET = my_seq
)
mean((st_1$theta-true_theta)^2)
mean((PL_BFGS$theta-true_theta)^2)

obj$fit$post_index

cpp_ctrl <- list(
    MAXT = K*20,
    BURN = 500,
    ETA = 7,
    PAIRS_PER_ITERATION = 100,
    PROB = 5/K
)
my_seq <- seq(0, cpp_ctrl$MAXT, 100)
st_100 <- fit_plCFA(
    DATA_LIST = list('DATA' = manifest, 'CONSTRMAT' = constrMat, 'CORRFLAG'=1),
    METHOD = 'st_hyper',
    CPP_CONTROL = cpp_ctrl,
    INIT = NULL,
    ITERATIONS_SUBSET = my_seq
)

cpp_ctrl <- list(
    MAXT = K*20,
    BURN = 500,
    ETA = 7,
    PAIRS_PER_ITERATION = 390,
    PROB = 5/K
)
my_seq <- seq(0, cpp_ctrl$MAXT, 100)
st_390 <- fit_plCFA(
    DATA_LIST = list('DATA' = manifest, 'CONSTRMAT' = constrMat, 'CORRFLAG'=1),
    METHOD = 'st_hyper',
    CPP_CONTROL = cpp_ctrl,
    INIT = NULL,
    ITERATIONS_SUBSET = my_seq
)

library(tidyverse)
par_tab <- get_tidy_path(st_1, 'path_av_theta') %>%
    mutate(mod = 'nu1') %>%
    bind_rows(
        get_tidy_path(st_100, 'path_av_theta') %>%
            mutate(mod = 'nu100')
    ) %>%
    bind_rows(
        get_tidy_path(st_390, 'path_av_theta') %>%
            mutate(mod = 'nu390')
    ) %>%
    mutate(
        thresholds = map_dbl(path_av_theta, ~mean((.x[1:length(true_tau)]-true_theta[1:length(true_tau)])^2)),
        loadings = map_dbl(path_av_theta, ~mean((get_lambda(.x, categories, constrMat)-get_lambda(true_theta, categories, constrMat))^2)),
        corr = map_dbl(path_av_theta, ~mean((get_corr(.x, constrMat)-get_corr(true_theta, constrMat))^2))
    ) %>%
    gather(key = 'par', value = 'mse', thresholds, loadings, corr)

bfgs_tab <- tibble(
    thresholds = mean((PL_BFGS$theta[1:length(true_tau)]-true_theta[1:length(true_tau)])^2),
    loadings = mean((get_lambda(PL_BFGS$theta, categories, constrMat)-get_lambda(true_theta, categories, constrMat))^2),
    corr = mean((get_corr(PL_BFGS$theta, constrMat)-get_corr(true_theta, constrMat))^2)
) %>%
    gather(key = 'par', value = 'mse', thresholds, loadings, corr)


gg <- ggplot()+
    geom_line(data = par_tab, aes(x = iter, y = mse, col = mod)) +
    geom_hline(data = bfgs_tab, aes(yintercept = mse), linetype = 'dashed')+
    facet_wrap(vars(par), scales = 'free')+
    theme_minimal()
plotly::ggplotly(gg, dynamicTicks = T)
get_Sigma_u(constrMat, obj$theta)
get_Sigma_u(constrMat, true_latent)

eig_list <- map( split(obj$fit$path_av_theta, seq(nrow(obj$fit$path_av_theta))) , ~eigen(get_Sigma_u(constrMat, .x))$values)
sum(reduce(map(eig_list, ~sum(.x<0)), c) != 0)
eig_list %>% filter(`1`<0)
matrixcalc::is.positive.definite(get_Sigma_u(constrMat, obj$theta))
PL_BFGS$fit
