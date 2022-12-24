
library(RcppClock)
seed <- 1
set.seed(seed)
p <- 40; q <- 4; n <- 10000
constrMat <- build_constrMat(p,q,'simple')
true_load <- gen_loadings(FIXED = NULL, CONSTRAINT_MAT = constrMat, SEED = seed)
true_tau <- c(-1.2, 0, 1.2)

mat <- matrix(runif(q*q, .3,.7),q,q)
mat[!lower.tri(mat)] <- 0
mat <- mat + t(mat)
diag(mat) <- 1
true_latent <- as.matrix(Matrix::nearPD(mat)$mat)
#true_latent <- clusterGeneration::rcorrmatrix(d = q)
#true_latent <- matrix(.6,q,q); diag(true_latent) <- 1
# true_latent <- matrix(c(1,0.5995658,0.4935413,
#   0.5995658,1,0.6684667,
#   0.4935413,0.6684667,1),q,q)
true_theta <- get_theta(rep(true_tau, p), true_load, true_latent, cat, constrMat, 0)
seed1 <- 1
manifest <- gen_URV_data(n, true_load, true_tau, true_latent,
                         SEED = seed1)

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
transformed_rhos_init = rep(0.55, q*(q-1)/2)
#get_Sigma_u(constrMat, transformed_rhos_init)
par_init <- c(lambda0_init, lambda_init, transformed_rhos_init)

PL_BFGS <- fit_plCFA(
    DATA_LIST = list('DATA' = manifest, 'CONSTRMAT' = constrMat, 'CORRFLAG'=1),
    METHOD = 'ucminf',
    CPP_CONTROL = cpp_ctrl,
    INIT = NULL
)
#PL_BFGS
# ############################
#
# par <- c(lambda0_init, lambda_init, transformed_rhos_init)
# freq_tab <- pairs_freq(manifest, categories)
# dim(freq_tab)
# freq_tab[,1:10]
#
# # function for nll
# pair_nll <- function(par_vec){
#     lambda0_ <- par_vec[1:length(lambda0_init)]
#     lambda_ <- par_vec[(length(lambda0_init)+1):(length(lambda0_init)+length(lambda_init))]
#     transformed_rhos_ <- par_vec[(length(lambda0_init)+1+length(lambda_init)):length(par_vec)]
#     mod <-multiThread_completePairwise(
#         Y = manifest,
#         C_VEC = categories,
#         A = CONSTRMAT,
#         FREQ = freq_tab,
#         TAU = lambda0_,
#         LAMBDA = lambda_,
#         TRANSFORMED_RHOS = transformed_rhos_,
#         CORRFLAG = 1,
#         GRFLAG = 0,
#         SILENTFLAG = 1
#     )
#     out <- mod$iter_nll
#     return(out)
# }
#
# pair_nll(par)
#
# # function for gradient
# pair_gr <- function(par_vec){
#     lambda0_ <- par_vec[1:length(lambda0_init)]
#     lambda_ <- par_vec[(length(lambda0_init)+1):(length(lambda0_init)+length(lambda_init))]
#     transformed_rhos_ <- par_vec[(length(lambda0_init)+1+length(lambda_init)):length(par_vec)]
#     mod <-multiThread_completePairwise(
#         Y = manifest,
#         C_VEC = categories,
#         A = CONSTRMAT,
#         FREQ = freq_tab,
#         TAU = lambda0_,
#         LAMBDA = lambda_,
#         TRANSFORMED_RHOS = transformed_rhos_,
#         CORRFLAG = 1,
#         GRFLAG = 1,
#         SILENTFLAG = 1
#     )
#
#     out <- mod$iter_ngradient
#     return(out)
# }
# pair_gr(par)
#
# pair_H <- function(par_vec){
#     lambda0_ <- par_vec[1:length(lambda0_init)]
#     lambda_ <- par_vec[(length(lambda0_init)+1):(length(lambda0_init)+length(lambda_init))]
#     transformed_rhos_ <- par_vec[(length(lambda0_init)+1+length(lambda_init)):length(par_vec)]
#     mod <- multiThread_completePairwise(
#         Y = manifest,
#         C_VEC = categories,
#         A = CONSTRMAT,
#         FREQ = freq_tab,
#         TAU = lambda0_,
#         LAMBDA = lambda_,
#         TRANSFORMED_RHOS = transformed_rhos_,
#         CORRFLAG = 1,
#         GRFLAG = 2,
#         SILENTFLAG = 1
#     )
#
#     out <- mod$H_approx
#     return(out)
# }
# pair_H(par)
#
# library(RcppClock)
# ###########
# mod <- plCFA(
#     DATA = manifest,
#     C_VEC = categories,
#     CONSTRMAT = constrMat,
#     TAU = lambda0_init,
#     LAMBDA = lambda_init,
#     TRANSFORMED_RHOS = transformed_rhos_init,
#     CORRFLAG = 1,
#     SILENTFLAG = 1,
#     METHODFLAG = 0,
#     PAIRS_PER_ITERATION = 1,
#     PROB = 10/100,
#     ETA = 7,
#     BURN = 200,
#     MAXT = 1000,
#     TOLGRAD = 5e-3,
#     TOLPAR = 1e-5,
#     TOLOBJ = 1e-3,
#     TOLCOUNT = 10,
#     SEED = 123
# )
# summary(clock, units = 's')
# mod
# mod$scale
# 1/(n*100)
# n_pairs <- (p*(p-1)/2)
# weights <- sampling_step(
#     FULL_POOL = 1:n_pairs,
#     METHODFLAG = 1,
#     PROB = .1,
#     PAIRS_PER_ITERATION = round(n_pairs/10,0),
#     N_ITEMS = p,
#     SEED = 123,
#     SILENTFLAG = 1,
#     ITER = 1
# )
# weights
# length(weights)
#
# stepsize_grid <- seq(3, 20, 1)
# nll_grid <- purrr::map_dbl(stepsize_grid, ~{
#     obj <- plCFA(
#         Y = manifest,
#         C_VEC = categories,
#         A = constrMat,
#         TAU = lambda0_init,
#         LAMBDA = lambda_init,
#         TRANSFORMED_RHOS = transformed_rhos_init,
#         CORRFLAG = 1,
#         SILENTFLAG = 1,
#         METHODFLAG = 0,
#         PAIRS_PER_ITERATION = 1,
#         PROB = 10/100,
#         ETA = .x,
#         BURN = 200,
#         MAXT = 100,
#         TOLGRAD = 5e-3,
#         TOLPAR = 1e-5,
#         TOLOBJ = 1e-3,
#         TOLCOUNT = 10,
#         SEED = 123
#     )
#
#     obj$path_nll[obj$last_iter]
# })
# nll_grid
# stepsize_grid[which.min(nll_grid)]

########


# cpp_ctrl <- list(
#     MAXT = K*20,
#     BURN = 1000,
#     PAIRS_PER_ITERATION = 1,
#     PROB = .1,
#     SEED = seed1
# )
# stepsize_tuning(
#     DATA_LIST = list('DATA' = manifest, 'CONSTRMAT' = constrMat, 'CORRFLAG'=1),
#     METHOD = 'st_hyper',
#     CPP_CONTROL = cpp_ctrl,
#     STEPSIZE_GRID = c(.5,.25,.1,.05,.025,.01, .005, .001),
#     INIT = NULL
# )
d <- length(true_theta)
K = p*(p-1)/2; K
K/8
nu <- 8;
sqrt(K/(nu*n))
1/sqrt(K/(nu*n))


# cpp_ctrl <- list(
#     MAXT = K*20,
#     BURN = K*5,
#     ETA = 3e-3,
#     PAIRS_PER_ITERATION = 1,
#     PROB = 5/K,
#     SEED = seed1
# )
# my_seq <- seq(0, cpp_ctrl$MAXT, 100)
# st_1 <- fit_plCFA(
#         DATA_LIST = list('DATA' = manifest, 'CONSTRMAT' = constrMat, 'CORRFLAG'=1),
#         METHOD = 'st_hyper',
#         CPP_CONTROL = cpp_ctrl,
#         INIT = NULL,
#         ITERATIONS_SUBSET = NULL
#     )
#
#
#
# mean((st_1$theta-true_theta)^2)
# mean((PL_BFGS$theta-true_theta)^2)
#
# cpp_ctrl <- list(
#     MAXT = K*20,
#     BURN = K*5,
#     ETA = 3e-3,
#     PAIRS_PER_ITERATION = 4,
#     PROB = 5/K,
#     SEED = seed1
# )
# my_seq <- seq(0, cpp_ctrl$MAXT, 100)
# st_4 <- fit_plCFA(
#     DATA_LIST = list('DATA' = manifest, 'CONSTRMAT' = constrMat, 'CORRFLAG'=1),
#     METHOD = 'st_hyper',
#     CPP_CONTROL = cpp_ctrl,
#     INIT = NULL,
#     ITERATIONS_SUBSET = my_seq
# )
#obj$fit$post_index
eta <- 3e-3
cpp_ctrl <- list(
    MAXT = K*20,
    BURN = K*5,
    ETA = eta,
    PAIRS_PER_ITERATION = 8,
    PROB = 5/K,
    SEED = seed1
)
my_seq <- seq(0, cpp_ctrl$MAXT, 100)
st_8 <- fit_plCFA(
    DATA_LIST = list('DATA' = manifest, 'CONSTRMAT' = constrMat, 'CORRFLAG'=1),
    METHOD = 'st_hyper',
    CPP_CONTROL = cpp_ctrl,
    INIT = NULL,
    ITERATIONS_SUBSET = NULL
)

cpp_ctrl <- list(
    MAXT = K*20,
    BURN = K*5,
    ETA = eta,
    PAIRS_PER_ITERATION = 16,
    PROB = 5/K,
    SEED = seed1
)
my_seq <- seq(0, cpp_ctrl$MAXT, 100)
st_16 <- fit_plCFA(
    DATA_LIST = list('DATA' = manifest, 'CONSTRMAT' = constrMat, 'CORRFLAG'=1),
    METHOD = 'st_hyper',
    CPP_CONTROL = cpp_ctrl,
    INIT = NULL,
    ITERATIONS_SUBSET = NULL
)



cpp_ctrl <- list(
    MAXT = K*20,
    BURN = K*5,
    ETA = eta,
    PAIRS_PER_ITERATION = 32,
    PROB = 5/K,
    SEED = seed1
)
my_seq <- seq(0, cpp_ctrl$MAXT, 100)
st_32 <- fit_plCFA(
    DATA_LIST = list('DATA' = manifest, 'CONSTRMAT' = constrMat, 'CORRFLAG'=1),
    METHOD = 'st_hyper',
    CPP_CONTROL = cpp_ctrl,
    INIT = NULL,
    ITERATIONS_SUBSET = NULL
)

# cpp_ctrl <- list(
#     MAXT = 2000,
#     BURN = 1000,
#     ETA = 3e-3,
#     PAIRS_PER_ITERATION = 256,
#     PROB = 5/K,
#     SEED = seed1
# )
# my_seq <- seq(0, cpp_ctrl$MAXT, 100)
# st_128 <- fit_plCFA(
#     DATA_LIST = list('DATA' = manifest, 'CONSTRMAT' = constrMat, 'CORRFLAG'=1),
#     METHOD = 'st_hyper',
#     CPP_CONTROL = cpp_ctrl,
#     INIT = NULL,
#     ITERATIONS_SUBSET = NULL
# )

# mean((st_1$theta-true_theta)^2)
# mean((st_4$theta-true_theta)^2)
mean((st_8$theta-true_theta)^2)
mean((st_16$theta-true_theta)^2)
mean((st_32$theta-true_theta)^2)
# mean((st_128$theta-true_theta)^2)
mean((PL_BFGS$theta-true_theta)^2)

# st_1$clock[1,2]*st_1$clock[1,6]
# st_4$clock[1,2]*st_4$clock[1,6]
st_8$clock[1,2]*st_8$clock[1,6]
st_16$clock[1,2]*st_16$clock[1,6]
st_32$clock[1,2]*st_32$clock[1,6]
# st_128$clock[1,2]*st_128$clock[1,6]
PL_BFGS$num_time

library(tidyverse)
lab <- 'path_av_theta'
par_tab <- get_tidy_path(st_8, lab) %>%
    mutate(mod = 'nu8') %>%
    # bind_rows(
    #     get_tidy_path(st_4, lab) %>%
    #         mutate(mod = 'nu4')
    # ) %>%
    # bind_rows(
    #     get_tidy_path(st_8, lab) %>%
    #         mutate(mod = 'nu8')
    # ) %>%
    bind_rows(
        get_tidy_path(st_16, lab) %>%
            mutate(mod = 'nu16')
    ) %>%
    bind_rows(
        get_tidy_path(st_32, lab) %>%
            mutate(mod = 'nu32')
    ) %>%
    # bind_rows(
    #     get_tidy_path(st_128, lab) %>%
    #         mutate(mod = 'nu64')
    # ) %>%
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
    geom_line(data = par_tab, aes(x = iter, y = log(mse), col = mod), alpha = .8, size = .5) +
    geom_hline(data = bfgs_tab, aes(yintercept = log(mse)), linetype = 'dashed')+
    facet_wrap(vars(par), scales = 'free')+
    theme_minimal()
plotly::ggplotly(gg, dynamicTicks = T)
# get_Sigma_u(constrMat, obj$theta)
# get_Sigma_u(constrMat, true_latent)

# eig_list <- map( split(obj$fit$path_av_theta, seq(nrow(obj$fit$path_av_theta))) , ~eigen(get_Sigma_u(constrMat, .x))$values)
# sum(reduce(map(eig_list, ~sum(.x<0)), c) != 0)
# eig_list %>% filter(`1`<0)
# matrixcalc::is.positive.definite(get_Sigma_u(constrMat, obj$theta))
# PL_BFGS$fit

############
library(tidyverse)
mvrnormR <- function(n, mu, sigma) {
    ncols <- ncol(sigma)
    mu <- rep(mu, each = n) ## not obliged to use a matrix (recycling)
    mu + matrix(rnorm(n * ncols), ncol = ncols) %*% chol(sigma)
}
sim_setting <- expand_grid(n = c(2000, 4000), p = p, sim_id = 1:2) %>%
    mutate(id = row_number()) %>%
    select(id, everything())
test <- pbapply::pblapply(purrr::transpose(.l = sim_setting), function(x){

    SAMPLE_SIZE <- x$n
    LOADINGS <- true_load
    THRESHOLDS <- true_tau
    LATENT_COV <- true_latent
    SEED <- x$sim_id

    p <- nrow(LOADINGS)
    error_variance <- diag(1, p, p) - diag(diag(LOADINGS%*%LATENT_COV%*%t(LOADINGS)),p,p)
    errors <- rmvn(SAMPLE_SIZE = SAMPLE_SIZE, VAR = error_variance)

    dim <- ncol(error_variance)
    #a <- matrix(rnorm(dim*SAMPLE_SIZE), dim, SAMPLE_SIZE)
    a <- matrix(rep(1,dim*SAMPLE_SIZE), dim, SAMPLE_SIZE)
    b <- t(chol(error_variance))
    VAR <- error_variance

    errors
    #matrixprod(t(b), a)
    #f <- crossprod(b, a)
    #mvrnormR(SAMPLE_SIZE, 0, VAR)
    #MASS::mvrnorm(n = SAMPLE_SIZE, mu = rep(0, ncol(VAR)), Sigma = VAR)
    #t(t(chol(VAR))%*%matrix(rnorm(dim*SAMPLE_SIZE), dim, SAMPLE_SIZE))
}, cl = 2)

test <- pbapply::pblapply(purrr::transpose(.l = sim_setting), function(x){

    data <- gen_URV_data(x$n, true_load, true_tau, true_latent, SEED = x$sim_id)

    suppressMessages(
        mod_obj <- fit_plCFA(
            DATA_LIST = list('DATA' = data, 'CONSTRMAT' = constrMat, 'CORRFLAG'=1),
            METHOD = 'ucminf',
            CPP_CONTROL = list(),
            INIT = NULL,
            ITERATIONS_SUBSET = NULL
        )
    )
    mod_obj
}, cl = 2)



##############
eta <- 7

K = p*(p-1)/2; K
cpp_ctrl <- list(
    MAXT = K*20,
    BURN = 500,
    ETA = eta,
    PAIRS_PER_ITERATION = 1,
    PROB = 5/K,
    SEED = seed1,
    CHECKCONVERGENCE = T
)
my_seq <- seq(0, cpp_ctrl$MAXT, 100)

st_fit <- fit_plCFA(
    DATA_LIST = list('DATA' = manifest, 'CONSTRMAT' = constrMat, 'CORRFLAG'=1),
    METHOD = 'st_hyper',
    CPP_CONTROL = cpp_ctrl,
    INIT = NULL,
    ITERATIONS_SUBSET = my_seq
)
path <- st_fit$fit$path_av_theta[-nrow(st_fit$fit$path_av_theta),]
path_shift <- st_fit$fit$path_av_theta[-1,]
diff <- (path_shift-path)
check <- apply(diff, MARGIN = 1, FUN = norm, type = "2")/apply(path, MARGIN = 1, FUN = norm, type = "2")
check

test <- diagnosis_pNorm(st_fit, LABEL = 'path_av_theta')
check_level <- 1e-6
gg_norm <- tibble(iter = st_fit$iterations_subset[-1], check = test) %>%
    ggplot(aes(x = iter, y = check))+
    geom_line()+
    geom_hline(yintercept = check_level, linetype = 'dashed')
plotly::ggplotly(gg_norm)

st_fit$iterations_subset[(length(st_fit$iterations_subset)-10):length(st_fit$iterations_subset)]
st_fit$fit$last_iter
st_fit$fit$path_av_theta
length(st_fit$fit$checkPar)
length(st_fit$iterations_subset)
length(test)


##### check freq ####à
seed <- 1
set.seed(seed)
p <- 10; q <- 2; n <- 100
constrMat <- build_constrMat(p,q,'simple')
true_load <- gen_loadings(FIXED = NULL, CONSTRAINT_MAT = constrMat, SEED = seed)
true_tau <- c(-1.2, 0, 1.2)

mat <- matrix(runif(q*q, 0.3,.7),q,q)
mat[!lower.tri(mat)] <- 0
mat <- mat + t(mat)
diag(mat) <- 1
true_latent <- mat
#true_latent <- matrix(.212,q,q); diag(true_latent) <- 1
# true_latent <- matrix(c(1,0.5995658,0.4935413,
#   0.5995658,1,0.6684667,
#   0.4935413,0.6684667,1),q,q)
true_theta <- get_theta(rep(true_tau, p), true_load, true_latent, cat, constrMat, 0)
seed1 <- 16
manifest <- gen_URV_data(n, true_load, true_tau, true_latent,
                         SEED = seed1)

categories <- apply(manifest, 2, max) + 1
freq_tab <- pairs_freq(manifest, categories)
freq_tab
freq_tab[,1:16]

manifest1 <- manifest
manifest1[1,1] <-NA
freq_tab1 <- pairs_freq(manifest1, categories)
freq_tab1[,1:16]
