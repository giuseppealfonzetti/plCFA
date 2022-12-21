
#'@export
diagnosis_pNorm <- function(MOD_OBJ, LABEL = 'path_av_theta'){
    path <- MOD_OBJ$fit[[LABEL]][-nrow(MOD_OBJ$fit[[LABEL]]),]
    shift <- MOD_OBJ$fit[[LABEL]][-1,]

    normFun <- function(x){sqrt(crossprod(x))}
    diff <- shift - path
    check <- apply(diff, MARGIN = 1, FUN = normFun)/apply(path, MARGIN = 1, FUN = normFun)

    return(check)
}
