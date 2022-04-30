#' simulate initial values for Rstan/stan 
#' 
#' @param param.names names of the parameters.
#' @param param.lower.
#' @param param.upper.
#' @param distribution.
#' @param number number of sets of start values.
#' @return a list of lists of stan start values.
#' @examples
#' sim.init(param.names=c("mu_beta0", "mu_beta1", "sigma_beta0", "sigma_beta1", "sigma_error"), param.lower=c(-0.1, -0.01, 0.001, 0.001, 0.001), param.upper=c(8., 0.0, 0.07, 0.02, 0.05), distribution=c("normal", "normal", "lognormal", "lognormal", "lognormal"), number=1)
#' 
#' 
#' @export
"sim.init" <- function(param.names, param.lower, param.upper, distribution = rep("uniform", 
    length = length(param.upper)), number, verbose=FALSE) {
###
###    simulate initial values for Stan (list of lists)
###          for given ranges and distribution (uniform is the default)
###
###
###        called from bayesprobplot, BatchALTBayesian, BayesianLFP
###            bayesian.groupm.Dest.Degrad, bayesian.single.RMD, ejeffreys.prior
###
###         note that sim.init02 samples from the specified prior distribution
###
    if (length(param.lower) != length(param.upper)) 
        stop("length param.lower must equal length param.upper")
    if (length(param.names) != length(param.upper)) 
        stop("length param.names must equal length param.upper")
    inits.list <- list()
###
### loop over the number of chains
###
    for (i in 1:number) {
        values <- list()
        for (j in 1:length(param.lower)) {
            if(distribution[j]!="uniform")param <- rangePrior(param.lower[j], param.upper[j], distribution[j])
            switch(distribution[j],
                   uniform = {
                value <- runif(1, param.lower[j], param.upper[j])
              },
                   lognormal = {
                value <- rlnorm(1, param[1], param[2])
            },
                   normal = {
                value <- rnorm(1, param[1], param[2])
            }, {
                stop("Distribution ", distribution[j], " not recognizded in sim.init")
            })
            values[[j]] <- value
        }
        names(values) <- param.names
        inits.list[[i]] <- values
    }
    if(verbose)print(inits.list)
    inits.list
}
