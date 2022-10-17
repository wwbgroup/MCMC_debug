####MCMC part####
####Including Prior/Sample/LL####
####By Suyi 2022/09/24####

library(BayesianTools)
####SEIRfitting

#Prior distribution
default_pars_density <- function(pars) {
  d_vec <- rep(NA, 4)
  ## PAAT1
  PAAT1 <- pars[1]
  d_vec[1] <- dnorm(PAAT1, PAAT_mean1, PAAT_sd1, log = T)
  ## PAAT2
  PAAT2 <- pars[2]
  d_vec[2] <- dnorm(PAAT2, PAAT_mean2, PAAT_sd2, log = T)
  ## PAAT3
  PAAT3 <- pars[3]
  d_vec[3] <- dnorm(PAAT3, PAAT_mean3, PAAT_sd3, log = T)
  ## PAAT4
  PAAT4 <- pars[4]
  d_vec[4] <- dnorm(PAAT4, PAAT_mean4, PAAT_sd4, log = T)
  ##
  return(sum(d_vec))
}

#Sample function
default_pars_sampler <- function(n = 1) {
  #n = 1
  s_vec <- matrix(NA, n, 4)
  ## PAAT1 
  PAAT1 <- rnorm(n, PAAT_mean1, PAAT_sd1)
  s_vec[, 1] <- PAAT1
  ## PAAT2
  PAAT2 <- rnorm(n, PAAT_mean2, PAAT_sd2)
  s_vec[, 2] <- PAAT2
  ## PAAT3
  PAAT3 <- rnorm(n, PAAT_mean3, PAAT_sd3)
  s_vec[, 3] <- PAAT3
  ## PAAT3
  PAAT4 <- rnorm(n, PAAT_mean4, PAAT_sd4)
  s_vec[, 4] <- PAAT4
  return(s_vec)
}

## wrapper for the analysis run
## Create BayesianSetup and settings, lower/upper for parameters: beta(1x10 vector), PAAT1, PAAT2, PAAT3, PAAT4
#' @param init_set_list   initial settings produced by Params.R
#' @param randomize_startValue    this function will randomly generate an initial condition if this argument is set to T; If you want to specify your own initial condition, set to F
#' @param startValue  If randomize_startValue is set to F, you can your own initial condition using this argument; If randomize_startValue is set to T, this argument will be ignored
#' @param output_ret  Whether to output parameter estimates output by MCMC
#' @param run_id this ID is meant to distinguish different runs (different parameters, random seeds, etc.). Run_ID will be included in the file names of all outputs
#' @param skip_MCMC This is meant for redrawing all results without rerunning MCMC
#' @param panel_B_R_ylim the y limit in panel B of the main result plot

randomize_startValue=T
startValue=NA
output_ret=T 
run_id=0
panel_B_R_ylim=4
plot_combined_fig=T
pars_density=default_pars_density
pars_sampler=default_pars_sampler
pars_name=c("PAAT1", "PAAT2", "PAAT3", "PAAT4")
calc_clearance=T
n_burn_in=0
n_iterations=3


onset_obs <- init_settings$daily_new_case[[1]]
init_states <- init_settings$init_states
n_pars = length(pars_name)
n_stage = length(init_settings$stage_intervals)

#try the model
# SEIRpred(pars = c(0.6, 0.7, 0.8, 0.6), init_settings = init_sets_list)[, "Onset_expect"]
################################################################################################################
## likelihood function
loglh_func <- function(pars){
  ypred <- as.matrix(SEIRpred(pars, init_settings = init_settings))
  ypred <- ypred[,c(102:111)]
  ypred_all <- as.numeric(apply(ypred, 1, sum))    ## sum Onset1 to Onset10 for different age groups
  # meant to suppress warnings when ypred is negative
  suppressWarnings(p <- dpois(onset_obs, ypred_all))
  
  if(any(p == 0) || any(is.nan(p))){
    logL <- -Inf
  }else{
    logL <- sum(log10(p))
  }
  return(logL)
}
## take a try
#loglh_func(pars = c(0.6, 0.7, 0.8, 0.6))

## Create BayesianSetup and settings, lower/upper for parameters: PAAT1, PAAT2, PAAT3, PAAT4

## take a try 
#pars_sampler(n = 1)
pars_prior <- createPrior(density = pars_density, sampler = pars_sampler, 
                          lower = init_settings$par_lower, upper = init_settings$par_upper)

  bayesSEIR <- createBayesianSetup(loglh_func, prior = pars_prior)
  startValue=pars_sampler()
  
  ## DRAM: Adaptive MCMC, prior optimization, delayed rejection
  mh_settings = list(startValue = startValue,
                     adapt = T, DRlevels = 2, iterations = n_iterations, thin = 10)
  mh_out <- runMCMC(bayesianSetup = bayesSEIR, sampler = "Metropolis", settings = mh_settings)
  #plot(mh_out)
  mcmc_pars_estimate <- getSample(mh_out, start = n_burn_in+2, thin = 1)  ## set start = 2002 as the burn in period
  mcmc_pars_estimate <- round(mcmc_pars_estimate, 3)
