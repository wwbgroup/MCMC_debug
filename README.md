# MCMC_debug
Error in function runMCMC
mh_out <- runMCMC(bayesianSetup = bayesSEIR, sampler = "Metropolis", settings = mh_settings)

BT runMCMC: trying to find optimal start and covariance valuesError in Matrix::nearPD(MASS::ginv(-hessian)) : 
  Matrix seems negative semi-definite
In addition: Warning message:
Start values outside prior range 
runMCMC terminated after 0.270000000018626seconds
