## Build default parameter vectors for "Student Error"
mcmc.defaultParams_student <- function()
{
  params        <- c(150000, 5000, 10, 0.5, 0.5, 2, 0.1, 0.1, 2, 0.01, 3.4, 0.15, 2, 6)
  names(params) <- c( "samples", "burn.in", "thin", "c", "d", "sigma.s", "a", "b", "a_exp", "b_exp", "a_deg", "b_deg", "sigma.mu", "fix.y.iter")
  params
}

## Build default parameter vectors for "Linear Model"
mcmc.defaultParams_Linear <- function()
{
  params        <- c(100000, 10000, 10, 0.5, 0.5, 2, 2, 0.01, 2)
  names(params) <- c( "samples", "burn.in", "thin", "c", "d", "sigma.s", "a", "b", "sigma.mu")
  params
}

## Build default parameter vectors for "Non Linear Model"
mcmc.defaultParams_nonLinear <- function()
{
  params        <- c(100000, 10000, 10, 0.5, 0.5, 10000, 0.25, 13, 2, 0.01, 2, 1.5)
  names(params) <- c( "samples", "burn.in", "thin", "c", "d", "trunc", "tau0",  "M", "a", "b", "sigma.mu", "a_pareto")
  params
}

## Build default parameter vectors for "Gaussian Error"
mcmc.defaultParams_gauss <- function()
{
  params        <- c(150000, 5000, 10, 0.5, 0.5, 2, 0.1, 0.1, 2, 0.01, 2, 6)
  names(params) <- c( "samples", "burn.in", "thin", "c", "d", "sigma.s", "a", "b", "a_exp", "b_exp", "sigma.mu", "fix.y.iter")
  params
}
