## plot priors
## Analyse mcmc output
plotPriors <- function(parameter.vec)
{
  old.par <- par(no.readonly = TRUE)
  param.type <- .paramVecType(parameter.vec)
  if(param.type == "Linear_Param")
  {
    .plotPriorsLinear(parameter.vec)
  }
  # Dodgy work around
  if(param.type == "nonLinear_Param")
  {
    .plotPriorsNonLinear(parameter.vec)
  }
  if(param.type == "RepsStudent_Param")
  {
    .plotPriorsErrorStudent(parameter.vec)
  }
  # Dodgy work around
  if(param.type == "RepsGauss_Param")
  {
    .plotPriorsErrorGauss(parameter.vec)
  }
  par(old.par)
}

.plotPriorsLinear <- function(parameter.vec)
{
  par(mfrow=c(3,2))
  # Plot sample, burnin and thin
  barplot(parameter.vec[1:3],  main = "Run parameters")
  # Plot rho prior
  x.vec <- 0.001*(0:1000)
  plot(y = dbeta(x.vec, parameter.vec[4], parameter.vec[5]), 
	x = x.vec, type = "l", col = "red",
	xlab = "Rho", ylab = "Probability density", main = "Prior distribution for \nRho")

  # Plot B prior
  x.vec <- seq(-3*parameter.vec[6], 3*parameter.vec[6], length.out = 1000)
  plot(y = dnorm(x.vec, 0, parameter.vec[6]), 
	x = x.vec, type = "l", col = "red",
	xlab = "B", ylab = "Probability density", main = "Prior distribution for \nCoefficients")

  # Plot lambda prior
  sd.aux <- sqrt(parameter.vec[7]/(parameter.vec[8])^2)
  x.vec <- seq(0, 5*sd.aux, length.out = 1000)
  plot(y = dgamma(x.vec, parameter.vec[7], parameter.vec[8]), 
	x = x.vec, type = "l", col = "red",
	xlab = "Lambda", ylab = "Probability density", main = "Prior distribution for \nPrecision")

  # Plot B prior
  x.vec <- seq(-3*parameter.vec[9], 3*parameter.vec[9], length.out = 1000)
  plot(y = dnorm(x.vec, 0, parameter.vec[9]), 
	x = x.vec, type = "l", col = "red",
	xlab = "Mu", ylab = "Probability density", main = "Prior distribution for \nIntercept")
}


.plotPriorsErrorGauss <- function(parameter.vec)
{
  par(mfrow=c(3,2))
  # Plot sample, burnin and thin
  barplot(parameter.vec[c(1:3)],  main = "Run parameters") #, 12
  # Plot rho prior
  x.vec <- 0.001*(0:1000)
  plot(y = dbeta(x.vec, parameter.vec[4], parameter.vec[5]), 
	x = x.vec, type = "l", col = "red",
	xlab = "Rho", ylab = "Probability density", main = "Prior distribution for \nRho")

  # Plot B prior
  x.vec <- seq(-3*parameter.vec[6], 3*parameter.vec[6], length.out = 1000)
  plot(y = dnorm(x.vec, 0, parameter.vec[6]), 
	x = x.vec, type = "l", col = "red",
	xlab = "B", ylab = "Probability density", main = "Prior distribution for \nCoefficients")

  # Plot lambda prior
  sd.aux <- sqrt(parameter.vec[7]/(parameter.vec[8])^2)
  x.vec <- seq(0, 5*sd.aux, length.out = 1000)
  plot(y = dgamma(x.vec, parameter.vec[7], parameter.vec[8]), 
	x = x.vec, type = "l", col = "red",
	xlab = "Lambda", ylab = "Probability density", main = "Prior distribution for \nPrecision")

  # Plot lambda prior
  sd.aux <- sqrt(parameter.vec[9]/(parameter.vec[10])^2)
  x.vec <- seq(0, 5*sd.aux, length.out = 1000)
  plot(y = dgamma(x.vec, parameter.vec[9], parameter.vec[10]), 
	x = x.vec, type = "l", col = "red",
	xlab = "LambdaExp", ylab = "Probability density", main = "Prior distribution for \nPrecision of Data Replicates")

  # Plot mu prior
  x.vec <- seq(-3*parameter.vec[11], 3*parameter.vec[11], length.out = 1000)
  plot(y = dnorm(x.vec, 0, parameter.vec[11]), 
	x = x.vec, type = "l", col = "red",
	xlab = "Mu", ylab = "Probability density", main = "Prior distribution for \nIntercept")
}

.plotPriorsErrorStudent <- function(parameter.vec)
{
  par(mfrow=c(3,3))
  # Plot sample, burnin and thin
  barplot(parameter.vec[c(1:3)],  main = "Run parameters") #, 14
  # Plot rho prior
  x.vec <- 0.001*(0:1000)
  plot(y = dbeta(x.vec, parameter.vec[4], parameter.vec[5]), 
	x = x.vec, type = "l", col = "red",
	xlab = "Rho", ylab = "Probability density", main = "Prior distribution for \nRho")

  # Plot B prior
  x.vec <- seq(-3*parameter.vec[6], 3*parameter.vec[6], length.out = 1000)
  plot(y = dnorm(x.vec, 0, parameter.vec[6]), 
	x = x.vec, type = "l", col = "red",
	xlab = "B", ylab = "Probability density", main = "Prior distribution for \nCoefficients")

  # Plot lambda prior
  sd.aux <- sqrt(parameter.vec[7]/(parameter.vec[8])^2)
  x.vec <- seq(0, 5*sd.aux, length.out = 1000)
  plot(y = dgamma(x.vec, parameter.vec[7], parameter.vec[8]), 
	x = x.vec, type = "l", col = "red",
	xlab = "Lambda", ylab = "Probability density", main = "Prior distribution for \nPrecision")

  # Plot lambda prior
  sd.aux <- sqrt(parameter.vec[9]/(parameter.vec[10])^2)
  x.vec <- seq(0, 5*sd.aux, length.out = 1000)
  plot(y = dgamma(x.vec, parameter.vec[9], parameter.vec[10]), 
	x = x.vec, type = "l", col = "red",
	xlab = "LambdaExp", ylab = "Probability density", main = "Prior distribution for \nPrecision of Data Replicates")

  # Plot lambda prior
  sd.aux <- sqrt(parameter.vec[11]/(parameter.vec[12])^2)
  x.vec <- seq(0, 5*sd.aux, length.out = 1000)
  plot(y = dgamma(x.vec, parameter.vec[11], parameter.vec[12]), 
	x = x.vec, type = "l", col = "red",
	xlab = "Deg Freedom", ylab = "Probability density", 
	main = "Prior distribution for \nDegrees of Freedom")


  # Plot mu prior
  x.vec <- seq(-3*parameter.vec[13], 3*parameter.vec[13], length.out = 1000)
  plot(y = dnorm(x.vec, 0, parameter.vec[13]), 
	x = x.vec, type = "l", col = "red",
	xlab = "Mu", ylab = "Probability density", main = "Prior distribution for \nIntercept")
  

}


.plotPriorsNonLinear <- function(parameter.vec)
{
  par(mfrow=c(3,2))
  # Plot sample, burnin and thin
  barplot(parameter.vec[1:3],  main = "Run parameters")

  # Plot rho prior
  x.vec <- 0.001*(0:1000)
  plot(y = dbeta(x.vec, parameter.vec[4], parameter.vec[5]), 
	x = x.vec, type = "l", col = "red",
	xlab = "Rho", ylab = "Probability density", main = "Prior distribution for \nRho")

  # Plot tau prior
  xy_vals <- .Fpareto(parameter.vec[12], parameter.vec[6], num.points=10000)
  plot(y = xy_vals$y, 
	x = xy_vals$x, type = "l", col = "red",
	xlab = "Tau", ylab = "Probability density", main = "Prior distribution for \nNonLinearity parameter")

  # Plot lambda prior
  sd.aux <- sqrt(parameter.vec[9]/(parameter.vec[10])^2)
  x.vec <- seq(0, 5*sd.aux, length.out = 1000)
  plot(y = dgamma(x.vec, parameter.vec[9], parameter.vec[10]), 
	x = x.vec, type = "l", col = "red",
	xlab = "Lambda", ylab = "Probability density", main = "Prior distribution \nPrecision")

  # Plot Mu prior
  x.vec <- seq(-3*parameter.vec[11], 3*parameter.vec[11], length.out = 1000)
  plot(y = dnorm(x.vec, 0, parameter.vec[9]), 
	x = x.vec, type = "l", col = "red",
	xlab = "Mu", ylab = "Probability density", main = "Prior distribution for \nIntercept")

  # Plot B1 and B2 prior
  sd.aux <- 1/sqrt(parameter.vec[7])
  x.vec <- seq(-3*sd.aux, 3*sd.aux, length.out = 1000)
  plot(y = dnorm(x.vec, 0, sd.aux), 
	x = x.vec, type = "l", col = "red",
	xlab = "B", ylab = "Probability density", main = "Prior distribution for \nTwo Spline Coefficients")

}
