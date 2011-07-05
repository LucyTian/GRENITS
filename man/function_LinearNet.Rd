\name{LinearNet}
\alias{LinearNet}
\title{Dynamic Bayesian Network Inference Using Linear Interactions}
\description{
Run Bayesian inference of linear interaction network. The function generates MCMC chains that can later be 
analysed.
}
\usage{
LinearNet(resultsFolder, timeSeries, ParamVec = NULL, chains = 2, user.seeds = NULL)
}
\arguments{
\item{resultsFolder}{Name of output folder. The folder will be created and the output of the run 
will be placed there.}
\item{timeSeries}{Data matrix containing gene expression time series. Where genes will be placed in rows 
and time points in columns. Gene names may be included as row names.}
\item{ParamVec}{A parameter vector created using "mcmc.defaultParams_Linear". If none is given, 
default parameters will be used. The vector contains parameters associated to the priors as well as MCMC
run length. (See mcmc.defaultParams_Linear)}
\item{chains}{Number of MCMC chains to run.}
\item{user.seeds}{An optional vector with seeds to use for MCMC chains.}
}
\value{
For each chain run, a folder (chain1, chain2, ...) will be created and the output of the MCMC run will
be placed there. The files will be B_mcmc (the coeffcients of the linear regression), Gamma_mcmc (the 
indicator variables of Gibbs variable selection), Lambda_mcmc (the precision of each regression),
Mu_mcmc (the intercept of each regression) and Rho_mcmc (the network connectivity parameter). 
}

\references{
Morrissey, E.R., Juarez, M.A., Denby, K.J. and Burroughs, N.J. 2010. On reverse engineering of gene interaction 
networks using time course data with repeated measurements.  Bioinformatics 2010; doi: 10.1093/bioinformatics/btq421

Morrissey, E.R., Juarez, M.A., Denby, K.J. and Burroughs, N.J. 2011.  Inferring the topology of a non-linear sparse 
gene regulatory network using fully Bayesian spline regression. (submitted)
}
\seealso{
\code{\link{mcmc.defaultParams_Linear}, \link{analyse.output}}.
}
\keyword{LinearNet}
\examples{
  # Load A. thaliana circadian clock ODE generated data
  data(Athaliana_ODE)
  # Folder where raw runs will be kept and analysed
  output.folder <- paste(tempdir(), "/Example_LinearNet",sep="")
  # Run network inference, place raw results in output.folder
  LinearNet(output.folder, Athaliana_ODE)
  # Analyse raw results, place analysis plots and files in output.folder
  analyse.output(output.folder)
}