\name{LinearNet}
\alias{LinearNet}
\title{Dynamic Bayesian Network Inference Using Linear Interactions}
\description{
Run Bayesian inference of linear interaction network. The function generates MCMC chains that can later be 
analysed.
}
\usage{
LinearNet(    resultsFolder,   timeSeries, ParamVec = NULL, chains = 2, user.seeds = NULL, 
          Regulators = NULL, fixMe = NULL)
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
\item{Regulators}{An optional vector with the indices of which genes are regulators. If
provided, all non-regulator genes will not be allowed to regulate.}
\item{fixMe}{An optional matrix of size genes x genes, where columns represent regulators and 
rows regulated genes. The matrix informs the model of network connections known to be present/absent.
For each position use either 0 (no regulation, fix off), 1 (known regulatory interaction, fix on) 
or NaN (no information, do not fix).}
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

Morrissey, E.R., Juarez, M.A., Denby, K.J. and Burroughs, N.J. 2011 Inferring the time-invariant 
topology of a nonlinear sparse gene regulatory network using fully Bayesian spline autoregression
Biostatistics 2011; doi: 10.1093/biostatistics/kxr009
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