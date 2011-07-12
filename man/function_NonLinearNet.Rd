\name{NonLinearNet}
\alias{NonLinearNet}
\title{Dynamic Bayesian Network Inference Using Non-Linear Interactions}
\description{
Run Bayesian inference of non-linear interaction network. Non linear interactions are modelled using 
Penalised Splines. The function generates MCMC chains that can later be analysed.
}
\usage{
NonLinearNet(     resultsFolder, timeSeries, ParamVec = NULL, chains = 2, user.seeds = NULL,
              Regulators = NULL, fixMe = NULL) 
}
\arguments{
\item{resultsFolder}{Name of output folder. The folder will be created and the output of the run will be placed there.}
\item{timeSeries}{Data matrix containing gene expression time series. Where genes will be placed in rows 
and time points in columns. Gene names may be included as row names.}
\item{ParamVec}{A parameter vector created using "mcmc.defaultParams_nonLinear". If none is given, 
default parameters will be used. The vector contains parameters associated to the priors as well as MCMC
run length. (See mcmc.defaultParams_nonLinear)}
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
be placed there. The files will be Gamma_mcmc (the indicator variables of Gibbs variable selection), 
Lambda_mcmc (the precision of each regression), Mu_mcmc (the intercept of each regression), 
Rho_mcmc (the network connectivity parameter), Tau_mcmc (the "smoothness parameter"), all_f (posterior mean
of all functions), all_f_sqr (posterior mean of the square of all functions) and Full_F_sqr (posterior 
mean of the square of the sum of all functions, for each regression). 
For the files all_f and all_f_sqr functions are placed in column-wise order. The file is filled by 
placing all interactions for each regression one after another.
}

\references{
Morrissey, E.R., Juarez, M.A., Denby, K.J. and Burroughs, N.J. 2011 Inferring the time-invariant 
topology of a nonlinear sparse gene regulatory network using fully Bayesian spline autoregression
Biostatistics 2011; doi: 10.1093/biostatistics/kxr009
}
\seealso{
\code{\link{mcmc.defaultParams_nonLinear}, \link{analyse.output}}.
}
\keyword{NonLinearNet}
\examples{
  # Synthetic data
  data(Athaliana_ODE)
  # Reduced data set to 3 genes 20 TP for faster run
  Athaliana_ODE.reduced <- Athaliana_ODE[c(1,3,5),1:20]
  # Folder where raw runs will be kept and later analysed
  output.folder <- paste(tempdir(), "/ExampleNonLinearNet", sep = "")
  # Run network inference and place raw results in output.folder
  NonLinearNet(output.folder , Athaliana_ODE.reduced)
  # Analyse raw results, place analysis plots and files in output.folder
  analyse.output(output.folder, Athaliana_ODE.reduced)
}