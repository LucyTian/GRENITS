\name{ReplicatesNet_student}
\alias{ReplicatesNet_student}
\title{Dynamic Bayesian Network Inference Using Linear Interactions and Student Experimental Noise}
\description{
Run Bayesian inference of linear interaction network on data with replicates. The replicates are assumed
to follow a Student distribution. The function generates MCMC chains that can later be analysed.
}
\usage{
ReplicatesNet_student(resultsFolder, timeSeries, numReps, ParamVec = NULL, chains = 2, user.seeds = NULL)
}
\arguments{
\item{resultsFolder}{Name of output folder. The folder will be created and the output of the run 
will be placed there.}
\item{timeSeries}{Data matrix containing gene expression time series. Where genes will be placed in rows 
and time points in columns. Each times series must be placed one after another, so that the first n 
columns correspond to time series replicate one, the next n columns to time series replicate two, etc.
Gene names may be included as row names.}
\item{numReps}{Number of replicate time series included in timeSeries matrix.}
\item{ParamVec}{A parameter vector created using "mcmc.defaultParams_student". If none is given, 
default parameters will be used. The vector contains parameters associated to the priors as 
well as MCMC run length. (See mcmc.defaultParams_student)}
\item{chains}{Number of MCMC chains to run.}
\item{user.seeds}{An optional vector with seeds to use for MCMC chains.}
}
\details{
The order in which the replicates are placed do not affect the output. In other words swapping timepoint
2 replicate 1 and timepoint 2 replicate 2 makes no difference. 
For the cases where a measurement is not available, an NaN may be used.
}
\value{
For each chain run, a folder (chain1, chain2, ...) will be created and the output of the MCMC run will
be placed there. The files will be B_mcmc (the coeffcients of the linear regression), Gamma_mcmc (the 
indicator variables of Gibbs variable selection), Lambda_mcmc (the precision of each regression),
Mu_mcmc (the intercept of each regression), Rho_mcmc (the network connectivity parameter), 
DataMean_standarised (a times series with the mean of the replicates), Lambda_exp_mcmc (the precision
of the replicate noise), Y_mean (the mean posterior value of the inferred "true mRNA"), DegFreedom_mcmc
(degrees of freedom of the student distributed replicate noise) and acceptanceRatio (acceptance ratio of the 
Metropolis-Hastings degrees of freedom update).
}

\references{
Morrissey, E.R., Juarez, M.A., Denby, K.J. and Burroughs, N.J. 2010. On reverse engineering of gene interaction 
networks using time course data with repeated measurements.  Bioinformatics 2010; doi: 10.1093/bioinformatics/btq421
}
\seealso{
\code{\link{mcmc.defaultParams_student}}, \code{\link{analyse.output}}.
}
\keyword{ReplicatesNet}
\examples{
  # Synthetic data
  data(Athaliana_ODE_4NoiseReps)
  # Folder where raw runs will be kept and later analysed
  output.folder <- paste(tempdir(), "/Example_ReplicatesNet_Student",sep="")
  # Run network inference and place raw results in output.folder
  ReplicatesNet_student(output.folder, Athaliana_ODE_4NoiseReps, numReps = 4)
  # Analyse raw results, place analysis plots and files in output.folder
  analyse.output(output.folder)
}