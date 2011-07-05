\name{read.chain}
\alias{read.chain}
\title{Read MCMC Chains}
\description{
Read MCMC chains for further analysis.
}
\usage{
read.chain(output.folder, chainNumber)
}
\details{
Read chains produced by NonLinearNet, LinearNet, ReplicatesNet_student and ReplicatesNet_gauss 
for further analysis.
}
\arguments{
\item{output.folder}{Name of folder (including path) where chains are kept}
\item{chainNumber}{Which of the chains will be read}
}
\value{
Returns a list of vectors/matrices with the value of the variables at each MCMC iteration.   
}

\references{
Morrissey, E.R., Juarez, M.A., Denby, K.J. and Burroughs, N.J. 2010. On reverse engineering of gene interaction 
networks using time course data with repeated measurements.  Bioinformatics 2010; doi: 10.1093/bioinformatics/btq421

Morrissey, E.R., Juarez, M.A., Denby, K.J. and Burroughs, N.J. 2011.  Inferring the topology of a non-linear sparse 
gene regulatory network using fully Bayesian spline regression. (submitted)
}

\seealso{
\code{\link{NonLinearNet}}, \code{\link{LinearNet}}, 
	\code{\link{ReplicatesNet_student} }, \code{\link{ReplicatesNet_gauss} }.
}
\keyword{ReadChains}
\examples{
  #############################################
  ## Run inference using one chain
  #############################################
  # Load A. thaliana circadian clock ODE generated data
  data(Athaliana_ODE)
  # Folder where raw runs will be kept and analysed
  output.folder <- paste(tempdir(), "/Example_LinearNet", sep="")
  # Run network inference, place raw results in output.folder
  # Run just one chain for example purpose
  LinearNet(output.folder, Athaliana_ODE, chains = 1)

  ###########################
  ## Read chain 
  ###########################
  chain1 <- read.chain(output.folder, 1)
}