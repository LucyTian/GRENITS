\name{mcmc.defaultParams_Linear}
\alias{mcmc.defaultParams_Linear}
\title{Default Parameters for Linear Model}
\description{
Create parameter vector with default parameters for LinearNet function
}
\usage{
mcmc.defaultParams_Linear()
}
\details{
Use this function to generate a template parameter vector to use non-default parameters for 
the LinearNet model.
}
\value{
Returns a single vector with the following elements (in this order):
\item{(1) samples}{Number of MCMC iterations to run.}
\item{(2) burn.in}{Number of initial iterations to discard as burn in.}
\item{(3) thin}{Subsampling frequency}
\item{(4) c}{Shape parameter 1 for Beta(c,d) prior on rho (connectivity parameter)}
\item{(5) d}{Shape parameter 2 for Beta(c,d) prior on rho (connectivity parameter)}
\item{(6) sigma.s}{Standard deviation parameter for N(0,sigma.s) prior on B (Regression coefficients)}
\item{(7) a}{Shape parameter for Gamma(a,b) prior on lambda (Regression precision)} 
\item{(8) b}{Rate parameter for Gamma(a,b) prior on lambda (Regression precision)} 
\item{(9) sigma.mu}{Standard deviation parameter for N(0,sigma.mu) prior on mu (Regression intercept)}
}

\references{
Morrissey, E.R., Juarez, M.A., Denby, K.J. and Burroughs, N.J. 2010. On reverse engineering of gene interaction 
networks using time course data with repeated measurements.  Bioinformatics 2010; doi: 10.1093/bioinformatics/btq421

Morrissey, E.R., Juarez, M.A., Denby, K.J. and Burroughs, N.J. 2011 Inferring the time-invariant 
topology of a nonlinear sparse gene regulatory network using fully Bayesian spline autoregression
Biostatistics 2011; doi: 10.1093/biostatistics/kxr009
}

\seealso{
\code{\link{plotPriors}}, \code{\link{LinearNet}}.
}
\keyword{LinearNet}
\examples{
    # Get default parameters
    linearNet.params <- mcmc.defaultParams_Linear()

    # Change run length
    linearNet.params[1] <- 150000

    # Change prior regression precision 
    linearNet.params[7] <- 0.001
    linearNet.params[8] <- 0.001

    # Plot to check changes
    plotPriors(linearNet.params)

    ## Use to run LinearNet ...
}
