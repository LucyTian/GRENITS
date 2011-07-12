\name{mcmc.defaultParams_nonLinear}
\alias{mcmc.defaultParams_nonLinear}
\title{Default Parameters for non-Linear Model}
\description{
Create parameter vector with default parameters for NonLinearNet function
}
\usage{
mcmc.defaultParams_nonLinear()
}
\details{
Use this function to generate a template parameter vector to use non-default parameters for 
the NonLinearNet model.
}
\value{
Returns a single vector with the following elements (in this order):
\item{(1) samples}{Number of MCMC iterations to run.}
\item{(2) burn.in}{Number of initial iterations to discard as burn in.}
\item{(3) thin}{Subsampling frequency}
\item{(4) c}{Shape parameter 1 for Beta(c,d) prior on rho (connectivity parameter)}
\item{(5) d}{Shape parameter 2 for Beta(c,d) prior on rho (connectivity parameter)}
\item{(6) trunc}{Truncation parameter for InvertedPareto prior on tau (smoothness parameter)}
\item{(7) tau0}{Precision parameter for N(0, tau0^(-0.5)) prior on B (first two coefficients)}
\item{(8) M}{Numer of knots used for each spline function}
\item{(9) a}{Shape parameter for Gamma(a,b) prior on lambda (Regression precision)} 
\item{(10) b}{Rate parameter for Gamma(a,b) prior on lambda (Regression precision)} 
\item{(11) sigma.mu}{Standard deviation parameter for N(0,sigma.mu) prior on mu (Regression intercept)}
\item{(12) a_pareto}{Pareto parameter for InvertedPareto prior on tau (smoothness parameter)}
}

\references{
Morrissey, E.R., Juarez, M.A., Denby, K.J. and Burroughs, N.J. 2011 Inferring the time-invariant 
topology of a nonlinear sparse gene regulatory network using fully Bayesian spline autoregression
Biostatistics 2011; doi: 10.1093/biostatistics/kxr009
}

\seealso{
\code{\link{plotPriors}}, \code{\link{NonLinearNet}}.
}
\keyword{NonLinearNet}
\examples{
    # Get default parameters
    nonLinearNet.params <- mcmc.defaultParams_nonLinear()

    # Change run length
    nonLinearNet.params[1] <- 150000

    # Change prior on smoothness parameter
    nonLinearNet.params[6] <- 30000 # Change truncation 
    nonLinearNet.params[12] <- 3 # Concentrate more mass close to linear region

    # Plot to check changes
    plotPriors(nonLinearNet.params)

    ## Use to run LinearNet ...
}
