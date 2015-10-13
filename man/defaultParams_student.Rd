\name{mcmc.defaultParams_student}
\alias{mcmc.defaultParams_student}
\title{Default Parameters for Linear Model with Student distributed replicates}
\description{
Create parameter vector with default parameters for ReplicatesNet_student function
}
\usage{
mcmc.defaultParams_student()
}
\details{
Use this function to generate a template parameter vector to use non-default parameters for 
the ReplicatesNet_student model.
}
\value{
Returns a single vector with the following elements (in this order):
\item{(1) samples}{Number of MCMC iterations to run}
\item{(2) burn.in}{Number of initial iterations to discard as burn in}
\item{(3) thin}{Subsampling frequency}
\item{(4) c}{Shape parameter 1 for Beta(c,d) prior on rho (connectivity parameter)}
\item{(5) d}{Shape parameter 2 for Beta(c,d) prior on rho (connectivity parameter)}
\item{(6) sigma.s}{Standard deviation parameter for N(0,sigma.s) prior on B (Regression coefficients)}
\item{(7) a}{Shape parameter for Gamma(a,b) prior on lambda (Regression precision)} 
\item{(8) b}{Rate parameter for Gamma(a,b) prior on lambda (Regression precision)} 
\item{(9) a_exp}{Shape parameter for Gamma(a_exp,b_exp) prior on tau (Replicates precision)} 
\item{(10) b_exp}{Rate parameter for Gamma(a_exp,b_exp) prior on tau (Replicates precision)} 
\item{(11) a_deg}{Shape parameter for Gamma(a_deg,b_deg) prior on nu (Student degrees of freedom)} 
\item{(12) b_deg}{Rate parameter for Gamma(a_deg,b_deg) prior on nu (Student degrees of freedom)} 
\item{(13) sigma.mu}{Standard deviation parameter for N(0,sigma.mu) prior on mu (Regression intercept)}
\item{(14) fix.y.iter}{Number of iterations for which sampled data Y is fixed}
}

\references{
Morrissey, E.R., Juarez, M.A., Denby, K.J. and Burroughs, N.J. 2010. On reverse engineering of gene interaction 
networks using time course data with repeated measurements.  Bioinformatics 2010; doi: 10.1093/bioinformatics/btq421
}

\seealso{
\code{\link{plotPriors}}, \code{\link{ReplicatesNet_student}}.
}
\keyword{ReplicatesNet_student}
\examples{
    # Get default parameters
    linearNet_Student.params <- mcmc.defaultParams_student()

    # Change run length
    linearNet_Student.params[1] <- 200000

    # Change prior regression precision 
    linearNet_Student.params[7] <- 0.001
    linearNet_Student.params[8] <- 0.001

    # Plot to visualise changes
    plotPriors(linearNet_Student.params)

    ## Use to run ReplicatesNet_student ...
}
