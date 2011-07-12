\name{analyse.output}
\alias{analyse.output}
\title{Analysis Plots}
\description{
Analyse output from network inference functions. Basic convergence and analysis plots.
}
\usage{
analyse.output(output.folder, timeSeries = NULL)
}
\details{
Read first two chains run and plot some basic convergence plots (ConvergencePlots.pdf), analysis
plots (AnalysisPlots.pdf), as well as inferred network probabilities in two formats 
(NetworkProbability_List.txt and NetworkProbability_Matrix.txt).
}
\arguments{
\item{output.folder}{Name of folder (including path) where chains are kept}
\item{timeSeries}{Only used by NonLinearModel analysis. Data matrix containing gene expression 
time series. Where genes will be placed in rows and time points in columns.}
}
\value{
The output of the analysis will be four files (five if nonLinearNet). The contents of the two plot files change
depending on the network inference function used. 
\item{ConvergencePlots.pdf}{Basic convergence plots. The posterior means of each variable are compared.}
\item{AnalysisPlots.pdf}{Heatmap plot of network link probabilities as well as marginal network uncertainty
plot. A plot of the number of links predicted by the model for a given probability threshold. For 
ReplicatesNet_student, the posterior distribution of the degrees of freedom
are also plotted. For NonLinearNet, the posterior of the smoothness parameter is plotted.}
\item{NetworkProbability_List.txt}{Posterior probabilities for each network connection in list format, 
including posterior interaction strength for linear models. Can be imported with network analysis software 
such as cytoscape.}
\item{NetworkProbability_Matrix.txt}{Posterior probabilities for each network connection in matrix format.}
\item{ProbNumParents.txt}{Posterior probabilities for number of regulators for each gene.}
\item{InferredFunctionPlots.pdf}{(Only for nonLinearNet) Posterior distribution of predicted functions. 
Data values are plotted as circles.}
}

\references{
Morrissey, E.R., Juarez, M.A., Denby, K.J. and Burroughs, N.J. 2010. On reverse engineering of gene interaction 
networks using time course data with repeated measurements.  Bioinformatics 2010; doi: 10.1093/bioinformatics/btq421

Morrissey, E.R., Juarez, M.A., Denby, K.J. and Burroughs, N.J. 2011 Inferring the time-invariant 
topology of a nonlinear sparse gene regulatory network using fully Bayesian spline autoregression
Biostatistics 2011; doi: 10.1093/biostatistics/kxr009
}

\seealso{
\code{\link{NonLinearNet}}, \code{\link{LinearNet}}, 
	\code{\link{ReplicatesNet_student} }, \code{\link{ReplicatesNet_gauss}}.
}
\keyword{AnalyseOutput}
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