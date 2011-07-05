\name{Athaliana_ODE_4NoiseReps}
\alias{Athaliana_ODE_4NoiseReps}
\docType{data}
\title{Gene expression time series generated with ODE model with added noise}
\description{
This data set was generated using a circadian clock regulatory network ODE model after which noise was added.
The data was simulated using COPASI, subsampled to produce hourly data and after this the log was caluclated. 
Further to this, four replicates were produced by adding student noise.
}
\usage{data(Athaliana_ODE_4NoiseReps)}
\format{A matrix containing genes in rows and (1h) time points in columns. The four replicates are appended (columns)}
\references{
Hoops, S., Sahle, S., Gauges, R., Lee, C., Pahle, J., Simus, N., Singhal, M., Xu, L., Mendes, P. and Kummer, U.
(2006) COPASI a COmplex PAthway SImulator. Bioinformatics, 22, 3067-3074.


Locke, J.C.W., Kozma-Bognar, L., Gould, P.D., Feher, B., Kevei, E., Nagy, F., Turner, M.S., Hall, A. and Millar,  A.J. (2006) 
Experimental validation of a predicted feedback loop in the multi-oscillator clock of Arabidopsis thaliana. Molecular Systems Biology
}
\keyword{datasets}
