.analyseNonLinear <- function(geneNames, output.folder, data.matrix)
{
  chain.names <- c("Chain 1", "Chain 2")
  file.name.1 <- paste(output.folder, "/chain1", sep="")
  file.name.2 <- paste(output.folder, "/chain2", sep="")

  ## Read data
  lambda.mean.1        <- .readLargeFileReturnMean_c(paste(file.name.1, "/Lambda_mcmc", sep=""))
  mu.mean.1            <- .readLargeFileReturnMean_c(paste(file.name.1, "/Mu_mcmc", sep=""))
  ro.mean.1            <- .readLargeFileReturnMean_c(paste(file.name.1, "/Rho_mcmc", sep=""))
  tau.mean.1           <- .readLargeFileReturnMean_c(paste(file.name.1, "/Tau_mcmc", sep=""))  
  all.f.1              <- read.table(paste(file.name.1, "/all_f", sep=""))      
  all.fsqr.1           <- read.table(paste(file.name.1, "/all_f_sqr", sep=""))      
  Full_F_sqr.1         <- read.table(paste(file.name.1, "/Full_F_sqr", sep=""))      
  # Fixed gamma files
  fixedGammaMat        <- as.matrix(read.table(paste(file.name.1, "/FixedGammaFile", sep="")))
  meanProbs_NumParents <- .readGammaFile_Return_MeanAndNumParents_c(paste(file.name.1, "/Gamma_mcmc", sep=""), fixedGammaMat)
  fixedGammaMat        <- as.vector(fixedGammaMat)
  meanProbs1           <- meanProbs_NumParents[[1]]
  numParents           <- meanProbs_NumParents[[2]]


  lambda.mean.2        <- .readLargeFileReturnMean_c(paste(file.name.2, "/Lambda_mcmc", sep=""))
  mu.mean.2            <- .readLargeFileReturnMean_c(paste(file.name.2, "/Mu_mcmc", sep=""))
  ro.mean.2            <- .readLargeFileReturnMean_c(paste(file.name.2, "/Rho_mcmc", sep=""))
  tau.mean.2           <- .readLargeFileReturnMean_c(paste(file.name.2, "/Tau_mcmc", sep=""))
  meanProbs2           <- .readLargeFileReturnMean_c(paste(file.name.2, "/Gamma_mcmc", sep=""))
  all.f.2              <- read.table(paste(file.name.2, "/all_f", sep=""))      

  genes <- length(lambda.mean.2)
  pdf(paste(output.folder, "/ConvergencePlots.pdf", sep=""), 9, 12, title = "ConvergencePlots.pdf")
    old.par <- par(no.readonly = TRUE)
    par(mfrow = c(3,2))
    plot(meanProbs1, meanProbs2, xlab = chain.names[1], ylab = chain.names[2], main = "Gammas")
    lines(c(-100,100), c(-100,100), col= "red")

    plot(tau.mean.1, tau.mean.2, xlab = chain.names[1], ylab = chain.names[2], main = "tau") #, col = col.i)
    lines(c(-100,100000), c(-100,100000), col= "red")

    plot(ro.mean.1, ro.mean.2, xlab = chain.names[1], ylab = chain.names[2], main = "Rho") #, col = col.i)
    lines(c(-100,100), c(-100,100), col= "red")

    plot(lambda.mean.1, lambda.mean.2, xlab = chain.names[1], ylab = chain.names[2], main = "Lambda")#, col = col.i)
    lines(c(-100,100000), c(-100,100000), col= "red")

    plot(mu.mean.1, mu.mean.2, xlab = chain.names[1], ylab = chain.names[2], main = "Mu")#, col = col.i)
    lines(c(-100,1000), c(-100,1000), col= "red")

    plot(as.vector(as.matrix(all.f.1)), as.vector(as.matrix(all.f.2)), xlab = chain.names[1], ylab = chain.names[2], main = "mean f(x_t)")#, col = col.i)
    lines(c(-100,1000), c(-100,1000), col= "red")
    par(old.par)
  dev.off()

#   ## Get value of tau when on
#   tau.on.1 <- .tauON(tau.1, as.matrix(gamma.1))
  ## Use fixed gamma info for gamma
  notfixed.indx                <- !is.finite(fixedGammaMat)
  fullMeanGamma                <- fixedGammaMat
  fullMeanGamma[notfixed.indx] <- meanProbs1

  net.prob           <- matrix(fullMeanGamma, genes, genes)
  rownames(net.prob) <- colnames(net.prob) <- geneNames  
  diag(net.prob)     <- 0

  probMat_withZeros                 <- net.prob
  probMat_withZeros[!notfixed.indx] <- 0

  ## Read run parameters, remove text and make numeric
  parameters.run <- read.table(paste(output.folder, "/runInfo.txt", sep=""), as.is = T)[,1]
  parameters.run <- as.numeric(parameters.run[2:length(parameters.run)])

  pdf(paste(output.folder, "/AnalysisPlots.pdf", sep=""), 9, 9, title = "AnalysisPlots.pdf")
    ## Network Heatmap
    .heatMap.ggplot(net.prob)
    ## Marginal uncertainty plot
    .plotCutOffGammas(as.vector(net.prob),  main.text = "Number of links included in model vs Threshold used")
    ## Plot links per cutoff
    .plotDistribParents.LargeMat(probMat_withZeros, numParents, geneNames)
    ## Plot prior and posterior tau
#     .plotDistribAndPriorTau.sepSelf(tau.on.1, colMeans(gamma.1), parameters.run)
  dev.off()

 ## library(GRENITS);  output.folder <- "ExampleNonLinearNet";   analyse.output(output.folder)

  pdf(paste(output.folder, "/InferredFunctionPlots.pdf", sep=""), 9, 12, title = "InferredFunctionPlots.pdf")
    ## Inferred functions
    .plotSplinesFunctions(all.f.1, all.fsqr.1, Full_F_sqr.1, as.vector(net.prob),  
			  as.matrix(data.matrix), geneNames, mu.mean.1)
  dev.off()

  rownames(numParents) <- geneNames
  colnames(numParents) <- paste(0:length(geneNames), "Regulators")
  ## Write to file num parents
  write.table(numParents, paste(output.folder, "/ProbNumParents.txt", sep=""))

  ## Check link marginals have "reasonably" converged 
  .linkConvergenceMessage(cbind(meanProbs1, meanProbs2))

  netLink  <- which(net.prob > -1, T)
  cytoNet  <- data.frame(geneNames[netLink[,2]], 
	      geneNames[netLink[,1]], net.prob[netLink]) #, mean.B[netLink])
  colnames(cytoNet) <- c("From", "To", "Probability") #, "Strength")
  write.table(cytoNet, paste(output.folder, "/NetworkProbability_List.txt", sep=""), row.names=F, quote=F, sep="\t")

  ## Output probabilities in matrix format
  write.table(net.prob, paste(output.folder, "/NetworkProbability_Matrix.txt", sep=""))
}
