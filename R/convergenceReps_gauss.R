.analyseReps_Gauss <- function(geneNames, output.folder)
{
  chain.names <- c("Chain 1", "Chain 2")
  file.name.1 <- paste(output.folder, "/chain1", sep="")
  file.name.2 <- paste(output.folder, "/chain2", sep="")

  ## Read data
  lambda.mean.1        <- .readLargeFileReturnMean_c(paste(file.name.1, "/Lambda_mcmc", sep=""))
  mu.mean.1            <- .readLargeFileReturnMean_c(paste(file.name.1, "/Mu_mcmc", sep=""))
  ro.mean.1            <- .readLargeFileReturnMean_c(paste(file.name.1, "/Rho_mcmc", sep=""))
  meancoeff1           <- .readLargeFileReturnMean_c(paste(file.name.1, "/B_mcmc", sep=""))  
  fixedGammaMat        <- as.matrix(read.table(paste(file.name.1, "/FixedGammaFile", sep="")))
  meanProbs_NumParents <- .readGammaFile_Return_MeanAndNumParents_c(paste(file.name.1, "/Gamma_mcmc", sep=""), fixedGammaMat)
  fixedGammaMat        <- as.vector(fixedGammaMat)
  meanProbs1           <- meanProbs_NumParents[[1]]
  numParents           <- meanProbs_NumParents[[2]]

  Y.1               <- read.table(paste(file.name.1, "/Y_mean", sep=""))
  lambda.exp.mean.1 <- .readLargeFileReturnMean_c(paste(file.name.1, "/Lambda_exp_mcmc", sep=""))
  X.1               <- as.vector(as.matrix(read.table(paste(file.name.1, "/DataMean_standarised", sep=""))))


  lambda.mean.2     <- .readLargeFileReturnMean_c(paste(file.name.2, "/Lambda_mcmc", sep=""))
  mu.mean.2         <- .readLargeFileReturnMean_c(paste(file.name.2, "/Mu_mcmc", sep=""))
  ro.mean.2         <- .readLargeFileReturnMean_c(paste(file.name.2, "/Rho_mcmc", sep=""))
  meanProbs2        <- .readLargeFileReturnMean_c(paste(file.name.2, "/Gamma_mcmc", sep=""))
  meancoeff2        <- .readLargeFileReturnMean_c(paste(file.name.2, "/B_mcmc", sep=""))  
  Y.2               <- read.table(paste(file.name.2, "/Y_mean", sep=""))
  lambda.exp.mean.2 <- .readLargeFileReturnMean_c(paste(file.name.2, "/Lambda_exp_mcmc", sep=""))


  genes <- length(lambda.mean.1)
  pdf(paste(output.folder, "/ConvergencePlots.pdf", sep=""), 9, 9, title = "ConvergencePlots.pdf")
    old.par <- par(no.readonly = TRUE)
    par(mfrow = c(3,3))
    plot(meanProbs1, meanProbs2, xlab = chain.names[1], ylab = chain.names[2], main = "Gammas")
    lines(c(-100,100), c(-100,100), col= "red")

    plot(meancoeff1, meancoeff2, xlab = chain.names[1], ylab = chain.names[2], main = "B") #, col = col.i)
    lines(c(-100,100), c(-100,100), col= "red")

    plot(ro.mean.1, ro.mean.2, xlab = chain.names[1], ylab = chain.names[2], main = "Rho") #, col = col.i)
    lines(c(-100,100), c(-100,100), col= "red")

    plot(lambda.mean.1, lambda.mean.2, xlab = chain.names[1], ylab = chain.names[2], main = "Lambda")#, col = col.i)
    lines(c(-100,100000), c(-100,100000), col= "red")

    plot(mu.mean.1, mu.mean.2, xlab = chain.names[1], ylab = chain.names[2], main = "Mu")#, col = col.i)
    lines(c(-100,1000), c(-100,1000), col= "red")

    plot(as.vector(as.matrix(Y.1)), as.vector(as.matrix(Y.2)), xlab = chain.names[1], ylab = chain.names[2], main = "Y_means")#, col = col.i)
    lines(c(-100,100000), c(-100,100000), col= "red")

    plot(lambda.exp.mean.1, lambda.exp.mean.2, xlab = chain.names[1], ylab = chain.names[2], main = "Lambda Exp")#, col = col.i)
    lines(c(-100,100000), c(-100,100000), col= "red")
    par(old.par)
  dev.off()

  ## Use fixed gamma info for gamma
  notfixed.indx                <- !is.finite(fixedGammaMat)
  fullMeanGamma                <- fixedGammaMat
  fullMeanGamma[notfixed.indx] <- meanProbs1
  # make a copy with fixed elements set to zero

  ## Use fixed gamma info for B
  reg.indx             <- which(notfixed.indx | fixedGammaMat == 1)
  fullMeanB            <- fixedGammaMat
  fullMeanB[reg.indx ] <- meancoeff1

#   mean.B             <- matrix(meancoeff1, genes, genes)
  mean.B             <- matrix(fullMeanB, genes, genes)
  rownames(mean.B )  <- colnames(mean.B) <- geneNames  

#   net.prob           <- matrix(meanProbs1, genes, genes)
  net.prob           <- matrix(fullMeanGamma, genes, genes)
  rownames(net.prob) <- colnames(net.prob) <- geneNames  
  diag(net.prob)     <- 0

  probMat_withZeros                 <- net.prob
  probMat_withZeros[!notfixed.indx] <- 0

  pdf(paste(output.folder, "/AnalysisPlots.pdf", sep=""), 9, 9, title = "AnalysisPlots.pdf")
    ## Network Heatmap
    .heatMap.ggplot(net.prob)
    ## Plot links per cutoff
    .plotCutOffGammas(as.vector(net.prob),  main.text = "Number of links included in model vs Threshold used")
    ## Marginal uncertainty plot
    .plotDistribParents.LargeMat(probMat_withZeros, numParents, geneNames)
    ## Plot inferred data vs mean replicates
    plot(as.vector(as.matrix(Y.1)), X.1, xlab = "Mean of Reps", ylab = "Inferred Value", main = "Mean of replicates Data vs Inferred Value of Data")
    lines(c(-100,1000), c(-100,1000), col= "red")
  dev.off()

  rownames(numParents) <- geneNames
  colnames(numParents) <- paste(0:length(geneNames), "Regulators")
  ## Write to file num parents
  write.table(numParents, paste(output.folder, "/ProbNumParents.txt", sep=""))

  ## Check link marginals have "reasonably" converged 
  .linkConvergenceMessage(cbind(colMeans(meanProbs1), colMeans(meanProbs2)))

  netLink  <- which(net.prob > -1, T)
  cytoNet  <- data.frame(geneNames[netLink[,2]], 
	      geneNames[netLink[,1]], net.prob[netLink], mean.B[netLink])
  colnames(cytoNet) <- c("From", "To", "Probability", "Strength")
  write.table(cytoNet, paste(output.folder, "/NetworkProbability_List.txt", sep=""), row.names=F, quote=F, sep="\t")

  ## Output probabilities in matrix format
  write.table(net.prob, paste(output.folder, "/NetworkProbability_Matrix.txt", sep=""))

}
