.analyseReps_Gauss <- function(geneNames, output.folder)
{
  chain.names <- c("Chain 1", "Chain 2")
  file.name.1 <- paste(output.folder, "/chain1", sep="")
  file.name.2 <- paste(output.folder, "/chain2", sep="")

  ## Read data
  lambda.mcmc.1     <- .readFast(paste(file.name.1, "/Lambda_mcmc", sep=""))
  mu.1              <- .readFast(paste(file.name.1, "/Mu_mcmc", sep=""))
  ro.1              <- .readFast(paste(file.name.1, "/Rho_mcmc", sep=""))
  gamma.1           <- .readFast(paste(file.name.1, "/Gamma_mcmc", sep="")) 
  coeff.1           <- .readFast(paste(file.name.1, "/B_mcmc", sep="")) 
  Y.1               <- read.table(paste(file.name.1, "/Y_mean", sep=""))
  lambda.exp.mcmc.1 <- .readFast(paste(file.name.1, "/Lambda_exp_mcmc", sep=""))
  X.1               <- as.vector(as.matrix(read.table(paste(file.name.1, "/DataMean_standarised", sep=""))))

  lambda.mcmc.2     <- .readFast(paste(file.name.2, "/Lambda_mcmc", sep=""))
  mu.2              <- .readFast(paste(file.name.2, "/Mu_mcmc", sep=""))
  ro.2              <- .readFast(paste(file.name.2, "/Rho_mcmc", sep=""))
  gamma.2           <- .readFast(paste(file.name.2, "/Gamma_mcmc", sep="")) 
  coeff.2           <- .readFast(paste(file.name.2, "/B_mcmc", sep="")) 
  Y.2               <- read.table(paste(file.name.2, "/Y_mean", sep=""))
  lambda.exp.mcmc.2 <- .readFast(paste(file.name.2, "/Lambda_exp_mcmc", sep=""))

  genes <- dim(lambda.mcmc.1)[2]
  pdf(paste(output.folder, "/ConvergencePlots.pdf", sep=""), 9, 9, title = "ConvergencePlots.pdf")
    old.par <- par(no.readonly = TRUE)
    par(mfrow = c(3,3))
    plot(colMeans(coeff.1), colMeans(coeff.2), xlab = chain.names[1], ylab = chain.names[2], main = "B") #, col = col.i)
    lines(c(-100,100), c(-100,100), col= "red")

    plot(colMeans(gamma.1), colMeans(gamma.2), xlab = chain.names[1], ylab = chain.names[2], main = "Gammas")
    lines(c(-100,100), c(-100,100), col= "red")

    plot(colMeans(ro.1), colMeans(ro.2), xlab = chain.names[1], ylab = chain.names[2], main = "Rho") #, col = col.i)
    lines(c(-100,100), c(-100,100), col= "red")

    plot(colMeans( lambda.mcmc.1), colMeans( lambda.mcmc.2), xlab = chain.names[1], ylab = chain.names[2], main = "Lambda")#, col = col.i)
    lines(c(-100,100000), c(-100,100000), col= "red")

    plot(colMeans( mu.1), colMeans( mu.2), xlab = chain.names[1], ylab = chain.names[2], main = "Mu")#, col = col.i)
    lines(c(-100,1000), c(-100,1000), col= "red")

    plot(as.vector(as.matrix(Y.1)), as.vector(as.matrix(Y.2)), xlab = chain.names[1], ylab = chain.names[2], main = "Y_means")#, col = col.i)
    lines(c(-100,100000), c(-100,100000), col= "red")

    plot(colMeans( lambda.exp.mcmc.1), colMeans( lambda.exp.mcmc.2), xlab = chain.names[1], ylab = chain.names[2], main = "Lambda Exp")#, col = col.i)
    lines(c(-100,100000), c(-100,100000), col= "red")
    par(old.par)
  dev.off()

  mean.B             <- matrix(colMeans(coeff.1), genes, genes)
  rownames(mean.B )  <- colnames(mean.B) <- geneNames  

  net.prob           <- matrix(colMeans(gamma.1), genes, genes)
  rownames(net.prob) <- colnames(net.prob) <- geneNames  
  diag(net.prob)     <- 0


  pdf(paste(output.folder, "/AnalysisPlots.pdf", sep=""), 9, 9, title = "AnalysisPlots.pdf")
    ## Network Heatmap
    .heatMap.ggplot(net.prob)
    ## Marginal uncertainty plot
    .marginalUncertaintyPlot(gamma.1, geneNames)
    ## Plot links per cutoff
    .plotCutOffGammas(colMeans(gamma.1),  main.text = "Number of links included in model vs Threshold used")
    ## Plot inferred data vs mean replicates
    plot(as.vector(as.matrix(Y.1)), X.1, xlab = "Mean of Reps", ylab = "Inferred Value", main = "Mean of replicates Data vs Inferred Value of Data")
    lines(c(-100,1000), c(-100,1000), col= "red")
  dev.off()

  ## Check link marginals have "reasonably" converged 
  .linkConvergenceMessage(cbind(colMeans(gamma.1), colMeans(gamma.2)))

  netLink  <- which(net.prob > -1, T)
  cytoNet  <- data.frame(geneNames[netLink[,2]], 
	      geneNames[netLink[,1]], net.prob[netLink], mean.B[netLink])
  colnames(cytoNet) <- c("From", "To", "Probability", "Strength")
  write.table(cytoNet, paste(output.folder, "/NetworkProbability_List.txt", sep=""), row.names=F, quote=F, sep="\t")

  ## Output probabilities in matrix format
  write.table(net.prob, paste(output.folder, "/NetworkProbability_Matrix.txt", sep=""))

}