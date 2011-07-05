.analyseNonLinear <- function(geneNames, output.folder, data.matrix)
{
  chain.names <- c("Chain 1", "Chain 2")
  file.name.1 <- paste(output.folder, "/chain1", sep="")
  file.name.2 <- paste(output.folder, "/chain2", sep="")

  ## Read data
  lambda.mcmc.1 <- .readFast(paste(file.name.1, "/Lambda_mcmc", sep=""))
  mu.1          <- .readFast(paste(file.name.1, "/Mu_mcmc", sep=""))
  ro.1          <- .readFast(paste(file.name.1, "/Rho_mcmc", sep=""))
  gamma.1       <- .readFast(paste(file.name.1, "/Gamma_mcmc", sep="")) 
  tau.1         <- .readFast(paste(file.name.1, "/Tau_mcmc", sep="")) 
  all.f.1       <- read.table(paste(file.name.1, "/all_f", sep=""))      
  all.fsqr.1    <- read.table(paste(file.name.1, "/all_f_sqr", sep=""))      
  Full_F_sqr.1  <- read.table(paste(file.name.1, "/Full_F_sqr", sep=""))      

  lambda.mcmc.2 <- .readFast(paste(file.name.2, "/Lambda_mcmc", sep=""))
  mu.2          <- .readFast(paste(file.name.2, "/Mu_mcmc", sep=""))
  ro.2          <- .readFast(paste(file.name.2, "/Rho_mcmc", sep=""))
  gamma.2       <- .readFast(paste(file.name.2, "/Gamma_mcmc", sep="")) 
  tau.2         <- .readFast(paste(file.name.2, "/Tau_mcmc", sep="")) 
  all.f.2       <- read.table(paste(file.name.2, "/all_f", sep=""))      


  genes <- dim(lambda.mcmc.1)[2]
  pdf(paste(output.folder, "/ConvergencePlots.pdf", sep=""), 9, 12, title = "ConvergencePlots.pdf")
    old.par <- par(no.readonly = TRUE)
    par(mfrow = c(3,2))
    plot(colMeans(gamma.1), colMeans(gamma.2), xlab = chain.names[1], ylab = chain.names[2], main = "Gammas")
    lines(c(-100,100), c(-100,100), col= "red")

    plot(colMeans(tau.1), colMeans(tau.2), xlab = chain.names[1], ylab = chain.names[2], main = "tau") #, col = col.i)
    lines(c(-100,100000), c(-100,100000), col= "red")

    plot(colMeans(ro.1), colMeans(ro.2), xlab = chain.names[1], ylab = chain.names[2], main = "Rho") #, col = col.i)
    lines(c(-100,100), c(-100,100), col= "red")

    plot(colMeans( lambda.mcmc.1), colMeans( lambda.mcmc.2), xlab = chain.names[1], ylab = chain.names[2], main = "Lambda")#, col = col.i)
    lines(c(-100,100000), c(-100,100000), col= "red")

    plot(colMeans( mu.1), colMeans( mu.2), xlab = chain.names[1], ylab = chain.names[2], main = "Mu")#, col = col.i)
    lines(c(-100,1000), c(-100,1000), col= "red")

    plot(as.vector(as.matrix(all.f.1)), as.vector(as.matrix(all.f.2)), xlab = chain.names[1], ylab = chain.names[2], main = "mean f(x_t)")#, col = col.i)
    lines(c(-100,1000), c(-100,1000), col= "red")
    par(old.par)
  dev.off()

#   mean.B             <- matrix(colMeans(coeff.1), genes, genes)
#   rownames(mean.B )  <- colnames(mean.B) <- geneNames  
  ## Get value of tau when on
  tau.on.1 <- .tauON(tau.1, as.matrix(gamma.1))

  net.prob           <- matrix(colMeans(gamma.1), genes, genes)
  rownames(net.prob) <- colnames(net.prob) <- geneNames  
  diag(net.prob)     <- 0

  ## Read run parameters, remove text and make numeric
  parameters.run <- read.table(paste(output.folder, "/runInfo.txt", sep=""), as.is = T)[,1]
  parameters.run <- as.numeric(parameters.run[2:length(parameters.run)])

  pdf(paste(output.folder, "/AnalysisPlots.pdf", sep=""), 9, 9, title = "AnalysisPlots.pdf")
    ## Network Heatmap
    .heatMap.ggplot(net.prob)
    ## Marginal uncertainty plot
    .marginalUncertaintyPlot(gamma.1, geneNames)
    ## Plot links per cutoff
    .plotCutOffGammas(colMeans(gamma.1),  main.text = "Number of links included in model vs Threshold used")
    ## Plot prior and posterior tau
    .plotDistribAndPriorTau.sepSelf(tau.on.1, colMeans(gamma.1), parameters.run)
  dev.off()

 ## library(GRENITS);  output.folder <- "ExampleNonLinearNet";   analyse.output(output.folder)

  pdf(paste(output.folder, "/InferredFunctionPlots.pdf", sep=""), 9, 12, title = "InferredFunctionPlots.pdf")
    ## Inferred functions
    .plotSplinesFunctions(all.f.1, all.fsqr.1, Full_F_sqr.1, colMeans(gamma.1),  
			  as.matrix(data.matrix), geneNames, colMeans(mu.1))
  dev.off()

  ## Check link marginals have "reasonably" converged 
  .linkConvergenceMessage(cbind(colMeans(gamma.1), colMeans(gamma.2)))

  netLink  <- which(net.prob > -1, T)
  cytoNet  <- data.frame(geneNames[netLink[,2]], 
	      geneNames[netLink[,1]], net.prob[netLink]) #, mean.B[netLink])
  colnames(cytoNet) <- c("From", "To", "Probability") #, "Strength")
  write.table(cytoNet, paste(output.folder, "/NetworkProbability_List.txt", sep=""), row.names=F, quote=F, sep="\t")

  ## Output probabilities in matrix format
  write.table(net.prob, paste(output.folder, "/NetworkProbability_Matrix.txt", sep=""))

}