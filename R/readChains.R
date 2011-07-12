## Analyse mcmc output
read.chain <- function(output.folder, chainNumber)
{
  runinfo   <- read.table(paste(output.folder, "/runInfo.txt", sep=""), as.is=T)[,1]
  geneNames <- read.table(paste(output.folder, "/geneNames.txt", sep=""), as.is=T)[,1]
  if (runinfo[1] == "LinearNet")
  {
    chainOut <- .read.linear.chain(output.folder, chainNumber, geneNames)
  }
  if (runinfo[1] == "NonLinearNet")
  {
    chainOut <- .read.splines.chain(output.folder, chainNumber, geneNames)
  }
  if (runinfo[1] == "ReplicatesNet_student")
  {
    chainOut <- .read.reps_student.chain(output.folder, chainNumber, geneNames)
  }
  if (runinfo[1] == "ReplicatesNet_gauss")
  {
    chainOut <- .read.reps_gauss.chain(output.folder, chainNumber, geneNames)
  }  
  chainOut
}

## Function to build vector of all-interactions lables 
.lablesGammaChain <- function(geneNames)
{
  fullLabels <- NULL
  for (i in geneNames)
  {
    fullLabels_i <- paste(i, "regs", geneNames, sep ="-")
    fullLabels <- c(fullLabels, fullLabels_i)
  }
  fullLabels 
}


.read.linear.chain <- function(output.folder, chainNumber, geneNames)
{
  all.interactions.Labels <- .lablesGammaChain(geneNames)
  chain.name  <- paste("chain", chainNumber, sep="")
  file.name.1 <- paste(output.folder, "/", chain.name, sep="")
  ## Check if chain exists, if not exit
  chainExists <- .checkChainExists(output.folder, chain.name)
  if(!chainExists){stop("Chain does not exist")}
  ## Read data
  fixedGammaMat  <- as.matrix(read.table(paste(file.name.1, "/FixedGammaFile", sep="")))
  lambda.mcmc.1  <- .readFast(paste(file.name.1, "/Lambda_mcmc", sep=""))
  colnames(lambda.mcmc.1) <- geneNames
  mu.1           <- .readFast(paste(file.name.1, "/Mu_mcmc", sep=""))
  colnames(mu.1) <- geneNames
  ro.1           <- .readFast(paste(file.name.1, "/Rho_mcmc", sep=""))
  gamma.1        <- .readFast(paste(file.name.1, "/Gamma_mcmc", sep="")) 
  notfixed.indx     <- !is.finite(fixedGammaMat)
  colnames(gamma.1) <- all.interactions.Labels[notfixed.indx]
  coeff.1       <- .readFast(paste(file.name.1, "/B_mcmc", sep="")) 
  reg.indx          <- which(notfixed.indx | fixedGammaMat == 1)
  colnames(coeff.1) <- all.interactions.Labels[reg.indx]

  list(lambda = lambda.mcmc.1,    mu =  mu.1,   rho =     ro.1, 
	gamma =       gamma.1, coeff = coeff.1)
}

.read.splines.chain <- function(output.folder, chainNumber, geneNames)
{
  all.interactions.Labels <- .lablesGammaChain(geneNames)
  chain.name  <- paste("chain", chainNumber, sep="")
  file.name.1 <- paste(output.folder, "/", chain.name, sep="")
  ## Check if chain exists, if not exit
  chainExists <- .checkChainExists(output.folder, chain.name)
  if(!chainExists){stop("Chain does not exist")}
  ## Read data
  fixedGammaMat     <- as.matrix(read.table(paste(file.name.1, "/FixedGammaFile", sep="")))
  lambda.mcmc.1     <- .readFast(paste(file.name.1, "/Lambda_mcmc", sep=""))
  colnames(lambda.mcmc.1) <- geneNames
  mu.1              <- .readFast(paste(file.name.1, "/Mu_mcmc", sep=""))
  colnames(mu.1)    <- geneNames
  ro.1              <- .readFast(paste(file.name.1, "/Rho_mcmc", sep=""))
  gamma.1           <- .readFast(paste(file.name.1, "/Gamma_mcmc", sep="")) 
  notfixed.indx     <- !is.finite(fixedGammaMat)
  colnames(gamma.1) <- all.interactions.Labels[notfixed.indx]
  tau.1             <- .readFast(paste(file.name.1, "/Tau_mcmc", sep="")) 
  reg.indx          <- which(notfixed.indx | fixedGammaMat == 1)
  colnames(tau.1)   <- all.interactions.Labels[reg.indx]
  all.f.1           <- read.table(paste(file.name.1, "/all_f", sep=""))  
  all.interactions.Labels.trans <- as.vector(t(matrix(all.interactions.Labels,length(geneNames),length(geneNames))))
  colnames(all.f.1) <- all.interactions.Labels.trans
  all.fsqr.1        <- read.table(paste(file.name.1, "/all_f_sqr", sep=""))      
  colnames(all.fsqr.1) <- all.interactions.Labels.trans
  list(lambda    = lambda.mcmc.1,  mu =  mu.1,   rho =     ro.1, 
	gamma    =       gamma.1, tau = tau.1, all.f =  all.f.1, 
	all.fsqr =    all.fsqr.1)
}

.read.reps_gauss.chain <- function(output.folder, chainNumber, geneNames)
{
  all.interactions.Labels <- .lablesGammaChain(geneNames)
  chain.name  <- paste("chain", chainNumber, sep="")
  file.name.1 <- paste(output.folder, "/", chain.name, sep="")
  ## Check if chain exists, if not exit
  chainExists <- .checkChainExists(output.folder, chain.name)
  if(!chainExists){stop("Chain does not exist")}
  ## Read data
  fixedGammaMat     <- as.matrix(read.table(paste(file.name.1, "/FixedGammaFile", sep="")))
  lambda.mcmc.1     <- .readFast(paste(file.name.1, "/Lambda_mcmc", sep=""))
  colnames(lambda.mcmc.1) <- geneNames
  mu.1              <- .readFast(paste(file.name.1, "/Mu_mcmc", sep=""))
  colnames(mu.1)    <- geneNames
  ro.1              <- .readFast(paste(file.name.1, "/Rho_mcmc", sep=""))
  gamma.1           <- .readFast(paste(file.name.1, "/Gamma_mcmc", sep="")) 
  notfixed.indx     <- !is.finite(fixedGammaMat)
  colnames(gamma.1) <- all.interactions.Labels[notfixed.indx]

  coeff.1           <- .readFast(paste(file.name.1, "/B_mcmc", sep=""))
  reg.indx          <- which(notfixed.indx | fixedGammaMat == 1)
  colnames(coeff.1) <- all.interactions.Labels[reg.indx]
 
  Y.1               <- read.table(paste(file.name.1, "/Y_mean", sep=""))
  rownames(Y.1)     <- geneNames
  lambda.exp.mcmc.1 <- .readFast(paste(file.name.1, "/Lambda_exp_mcmc", sep=""))
  colnames(lambda.exp.mcmc.1) <- geneNames
#   X.1               <- as.vector(as.matrix(read.table(paste(file.name.1, "/DataMean_standarised", sep=""))))
  X.1               <- as.matrix(read.table(paste(file.name.1, "/DataMean_standarised", sep="")))
  rownames(X.1)     <- geneNames
  list(    lambda =     lambda.mcmc.1,     mu =    mu.1, rho = ro.1, 
	    gamma =           gamma.1,  coeff = coeff.1,   Y =  Y.1, 
       lambda.exp = lambda.exp.mcmc.1, X_mean =     X.1)
}


.read.reps_student.chain <- function(output.folder, chainNumber, geneNames)
{
  all.interactions.Labels <- .lablesGammaChain(geneNames)
  chain.name  <- paste("chain", chainNumber, sep="")
  file.name.1 <- paste(output.folder, "/", chain.name, sep="")
  ## Check if chain exists, if not exit
  chainExists <- .checkChainExists(output.folder, chain.name)
  if(!chainExists){stop("Chain does not exist")}
  ## Read data
  fixedGammaMat     <- as.matrix(read.table(paste(file.name.1, "/FixedGammaFile", sep="")))
  lambda.mcmc.1     <- .readFast(paste(file.name.1, "/Lambda_mcmc", sep=""))
  colnames(lambda.mcmc.1) <- geneNames
  mu.1              <- .readFast(paste(file.name.1, "/Mu_mcmc", sep=""))
  colnames(mu.1)    <- geneNames
  ro.1              <- .readFast(paste(file.name.1, "/Rho_mcmc", sep=""))
  gamma.1           <- .readFast(paste(file.name.1, "/Gamma_mcmc", sep="")) 
  notfixed.indx     <- !is.finite(fixedGammaMat)
  colnames(gamma.1) <- all.interactions.Labels[notfixed.indx]
  coeff.1           <- .readFast(paste(file.name.1, "/B_mcmc", sep="")) 
  reg.indx          <- which(notfixed.indx | fixedGammaMat == 1)
  colnames(coeff.1) <- all.interactions.Labels[reg.indx]
  Y.1               <- read.table(paste(file.name.1, "/Y_mean", sep=""))
  rownames(Y.1)     <- geneNames
  lambda.exp.mcmc.1 <- .readFast(paste(file.name.1, "/Lambda_exp_mcmc", sep=""))
  colnames(lambda.exp.mcmc.1) <- geneNames
  X.1               <- as.vector(as.matrix(read.table(paste(file.name.1, "/DataMean_standarised", sep=""))))
  rownames(X.1)     <- geneNames
  deg.1             <- .readFast(paste(file.name.1, "/DegFreedom_mcmc", sep=""))
  colnames(deg.1) <- geneNames
  acceptanceRatio.1 <- read.table(paste(file.name.1, "/acceptanceRatio", sep=""))[,1]
  list(    lambda =     lambda.mcmc.1,     mu =    mu.1,         rho =  ro.1, 
	    gamma =           gamma.1,  coeff = coeff.1,           Y =   Y.1, 
       lambda.exp = lambda.exp.mcmc.1, X_mean =     X.1,  degFreedom = deg.1,
  acceptanceRatio = acceptanceRatio.1)
}


.checkChainExists <- function(output.folder, chain.name)
{
   dir.out    <- dir(output.folder)
   indx.chain <- grep(chain.name, dir.out)   
   return(length(indx.chain) > 0)  
}