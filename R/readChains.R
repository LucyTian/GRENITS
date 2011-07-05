## Analyse mcmc output
read.chain <- function(output.folder, chainNumber)
{
  runinfo   <- read.table(paste(output.folder, "/runInfo.txt", sep=""), as.is=T)[,1]
  if (runinfo[1] == "LinearNet")
  {
    chainOut <- .read.linear.chain(output.folder, chainNumber)
  }
  if (runinfo[1] == "NonLinearNet")
  {
    chainOut <- .read.splines.chain(output.folder, chainNumber)
  }
  if (runinfo[1] == "ReplicatesNet_student")
  {
    chainOut <- .read.reps_student.chain(output.folder, chainNumber)
  }
  if (runinfo[1] == "ReplicatesNet_gauss")
  {
    chainOut <- .read.reps_gauss.chain(output.folder, chainNumber)
  }  
  chainOut
}

.read.linear.chain <- function(output.folder, chainNumber)
{
  chain.name  <- paste("chain", chainNumber, sep="")
  file.name.1 <- paste(output.folder, "/", chain.name, sep="")
  ## Check if chain exists, if not exit
  chainExists <- .checkChainExists(output.folder, chain.name)
  if(!chainExists){stop("Chain does not exist")}
  ## Read data
  lambda.mcmc.1 <- .readFast(paste(file.name.1, "/Lambda_mcmc", sep=""))
  mu.1          <- .readFast(paste(file.name.1, "/Mu_mcmc", sep=""))
  ro.1          <- .readFast(paste(file.name.1, "/Rho_mcmc", sep=""))
  gamma.1       <- .readFast(paste(file.name.1, "/Gamma_mcmc", sep="")) 
  coeff.1       <- .readFast(paste(file.name.1, "/B_mcmc", sep="")) 
  list(lambda = lambda.mcmc.1,  mu =  mu.1,   rho =     ro.1, 
	gamma =       gamma.1, coeff = coeff.1)
}

.read.splines.chain <- function(output.folder, chainNumber)
{
  chain.name  <- paste("chain", chainNumber, sep="")
  file.name.1 <- paste(output.folder, "/", chain.name, sep="")
  ## Check if chain exists, if not exit
  chainExists <- .checkChainExists(output.folder, chain.name)
  if(!chainExists){stop("Chain does not exist")}
  ## Read data
  lambda.mcmc.1 <- .readFast(paste(file.name.1, "/Lambda_mcmc", sep=""))
  mu.1          <- .readFast(paste(file.name.1, "/Mu_mcmc", sep=""))
  ro.1          <- .readFast(paste(file.name.1, "/Rho_mcmc", sep=""))
  gamma.1       <- .readFast(paste(file.name.1, "/Gamma_mcmc", sep="")) 
  tau.1         <- .readFast(paste(file.name.1, "/Tau_mcmc", sep="")) 
  all.f.1       <- read.table(paste(file.name.1, "/all_f", sep=""))      
  all.fsqr.1    <- read.table(paste(file.name.1, "/all_f_sqr", sep=""))      
  list(lambda    = lambda.mcmc.1,  mu =  mu.1,   rho =     ro.1, 
	gamma    =       gamma.1, tau = tau.1, all.f =  all.f.1, 
	all.fsqr =    all.fsqr.1)
}

.read.reps_gauss.chain <- function(output.folder, chainNumber)
{
  chain.name  <- paste("chain", chainNumber, sep="")
  file.name.1 <- paste(output.folder, "/", chain.name, sep="")
  ## Check if chain exists, if not exit
  chainExists <- .checkChainExists(output.folder, chain.name)
  if(!chainExists){stop("Chain does not exist")}
  ## Read data
  lambda.mcmc.1     <- .readFast(paste(file.name.1, "/Lambda_mcmc", sep=""))
  mu.1              <- .readFast(paste(file.name.1, "/Mu_mcmc", sep=""))
  ro.1              <- .readFast(paste(file.name.1, "/Rho_mcmc", sep=""))
  gamma.1           <- .readFast(paste(file.name.1, "/Gamma_mcmc", sep="")) 
  coeff.1           <- .readFast(paste(file.name.1, "/B_mcmc", sep="")) 
  Y.1               <- read.table(paste(file.name.1, "/Y_mean", sep=""))
  lambda.exp.mcmc.1 <- .readFast(paste(file.name.1, "/Lambda_exp_mcmc", sep=""))
  X.1               <- as.vector(as.matrix(read.table(paste(file.name.1, "/DataMean_standarised", sep=""))))
  list(    lambda =     lambda.mcmc.1,     mu =    mu.1, rho = ro.1, 
	    gamma =           gamma.1,  coeff = coeff.1,   Y =  Y.1, 
       lambda.exp = lambda.exp.mcmc.1, X_mean =     X.1)
}


.read.reps_student.chain <- function(output.folder, chainNumber)
{
  chain.name  <- paste("chain", chainNumber, sep="")
  file.name.1 <- paste(output.folder, "/", chain.name, sep="")
  ## Check if chain exists, if not exit
  chainExists <- .checkChainExists(output.folder, chain.name)
  if(!chainExists){stop("Chain does not exist")}
  ## Read data
  lambda.mcmc.1     <- .readFast(paste(file.name.1, "/Lambda_mcmc", sep=""))
  mu.1              <- .readFast(paste(file.name.1, "/Mu_mcmc", sep=""))
  ro.1              <- .readFast(paste(file.name.1, "/Rho_mcmc", sep=""))
  gamma.1           <- .readFast(paste(file.name.1, "/Gamma_mcmc", sep="")) 
  coeff.1           <- .readFast(paste(file.name.1, "/B_mcmc", sep="")) 
  Y.1               <- read.table(paste(file.name.1, "/Y_mean", sep=""))
  lambda.exp.mcmc.1 <- .readFast(paste(file.name.1, "/Lambda_exp_mcmc", sep=""))
  X.1               <- as.vector(as.matrix(read.table(paste(file.name.1, "/DataMean_standarised", sep=""))))
  deg.1             <- .readFast(paste(file.name.1, "/DegFreedom_mcmc", sep=""))
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