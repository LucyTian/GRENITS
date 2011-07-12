NonLinearNet <- function(resultsFolder, timeSeries, ParamVec = NULL, chains = 2, user.seeds = NULL, 
			  Regulators = NULL, fixMe = NULL)
{
  # If no parameters given use default
  if (is.null(ParamVec)){ParamVec <- mcmc.defaultParams_nonLinear()}
  # Check to see if type of ParamVec correct
  param.type <- .paramVecType(ParamVec)
  if(param.type != "nonLinear_Param"){stop(paste("Parameter vector should be of type \"nonLinear_Param\", not \"", param.type, "\"", sep =""))}

  timeSeries <- as.matrix(timeSeries)
  resultsFolder <- .remove_trailingSlash(resultsFolder)    
  message(paste("Created folder", resultsFolder))
  dir.create(resultsFolder)

  # Write parameters to file
  ParamVec.ForFile <- c("NonLinearNet", ParamVec)
  names(ParamVec.ForFile)[1] <- "runType"
  # Write parameters to file
  write.table(ParamVec.ForFile, paste(resultsFolder, "/runInfo.txt", sep=""))
  # Write gene names to file
  .writeGeneNames(timeSeries, resultsFolder)

  # If no fixed values given, fix diag on
  fixMe <- .makeFixedGammasMat(dim(timeSeries)[1], Regulators, fixMe)

  if (!is.null(user.seeds)){chains <- length(user.seeds)}
  for (i in 1:chains)
  {
    if (!is.null(user.seeds)){set.seed(user.seeds[i])}
    resultsFolder.i <- paste(resultsFolder, "/chain", i, "/", sep = "")
    dir.create(resultsFolder.i)
    message(paste("Started MCMC chain", i, " ============= "))
    .Call("callSplines", timeSeries, resultsFolder.i, ParamVec, fixMe, PACKAGE = "GRENITS")
    message(paste("MCMC chain", i, "finished!"))
    write.table(fixMe, paste(resultsFolder.i, "/FixedGammaFile", sep=""), col.names = F, row.names = F)
  }  
}

