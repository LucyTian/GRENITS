LinearNet <- function(resultsFolder, timeSeries, ParamVec = NULL, chains = 2, user.seeds = NULL)
{
  # If no parameters given use default
  if (is.null(ParamVec)){  ParamVec <- mcmc.defaultParams_Linear()}
  # Check to see if type of ParamVec correct
  param.type <- .paramVecType(ParamVec)
  if(param.type != "Linear_Param"){stop(paste("Parameter vector should be of type \"Linear_Param\", not \"", param.type, "\"", sep =""))}

  resultsFolder <- .remove_trailingSlash(resultsFolder)
  timeSeries <- as.matrix(timeSeries)
  message(paste("Created folder", resultsFolder))
  dir.create(resultsFolder)

  # Write parameters to file
  write.table(c("LinearNet", ParamVec), paste(resultsFolder, "/runInfo.txt", sep=""))
  # Write gene names to file
  .writeGeneNames(timeSeries, resultsFolder)

  if (!is.null(user.seeds)){chains <- length(user.seeds)}
  for (i in 1:chains)
  {
    if (!is.null(user.seeds)){set.seed(user.seeds[i])}
    resultsFolder.i <- paste(resultsFolder, "/chain", i, "/", sep = "")
    dir.create(resultsFolder.i)
    message(paste( "Started MCMC chain", i, " ============= "))
    .Call("callAR1", timeSeries, resultsFolder.i, ParamVec, PACKAGE = "GRENITS")
    message(paste("MCMC chain", i, "finished!"))
  }  
}

