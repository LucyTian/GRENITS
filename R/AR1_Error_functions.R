ReplicatesNet_gauss <- function(resultsFolder, timeSeries, numReps, ParamVec = NULL, chains = 2, user.seeds = NULL)
{
  # If no parameters given use default
  if (is.null(ParamVec)){ParamVec <- mcmc.defaultParams_gauss()}
  # Check to see if type of ParamVec correct
  param.type <- .paramVecType(ParamVec)
  if(param.type != "RepsGauss_Param"){stop(paste("Parameter vector should be of type \"RepsGauss_Param\", not \"", param.type, "\"", sep =""))}

  resultsFolder <- .remove_trailingSlash(resultsFolder)    
  timeSeries <- as.matrix(timeSeries)
  message(paste("Created folder", resultsFolder))
  dir.create(resultsFolder)

  # Add name and reps
  ParamVec <- c(ParamVec, numReps)
  ParamVec.ForFile <- c("ReplicatesNet_gauss", ParamVec)
  names(ParamVec.ForFile)[c(1, length(ParamVec.ForFile))] <- c("runType","NumReplicates")
  # Write parameters to file
  write.table(ParamVec.ForFile, paste(resultsFolder, "/runInfo.txt", sep=""))
  # Write gene names to file
  .writeGeneNames(timeSeries, resultsFolder)

  if (!is.null(user.seeds)){chains <- length(user.seeds)}
  for (i in 1:chains)
  {
    ## Set seed if given
    if (!is.null(user.seeds)){set.seed(user.seeds[i])}
    resultsFolder.i <- paste(resultsFolder, "/chain", i, "/",sep ="")
    dir.create(resultsFolder.i)
    message(paste("Started MCMC chain", i, " ============= "))
    .Call("callGauss_Error", timeSeries, resultsFolder.i, ParamVec, PACKAGE = "GRENITS")
    message(paste("MCMC chain", i, "finished!"))
  }
}

