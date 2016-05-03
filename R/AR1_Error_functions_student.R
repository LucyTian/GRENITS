ReplicatesNet_student <- function(resultsFolder, timeSeries, numReps, ParamVec = NULL, chains = 2, user.seeds = NULL, Regulators = NULL, fixMe = NULL)
{
  # If no parameters given use default
  if (is.null(ParamVec)){ParamVec <- mcmc.defaultParams_student()}
  # Check to see if type of ParamVec correct
  param.type <- .paramVecType(ParamVec)
  if(param.type != "RepsStudent_Param"){stop(paste("Parameter vector should be of type \"RepsStudent_Param\", not \"", param.type, "\"", sep =""))}

  timeSeries <- as.matrix(timeSeries)
  resultsFolder <- .remove_trailingSlash(resultsFolder)    
  message(paste("Created folder", resultsFolder))
  dir.create(resultsFolder)

  # Add name and reps
  ParamVec <- c(ParamVec, numReps)
  ParamVec.ForFile <- c("ReplicatesNet_student", ParamVec)
  names(ParamVec.ForFile)[c(1, length(ParamVec.ForFile))] <- c("runType","NumReplicates")
  # Write parameters to file
  write.table(ParamVec.ForFile, paste(resultsFolder, "/runInfo.txt", sep=""))
  # Write gene names to file
  .writeGeneNames(timeSeries, resultsFolder)

  # If no fixed values given, fix diag on
  fixMe <- .makeFixedGammasMat(dim(timeSeries)[1], Regulators, fixMe)

  if (!is.null(user.seeds)){chains <- length(user.seeds)}
  for (i in 1:chains)
  {
    ## Set seed if given
    if (!is.null(user.seeds)){set.seed(user.seeds[i])}
    resultsFolder.i <- paste(resultsFolder, "/chain", i, "/",sep ="")
    dir.create(resultsFolder.i)
    message(paste("Started MCMC chain", i, " ============= "))
    .Call("callStudent_Error", timeSeries, resultsFolder.i, ParamVec, fixMe,PACKAGE = "GRENITS")
    message(paste("MCMC chain", i, "finished!"))
    write.table(fixMe, paste(resultsFolder.i, "/FixedGammaFile", sep=""), col.names = F, row.names = F)
  }  
}

