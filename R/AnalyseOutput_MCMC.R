## Analyse mcmc output
analyse.output <- function(output.folder, timeSeries = NULL)
{
  output.folder <- .remove_trailingSlash(output.folder)
  runinfo   <- read.table(paste(output.folder, "/runInfo.txt", sep=""), as.is=T)[,1]
  geneNames <- read.table(paste(output.folder, "/geneNames.txt", sep=""), as.is=T)[,1]
  if (runinfo[1] == "LinearNet")
  {
    .analyseLinear(geneNames, output.folder)    
  }
  if (runinfo[1] == "NonLinearNet")
  {
    if(is.null(timeSeries))
    {
      message("For analysis of the NonLinear run, the data set used is needed.")
    }else{
      .analyseNonLinear(geneNames, output.folder, timeSeries)
    }
  }
  if (runinfo[1] == "ReplicatesNet_student")
  {
    .analyseReps_Student(geneNames, output.folder)
  }
  if (runinfo[1] == "ReplicatesNet_gauss")
  {
    .analyseReps_Gauss(geneNames, output.folder)
  }
  message(paste("Analysis finished. Output plots can be found in folder: \"", output.folder, "\"", sep = ""))
}
