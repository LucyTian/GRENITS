.getSquareMatIndicesFromFlat <- function(flatIndex, matSize)
{
  colIndex              <- ceiling(flatIndex/matSize)
  rowIndex              <- flatIndex%%matSize
  rowIndex[rowIndex==0] <- matSize
  cbind(rowIndex,colIndex)
}

.makeFixedGammasMat <- function(genes, Regulators.vec, fixMe)
{
  if (is.null(fixMe)){   
    fixMe       <- matrix(NaN, genes, genes)
    diag(fixMe) <- 1
  }

  if (!is.null(Regulators.vec)){
    # Make non regulators index vector
    nonRegs.indx         <- 1:genes
    nonRegs.indx         <- nonRegs.indx[-1*Regulators.vec]
    # Store self interaction element
    diag.aux             <- diag(fixMe)    
    # Fix non regulators regulation elems to off
    fixMe[,nonRegs.indx] <- 0
    # Restore self interaction to initial value
    diag(fixMe)          <- diag.aux
  } 
  fixMe
}

.readLargeFileReturnMean_c <- function(file.name)
{
  aa <- .Call("readLargeFileGetMean", file.name, PACKAGE = "GRENITS")
  aa
}

.readGammaFile_Return_MeanAndNumParents_c <- function(file.name, linksFixed)
{  
  aa <- .Call("readGamma_getMean_numParents", file.name, as.matrix(linksFixed), PACKAGE = "GRENITS")
  aa
}


.readLargeFileReturnMean <- function(Bfile, chain.len)
{
  con_B         <- file(Bfile, "r")
  B.j           <- scan(con_B, nlines  =1, sep=",", quiet = T) #, flush = T)
  sumB          <- B.j
  for(j in 2:chain.len)
  {    
#     B.j <- try(scan(con_B, nlines  =1, sep=",", quiet = T), silent = T)
    B.j  <- scan(con_B, nlines  =1, sep=",", quiet = T)
    sumB <- sumB + B.j
  }
  close(con_B)
  # Return mean
  sumB/chain.len
}


.paramVecType <- function(parameter.vec)
{
  length.priors <- length(parameter.vec)
  if(length.priors == 9)
  {
    param.type <- "Linear_Param"
  }
  # Dodgy work around
  if(names(parameter.vec)[6] == "trunc")
  {
    param.type <- "nonLinear_Param"
  }
  if(length.priors == 14)
  {
    param.type <- "RepsStudent_Param"
  }
  # Dodgy work around
  if((length.priors == 12) && (names(parameter.vec)[6] != "trunc"))
  {
    param.type <- "RepsGauss_Param"
  }
  param.type
}

.remove_trailingSlash <- function(folderName)
{
    len.name <- nchar(folderName)
    lastChar <- substr(folderName, len.name, len.name)
    # Remove last char
    if (lastChar == "/" | lastChar=="\\"){folderName <- substr(folderName, 1, len.name-1)}
    return(folderName)
}

.readFast <- function(file.name)
{
  numElems <- length(scan(file.name, sep =",", nlines=1, quiet = T))
  aa       <- matrix(scan(file.name, sep =",", quiet = T), ncol = numElems, byrow = T)
  if (sum(is.nan(aa)) !=0)
  {
    print(paste("Warning found NANs in ", file.name))
  }
  aa
}

.writeGeneNames <- function(dataMatrix, resultsFolder)
{
  numGenes  <- dim(dataMatrix)[1]
  geneNames <- rownames(dataMatrix)
  ## Check to see if they are numbers or if any names (matrix problem)  
  num.geneNames <- suppressWarnings(as.numeric(geneNames))
  if(!any(is.na(num.geneNames)) | (length(geneNames) != numGenes))
  {
    geneNames <- paste("G", 1:numGenes)
  }
  write.table(geneNames, paste(resultsFolder, "/geneNames.txt", sep=""))
}



## If any links are more than problem.thresh between chains report
.linkConvergenceMessage <- function(link.probs, problem.thresh = 0.2)
{
  link.diff <- abs(link.probs[,1] - link.probs[,2])
  if (sum(link.diff>problem.thresh) != 0)
  {
    sink("CONVERGENCE_WARNING")
    print(paste("There are", sum(link.diff>problem.thresh), "links with a differnce of more than", problem.thresh))
    print("Check ConvergencePlots.pdf and consider a longer run.")
    sink()
    message(paste("There are", sum(link.diff>problem.thresh), "link probabilities that have not fully converged."))
    message("Check ConvergencePlots.pdf and consider a longer run.")
  }
}


## -----------------------------------------------------------------------------------------------
##  Cut off gammas   ------------------------------------------------------------
## -----------------------------------------------------------------------------------------------
.plotCutOffGammas <- function(link.prob.1, chainName="", true.net = NULL, main.text = "", nonTF = 0)
{
  genes     <- sqrt(length(link.prob.1))
  if ( nonTF == 0 | genes == nonTF)
  {
    diag.elem <- diag(matrix(1:length(link.prob.1), genes, genes))
  }else{
    diag.elem <- matrix(1:length(link.prob.1), nrow = nonTF)[,1]
  }
  
  remove.me <- NULL
  ## If large, remove links (<0.01) 
  if(length(link.prob.1)> 1000)
  {	
    remove.me <- which(link.prob.1 < 0.01)
  }
  ## Remove diagonals
#   links <- sort(link.prob.1[-diag.elem], decreasing =T)
  link.noDiag <- link.prob.1[-c(diag.elem, remove.me)]
  new.order   <- order(link.noDiag , decreasing =T)
  links       <- link.noDiag[new.order]
  numLinks    <- length(links)
  col.all     <- rep("black", numLinks)
  pch.all     <- rep(1, numLinks)
  if (length(true.net)!=0)
  {
    true.vec     <- as.vector(true.net)[-diag.elem]
    ordered.true <- true.vec[new.order] != 0
    col.all[ordered.true] <- "red"
    pch.all[ordered.true] <- 4
    
  }
  plot(x = links, y = 1:length(links),
       xlab = "Link probability threshold",ylab = "Number of links above threshold", 
       xlim = c(0,1), col = col.all, pch = pch.all, main = main.text, cex = 2.5 )#, cex.lab = 1.7) 
#, main =paste("Number of links for given threshold", chainName))
  
}


## -----------------------------------------------------------------------------------------------
##  Heatmap using ggplot              ------------------------------------------------------------
## -----------------------------------------------------------------------------------------------
.heatMap.ggplot <- function(prob.mat, size.max = 80)
{
      title.string <- "Link Probability"
      ## If too large, sub select
      if(dim(prob.mat)[1]> size.max)
      {
	topNumLinks <- floor(size.max/2)
	prob.mat.largest <- order(prob.mat, decreasing= T)
	rowsAndCols      <- as.vector(.getSquareMatIndicesFromFlat(prob.mat.largest[1:topNumLinks], dim(prob.mat)[1]))
	prob.mat         <- prob.mat[rowsAndCols, rowsAndCols]
	title.string     <- paste(title.string, "(Selected genes)")
      }
      all.m <- melt(as.matrix(prob.mat))
      (p <- ggplot(all.m, aes(X2, X1)) + 
      geom_tile(aes(fill = value),colour = "white") + 
      scale_fill_gradient(low = "white",     high = "steelblue", limits = c(0,1)) +
      opts(axis.ticks = theme_blank(), axis.text.x = theme_text(angle = 30, hjust = 1), 
      title = title.string) +
      labs(x = "Regulator", y = "Regulated") +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) 
    )
    print(p)
}

##-----------------------------------------------------------------------------------------------
## Calculate and plot distributions on num of parents     ---------------------------------------
## excluding self interactions                            ---------------------------------------
##-----------------------------------------------------------------------------------------------
.plotDistribParents.LargeMat <- function(probMat, numParentsMat, geneNames, numRegs = NULL, max.parents = 4, numHighLinks = 4)
{
  nonTF.net <- T
  if (is.null(numRegs)){
    numRegs   <- length(geneNames)
    nonTF.net <- F
  }
  numNonReg <- length(geneNames)

  numberParents <- dim(numParentsMat)[2] - 1
  ## Reduce numParentsMat to max number of parents if necessary
  if (numRegs - (max.parents+1) > 1)
  {
    max.parentsOrMore               <- rowSums(numParentsMat[,(max.parents+1):numRegs])
    numParentsProbs.final           <- cbind(numParentsMat[,1:max.parents], max.parentsOrMore)
    col.names.aux                   <- c(0:(max.parents-1), paste(">", max.parents-1))
#     colnames(numParentsProbs.final) <- col.names.aux
  }else{
    numParentsProbs.final           <- numParentsMat
    col.names.aux                   <- as.character(c(0:numberParents))
    colnames(numParentsProbs.final) <- as.character(c(0:numberParents))
  }
  rownames(numParentsProbs.final) <- geneNames
  numParentsProbs.final           <- data.frame(rownames(numParentsProbs.final), rep("Same", numNonReg), rep("Posterior Probability # Parents", numNonReg), numParentsProbs.final)
  colnames(numParentsProbs.final) <- c("GeneNames", "Same", "ProbType", col.names.aux )
  
  matHighLinks <- matrix(0, numNonReg, numHighLinks)
  if (!nonTF.net){diag(probMat) <- 0}else{probMat[,1] <- 0}
  
  for(i in 1:numNonReg)
  {
    matHighLinks[i,] <- sort(probMat[i,], T)[1:numHighLinks]
  }

  matHighLinks <- data.frame(geneNames,  rep("Same", numNonReg), rep(paste("Top", numHighLinks, "Regulator Probabilities"), numNonReg), matHighLinks)
  colnames(matHighLinks) <- c("GeneNames", "Same", "ProbType",paste("#", 1:numHighLinks))

  ## Filter out those with prob[0 parents] > 0.6
  remove.these <- which(numParentsProbs.final[,4] > 0.5)
  # Remove only if more than 20 would be left
  if (numNonReg-length(remove.these) > 20)
  {
     matHighLinks          <- matHighLinks[-1*remove.these,]
     numParentsProbs.final <- numParentsProbs.final[-1*remove.these,]
  }

  ## How many plots?
  numRowsToPlot  <- dim(matHighLinks)[1]
  max.chunk.size <- 30
  numChunks      <- ceiling(numRowsToPlot/max.chunk.size)
  for (i in 1:numChunks)
  {
    first.elem <- max.chunk.size*(i-1) + 1
    last.elem  <- max.chunk.size*i
    if (i == numChunks){last.elem <- numRowsToPlot}
    numParentsProbs.final_i <- numParentsProbs.final[first.elem:last.elem,]
    matHighLinks_i          <- matHighLinks[first.elem:last.elem,]
    
    all.m <- suppressMessages(melt(list(numParentsProbs.final_i, matHighLinks_i )))
    (p <- ggplot(all.m, aes(variable, GeneNames)) + 
      geom_tile(aes(fill = value),colour = "white") + 
      scale_fill_gradient(low = "white",     high = "steelblue", limits = c(0,1)) +
      opts(title = "Network Uncertainty", axis.ticks = theme_blank(), axis.text.x = theme_text(angle = 30, hjust = 1)) +
      labs(x = "",y="") +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +  
      facet_grid( Same ~ ProbType, scales = "free")
    )
    ## Save
    print(p)
  }
}

##-----------------------------------------------------------------------------------------------
## Calculate and plot distributions on num of parents     ---------------------------------------
## excluding self interactions                            ---------------------------------------
##-----------------------------------------------------------------------------------------------
.marginalUncertaintyPlot <- function(gamma.chain, geneNames, numRegs = NULL, max.parents = 4, numHighLinks = 4)
{
  nonTF.net <- T
  if (is.null(numRegs)){
    numRegs   <- length(geneNames)
    nonTF.net <- F
  }
  numNonReg <- length(geneNames)
  lenChain <- dim(gamma.chain)[1]
  ## Make count number of parents matrix 
  numParentsMat <- matrix(0, numNonReg, numRegs)  
  for (i in 1:lenChain)
  {
    net.i <- matrix(gamma.chain[i,], nrow = numNonReg)
    ## Remove self interactions
    if (!nonTF.net){diag(net.i) <- 0}else{net.i[,1] <- 0}
    parents.i   <- rowSums(net.i)
    for (j in 1:numNonReg)
    {
	numParentsMat[j, parents.i[j] + 1] <- numParentsMat[j, parents.i[j] + 1] + 1
    }
  }
  ## Divide by chain length
  numParentsMat <- numParentsMat/lenChain 
#   write.table(numParentsMat, "ProbNumParents_Full")
  probMat <- matrix(colMeans(gamma.chain), ncol = numRegs)
  if(!nonTF.net){numRegs <- NULL}
  .plotDistribParents.LargeMat(probMat, numParentsMat, geneNames, numRegs , max.parents , numHighLinks)
}


## -----------------------------------------------------------------------------------------------
##  Plot distribution and prior lambda --------------------------------------------------
## -----------------------------------------------------------------------------------------------
.plotDistribWithGammaPrior <- function(lambda, prior.a, prior.b, numGenes, main.text ="Prior and posterior lambda", max.x = 500)
{
  all.dens <- matrix(0, numGenes, 512)
  all.x    <- matrix(0, numGenes, 512)
  ## Loop through genes plotting densities
  for (i in 1:numGenes)
  {
    dens.i       <- density(lambda[,i], na.rm=T, from = 0, n = 512)
    all.dens[i,] <- dens.i$y
    all.x[i,]    <- dens.i$x
  } 

  #max.x  <- max(all.x)*1.1
  # max.x  <- 500
  max.y  <- max(all.dens)*1.2
  x.seq  <- seq(from=0, to = max.x, length.out = 10000)
  prior1 <- list(x = x.seq ,
		  y = dgamma(x.seq, shape = prior.a, rate = prior.b))

  ## Plot priors
  plot(x = prior1$x, y = prior1$y, col="red", lty=2, type="l", 
	ylim = c(0,max.y), lwd = 2, main=main.text, xlab = "Degrees of freedom", 
	ylab = "Probability density" )

  ## Loop through genes plotting densities
  for (i in 1:numGenes)
  {
    lines(y = all.dens[i,], x=  all.x[i,], lwd = 0.7)
  } 
  ## Add legend 
  legend("topright", c(paste("prior ", paste(prior.a, prior.b)), "posterior"), col=c("red", "black"), lty=c(2,1)) 
}

##-----------------------------------------------------------------------------------------------
## Pareto  pdf     ------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------
.dpareto <- function(tau,a, b)
{
  a*b^(-a)*tau^(a-1)
}

.Fpareto <- function(a, b, num.points=10000)
{
  tau.vals <- matrix(seq(from=0, to = b, length.out = num.points),nrow=1)
  #print(tau.vals)
  list(y=apply(tau.vals, 2, .dpareto, a=a, b=b),x=tau.vals)
}

## -----------------------------------------------------------------------------------------------
##  Calculate tau when link on ------------------------------------------------------------
## -----------------------------------------------------------------------------------------------
.tauON <- function(tau, gamma.vals)
{
  tau.on     <- tau*gamma.vals
  tau.on[tau.on ==0 ]<-NA 
  tau.on
}

## -----------------------------------------------------------------------------------------------
##  Plot distribution and prior tau ------------------------------------------------------------
## -----------------------------------------------------------------------------------------------
.plotDistribAndPriorTau.sepSelf <- function(tau1, link.prob.1, prior.param)
{
  ## Prior parameter values
  a_pareto <- prior.param[12]
  truncMe  <- prior.param[ 6]

  ## Separate self interaction elements 
  num.links     <- length(link.prob.1)
  numGenes      <- sqrt(num.links)
  self.interact <- diag(matrix(1:num.links, numGenes, numGenes))
  no.self       <- 1:num.links
  no.self       <- no.self[-self.interact]

  ## Find max value of link prob to use as plotting limits for tau
  non.diag.links  <- link.prob.1[no.self]
  pos.max.nonSelf <- which(non.diag.links == max(non.diag.links))[1]
  i.max           <- no.self[pos.max.nonSelf]

  ### Plot NON self interactions -------------
  dens.kk <- density(tau1[,i.max ], na.rm=T)
  max.x  <- max(dens.kk$x)
  max.y  <- max(dens.kk$y)*1.2
  # prior1 <- gaga(prior.param$a_spline, prior.param$c_spline, prior.param$d_spline,  interval=c(0,max.x))
  prior1 <- .Fpareto(a_pareto, truncMe)
  ## Plot priors
  # plot(x = c(0, truncMe), y = c(1/truncMe,1/truncMe), col="green",lty=2,type="l", xlim=c(0,truncMe),ylim = c(0,max.y),lwd = 2, main="Prior and posterior tau (with prob>0.35) NON-SELF")
  plot(x = prior1$x, y = prior1$y, xlab = "tau", ylab = "Probability Density", col="red",lty=2,type="l", ylim = c(0,max.y),lwd = 2, main="Prior and posterior tau (with prob>0.35)")

  ## Loop through genes plotting densities
  for (i in no.self )
  {
    ## tau1
    if (link.prob.1[i]>0.35)
    {
      dens.i <- density(tau1[,i], na.rm=T, from=0)
      lines(dens.i, lwd = 0.7)
    }
  } 
  # legend("topright", c(paste("prior ", paste(prior.param$a_spline, prior.param$c_spline, prior.param$d_spline)), "posterior"), col=c("green", "black"), lty=c(2,1)) 
  legend("topright", c(paste("prior pareto a =", a_pareto, ", trunc =", truncMe), "posterior"), col=c("red", "black"), lty=c(2,1)) 


  ### Plot SELF interactions if on -----------------
  if (link.prob.1[1]!=0){
      dens.kk <- density(tau1[,1], na.rm=T)
      max.x  <- max(dens.kk$x)
      max.y  <- max(dens.kk$y)*1.5
      # prior1 <- gaga(prior.param$a_spline, prior.param$c_spline, prior.param$d_spline,  interval=c(0,max.x))
      prior1 <- .Fpareto(3*a_pareto, truncMe)

      ## Plot priors
    #   plot(x = c(0, truncMe), y = c(1/truncMe,1/truncMe), col="green",lty=2,type="l", xlim=c(0,truncMe),ylim = c(0,max.y),lwd = 2, main="Prior and posterior tau SELF")
      plot(x = prior1$x, y = prior1$y, xlab = "tau", ylab = "Probability Density", col="red",lty=2,type="l", ylim = c(0,max.y),lwd = 2, main="Prior and posterior tau SELF")

      ## Loop through genes plotting densities
      for (i in self.interact )
      {
	dens.i <- density(tau1[,i], na.rm=T, from=0)
	lines(dens.i, lwd = 0.7)
      } 
      legend("topright", c(paste("prior pareto a = 3x", a_pareto, ", trunc =", truncMe), "posterior"), col=c("red", "black"), lty=c(2,1)) 
  }

}

## -----------------------------------------------------------------------------------------------
##  Reconstruction   -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------
.plotSplinesFunctions <- function(     allfs1, all.fsqr.1,       Full_F_sqr.1, link.probs1, 
				  data.matrix,  GeneNames, mu.mean, plot.thresh = 0.35) 
# tau1, tau2, allfs1, allfs2, link.probs1, link.probs2, plot.thresh, scaled.data, chainNames, GeneNames, op)
{
  sd.f        <- sqrt(all.fsqr.1 - allfs1^2)
  scaled.data <- t(scale(t(data.matrix)))
  time.m      <- dim(scaled.data)[2] -1
  genes       <- dim(scaled.data)[1]
  probs1      <- matrix(link.probs1, genes, genes)
  index.flat  <- matrix(  1:genes^2, genes, genes)
  diag(probs1) <-1 
  for (i in 1:genes)
  {  
    whichParents <- which(probs1[i,]>plot.thresh)
    numPlots     <- length(whichParents)
#     total.row.plot <- ceiling((length(whichParents) + 1) /2)  
    ## Need either 3 or 4 rows
    if (numPlots < 4){plot.rows <- 2}else{plot.rows <- 3}
#     par(mfrow=c(total.row.plot, 4), mar = c(2,2,4,2), oma = c(0.1,0.1,3.5,0.1))
    par(mfrow=c(plot.rows, plot.rows), mar = c(4,4,3.5,0.6), oma = c(1,1,3.5,1))
    ## Place intersection first
    plotMe.indx <- sort(which(probs1[i,]>plot.thresh), T)
    ## Place self interaction first 
    plotMe.indx <- c(i, plotMe.indx[-1*which(plotMe.indx == i)])
    for (j in plotMe.indx)
    {
      tot.index <- genes*(i-1) + j   
      .plot.reconst(probs1[i,j], allfs1[,tot.index], sd.f[,tot.index], scaled.data, i, j, GeneNames, plot.thresh)      
    }
    meanF     <- rowSums(allfs1[, genes*(i-1) + (1:genes)])
    sd.F      <- sqrt(Full_F_sqr.1[,i] - meanF^2)
    ## Add posterior mean intercept
    meanF     <- meanF + mu.mean[i]
    ## 2 sd lines
    upperLine <- meanF + 2*sd.F
    lowerLine <- meanF - 2*sd.F
    ## Plot y lims
    y.lims <- range(scaled.data[i,-(time.m+1)]) + c(-0.18,0.1)
#     x.lims <- range(scaled.data[j,-total.time]) + c(0,0.18) 
#  plotSplinesFunctions(allfs1, all.fsqr.1,       Full_F_sqr.1, link.probs1, Locke1hLoggedNoNoise.data,  rownames(Locke1hLoggedNoNoise.data))

    # Full reconstruction
    plot(y = upperLine, x = 1:time.m,col = "red", 
    lty = 2, lwd = 2, type = "l", ylim = y.lims, main="", 
    xlab = "time", 
    ylab = paste(GeneNames[i],"(t + 1)") )   

    lines(y = lowerLine, x = 1:time.m, col = "red", lty = 2, lwd = 2)   
    polygon(c(1:time.m, rev(1:time.m)), 
	    c(upperLine, rev(lowerLine)), col = "mistyrose", border = NA)
    lines(y = meanF, x = 1:time.m, col = "red", lty = 1, lwd = 2)   
    points(scaled.data[i,-1])

#     title(main.text, cex.main = 1,   font.main= 2, col.main= "blue")
    legend("bottomright", c(paste("mean F(t)",sep=""), "2xsd"), 
	    col=c("red"), bty = "n", lty = c(1,2), lwd = c(2,2.5))
    title("Full reconstruction",  cex.main = 2,   font.main= 3, col.main= 5)
    mtext(paste("Regulators of:", GeneNames[i]), side=3, line=1.5, font=2, cex=2, col='red',  outer=TRUE)

#     plot(scaled.data[i,-1], main = "", xlab = "Time", ylab = GeneNames[i])
#     meanF <- rowSums(allfs1[, genes*(i-1) + (1:genes)])
#     sd.F  <- sqrt(Full_F_sqr.1 - meanF^2)
#     lines(meanF, col="red", lwd = 1.8)
#     par(op)
  }
}


# Plot data and spline function
.plot.reconst <- function(probs1.ij, allfs1.ij, sd.f.ij, scaled.data, i, j, GeneNames, plot.thresh)
{
  total.time <- dim(scaled.data)[2]
  if (i==j)
  {
    main.text <- "Self Interaction"  
    sub.title <- ""
  }else{
    main.text <- paste(GeneNames[j]," regulates", GeneNames[i], "\n (prob = ", round(probs1.ij,1), ")", sep="")
  }

  kk <- order(scaled.data[j,-total.time])
  ## 2 sd lines
  upperLine <- allfs1.ij[kk] + 2*sd.f.ij[kk]
  lowerLine <- allfs1.ij[kk] - 2*sd.f.ij[kk]
  ## Plot y lims
  y.lims <- range(scaled.data[i,-total.time]) + c(-0.18,0.1)
  x.lims <- range(scaled.data[j,-total.time]) + c(0,0.18) 

  plot(y = upperLine, x = t(scaled.data[j,kk]), col = "red", 
  lty = 2, lwd = 2, type = "l", ylim = y.lims, , xlim = x.lims, main="", 
  xlab = paste(GeneNames[j],"(t)"), 
  ylab = paste(GeneNames[i],"(t + 1)") )   

  lines(y = lowerLine, x = t(scaled.data[j,kk]), col = "red", lty = 2, lwd = 2)   
  polygon(c(t(scaled.data[j,kk]), rev(t(scaled.data[j,kk]))), 
	  c(upperLine, rev(lowerLine)), col = "mistyrose", border = NA)
  lines(y = allfs1.ij[kk], x = t(scaled.data[j,kk]), col = "red", lty = 1, lwd = 2)   
  points(y= t(scaled.data[i,-1]), x=t(scaled.data[j,-total.time]))

  title(main.text, cex.main = 1,   font.main= 2, col.main= "blue")
  legend("bottomright", c(paste("mean f(", GeneNames[j],")",sep=""), "2xsd"), 
	  col=c("red"), bty = "n", lty = c(1,2), lwd = c(2,2.5))
}

