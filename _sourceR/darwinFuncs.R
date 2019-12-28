# --- function to regress responses for an item onto factor scores for g -----------
regong <- function(anitem, gfascores) {
  reg1 <- lm(anitem~gfascores)
  resids <- reg1$residuals
  return(resids+mean(anitem))
}

# --- function for exploration plot, for each political groups separately ----------
plotexploration <- function(topcons, x, y, z, xlabarg, ylabarg, zlabarg, colpal, paldir, t4plot) {
  tcp <- ggplot(topcons, aes(x,y)) +
    geom_point(aes(colour=z),size=4,shape = 19,
               position=position_jitter(width=0.01, height=0.01)) +
    scale_color_distiller(palette=colpal,direction=paldir) +
    theme(plot.title = element_text(size = 18),
          panel.background = element_rect(fill = "white", colour = "grey50"),
          text = element_text(size=14)) +
    geom_hline(yintercept=mean(y), linetype="dashed", color = "darkgray") +
    geom_vline(xintercept=mean(x), linetype="dashed", color = "darkgray") +
    xlim(0.8,7.2) + ylim(0.8,7.2) +
    labs(colour=paste(zlabarg,"\n",sep=""), x=xlabarg, y=ylabarg, title=t4plot)
  invisible(print(tcp))
  return(tcp)
}

#function to fit and print a measurement model ------------------------------------------------------------
fitMM <- function (humspecmod,this1) {
  fithumspec <- sem(humspecmod, data=this1)
  summary(fithumspec, standardized=TRUE, fit.measures=TRUE)
  mod_ind <- modificationindices(fithumspec)
  print(head(mod_ind[order(mod_ind$mi, decreasing=TRUE), ], 4))
  return(fitMeasures(fithumspec))
}

#function to create and run a single indicator model and return omega and so on --------------------------
#  indicators are provided in character vector
#  pairwise correlations are optionally provided as a list of pairwise character vectors
singleindmodel <- function(indicators,paircorrs,dat) {
  pcs <- paste(unlist(lapply(paircorrs,FUN=function(x) paste(x[1], "~~", x[2], "\n"))),collapse=" ")
  SImod <- paste("factor =~", paste(indicators, collapse=" + "), "\n",
                 pcs, 
                 "phantfac <~", paste(paste("1*",indicators,sep=""), collapse=" + "), "\n",
                 "factor ~~ 0*phantfac")
  
  fitSImod <- sem(SImod, dat[,indicators], estimator="ML")
  ParSImod <- parameterEstimates(fitSImod, standardized=TRUE)
  LoadingsSImod <- ParSImod[1:length(indicators), "std.all"]
  ErrorvarSImod <- 1 - LoadingsSImod^2
  ImpliedCorrSImod <- lavTech(fitSImod, what='cor.lv')
  OmegaSImod <- ImpliedCorrSImod[[1]][1,2]^2            #Squared correlation between factor and phantom variable
  
  varSImod <- var(apply(dat[,indicators], MARGIN=1, FUN=mean)) #variance of composite
  SDSImod <- sqrt(varSImod)               #SD of composite, to pass back for analysis
  eSImod <- (1-OmegaSImod)*varSImod       #error term of single indicator
  return(listN(OmegaSImod,eSImod,SDSImod))
}

# Miscellaneous functions -------------------------------------------------------------------------------
#http://stackoverflow.com/questions/21011672/automatically-add-variable-names-to-elements-of-a-list
listN <- function(...){
  anonList <- list(...)
  names(anonList) <- as.character(substitute(list(...)))[-1]
  anonList
}
revscore <- function (x,mm) {  #this reverse scores a scale for items of reversed polarity
  return ((mm+1)-x)            #same as fixscore but to be conceptually clear it's a separate function
}
fixscore <- function (x,mm) {  #this fixes scale so SD=1 and SA=7 irrespective of polarity
  return ((mm+1)-x)
}
interleave <- function(x){     #function to grab a row and interleave with parentheses for printing
  hlx <- length(x)/2
  retstr <- paste(sapply(1:hlx, FUN=function(i) paste(" & ",as.character(x[i])," & (",as.character(x[i+hlx]),")",sep="")),collapse="",sep="")
}
