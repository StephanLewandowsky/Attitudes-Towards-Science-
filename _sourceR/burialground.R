

smallCorrel <- c("
                 evoFac   =~ evo
                 mwevoFac =~ mwevo
                 mwnatFac =~ mwnat        
                 mwequFac =~ mwequ
                 religFac =~ relig
                 consFac  =~ cons

                  evo ~~ ",  evoSI$eSImod,   "*evo",
                 "mwevo ~~", mwevoSI$eSImod, "*mwevo",
                 "mwnat ~~", mwnatSI$eSImod, "*mwnat",
                 "mwequ ~~", mwequSI$eSImod, "*mwequ",
                 "relig ~~", religSI$eSImod, "*relig",
                 "cons ~~",  consSI$eSImod, "*cons"
)
fitsmallCorrel <- sem(smallCorrel, compositedarwin, std.lv=TRUE, estimator="ML")
summary(fitsmallCorrel,standardized=TRUE, fit.measures=TRUE)
lavInspect(fitsmallCorrel, what = "cor.lv")
#fitsmallCorrelFM <- measEq.syntax(smallCorrel, compositedarwin, std.lv=TRUE, estimator="ML",group="gender") # http://lavaan.ugent.be/tutorial/groups.html
x11()
semPaths(fitsmallCorrel, what="std", title =FALSE, curvePivot = TRUE,residuals=FALSE, structural=TRUE, layout="circle")


#* ------------- Now model various scientific propositions, each on its own first -------------
# Rejection of CAM : with mwnatFac as well as mwevoFac as predictor, cov matrix not positive definite.
# cov matrix not PD even without mw constructs as predictors
modelCAM <- c("  camFac   ~ fmFac +  religFac 
                  mwevoFac ~~ 0*consFac
                  mwevoFac ~~ 0*vaxFac
                  consFac  ~~ 0*vaxFac
                  mwnatFac =~ mwnat
                  fmFac    =~ fm
                  evoFac   =~ evo
                  camFac   =~ cam
                  mwevoFac =~ mwevo

                  consFac  =~ cons
                  mwequFac =~ mwequ
                  religFac =~ relig
                  vaxFac   =~ vax

                  fm ~~",     fmSI$eSImod,    "*fm",
              "evo ~~ ",  evoSI$eSImod,   "*evo",
              "cam ~~",   camSI$eSImod,   "*cam",
              "mwevo ~~", mwevoSI$eSImod, "*mwevo",
              "mwnat ~~", mwnatSI$eSImod, "*mwnat",
              "cons ~~",  consSI$eSImod, "*cons",
              "mwequ ~~", mwequSI$eSImod, "*mwequ",
              "relig ~~", religSI$eSImod, "*relig",
              "vax ~~",   vaxSI$eSImod,    "*vax"
)
fitCAM <- sem(modelCAM, compositedarwin, std.lv=TRUE, estimator="ML")
summary(fitCAM,standardized=TRUE, fit.measures=TRUE)
lavInspect(fitCAM, "cov.lv")
fitCAMgof <- fitMeasures(fitCAM)
modificationIndices(fitCAM, sort. = TRUE, maximum.number = 4) 

# Evolution : with mwnatFac as well as mwevoFac as predictor, cov matrix not positive definite.
modelEvo <- c("  evoFac   ~ fmFac +  religFac + mwequFac + mwevoFac + camFac + consFac 
                  mwevoFac ~~ 0*religFac
                  mwevoFac ~~ 0*vaxFac
                  mwevoFac ~~ 0*fmFac
                  consFac  ~~ 0*vaxFac
                  mwnatFac =~ mwnat
                  fmFac    =~ fm
                  evoFac   =~ evo
                  camFac   =~ cam
                  mwevoFac =~ mwevo

                  consFac  =~ bc * cons
                  mwequFac =~ mwequ
                  religFac =~ relig
                  vaxFac   =~ vax
                  bc > 0  

                  fm ~~",     fmSI$eSImod,    "*fm",
              "evo ~~ ",  evoSI$eSImod,   "*evo",
              "cam ~~",   camSI$eSImod,   "*cam",
              "mwevo ~~", mwevoSI$eSImod, "*mwevo",
              "mwnat ~~", mwnatSI$eSImod, "*mwnat",
              "cons ~~",  consSI$eSImod, "*cons",
              "mwequ ~~", mwequSI$eSImod, "*mwequ",
              "relig ~~", religSI$eSImod, "*relig",
              "vax ~~",   vaxSI$eSImod,    "*vax"
)
fitEvo <- sem(modelEvo, compositedarwin, std.lv=TRUE, estimator="ML")
summary(fitEvo,standardized=TRUE, fit.measures=TRUE)
lavInspect(fitEvo, "cov.lv")
fitEvogof <- fitMeasures(fitEvo)
modificationIndices(fitEvo, sort. = TRUE, maximum.number = 4)


# Vaccination : with mwnatFac as well as mwevoFac as predictor, cov matrix not positive definite.
modelVax <- c("  vaxFac   ~ fmFac +  mwequFac + consFac 
                  mwevoFac ~~ 0*religFac
                  mwevoFac ~~ 0*vaxFac
                  mwevoFac ~~ 0*fmFac
                  consFac  ~~ 0*vaxFac
                  mwnatFac =~ mwnat
                  fmFac    =~ fm
                  evoFac   =~ evo
                  camFac   =~ cam
                  mwevoFac =~ mwevo
                  consFac  =~ bc * cons
                  mwequFac =~ mwequ
                  religFac =~ relig
                  vaxFac   =~ vax
                  bc > 0

                  fm ~~",     fmSI$eSImod,    "*fm",
              "evo ~~ ",  evoSI$eSImod,   "*evo",
              "cam ~~",   camSI$eSImod,   "*cam",
              "mwevo ~~", mwevoSI$eSImod, "*mwevo",
              "mwnat ~~", mwnatSI$eSImod, "*mwnat",
              "cons ~~",  consSI$eSImod, "*cons",
              "mwequ ~~", mwequSI$eSImod, "*mwequ",
              "relig ~~", religSI$eSImod, "*relig",
              "vax ~~",   vaxSI$eSImod,    "*vax"
)
fitVax <- sem(modelVax, compositedarwin, std.lv=TRUE, estimator="ML")
summary(fitVax,standardized=TRUE, fit.measures=TRUE)
fitVaxgof <- fitMeasures(fitVax)
modificationIndices(fitVax, sort. = TRUE, maximum.number = 4)
#Owing to warning "covariance matrix of latent variables is not positive definite", the mwnat LV had to be removed
lavInspect(fitVax, what = "cor.lv")
eigen( inspect(fitVax, "cov.lv") )$values 
x11()
semPaths(fitVax, what="std", title =FALSE, curvePivot = TRUE,residuals=FALSE,
         structural=TRUE, layout="tree2",rotation=2)


#2nd order factor vax as function of all of politics combined
modelVax2O <- c("vaxFac  ~ allPolFac + fmFac
                  allPolFac =~ fmFac + consFac + religFac
                  
                  fmFac    =~ bc*fm
                  consFac  =~ bc*cons
                  bc > 0
                  religFac =~ relig
                  vaxFac   =~ vax

                  fm ~~",     fmSI$eSImod,    "*fm",
                "cons ~~",  consSI$eSImod, "*cons",
                "relig ~~", religSI$eSImod, "*relig",
                
                "vax ~~",   vaxSI$eSImod,    "*vax"
)
fitVax2O <- sem(modelVax2O, compositedarwin, std.lv=TRUE, estimator="ML")
summary(fitVax2O,standardized=TRUE, fit.measures=TRUE)
fitVax2Ogof <- fitMeasures(fitVax2O)
#modificationIndices(fitVax2O, sort. = TRUE, maximum.number = 4)

#reliability(fitVax2O) # Should provide a warning for the endogenous variables
reliabilityL2(fitVax2O, "allPolFac")




#*-----------------------------------------------------------------------------------
# Vax & evolution & CAM rejection: first full then slightly reduced (i.e., some covariances reduced to 0)  
modelVaxEvoCAMfull <- c("  vaxFac   ~ fmFac +  mwequFac + consFac
                            evoFac   ~ fmFac +  religFac + mwequFac + mwevoFac  + consFac 
                            camFac   ~ fmFac +  religFac
                            
                      mwnatFac =~ mwnat
                      fmFac    =~ fm
                      evoFac   =~ evo
                      camFac   =~ cam
                      mwevoFac =~ mwevo
                      mwequFac =~ mwequ
                      religFac =~ relig
                      vaxFac   =~ vax
                      consFac  =~ cons
                      
                      fm ~~",     fmSI$eSImod,    "*fm",
                        "evo ~~ ",  evoSI$eSImod,   "*evo",
                        "cam ~~",   camSI$eSImod,   "*cam",
                        "mwevo ~~", mwevoSI$eSImod, "*mwevo",
                        "mwnat ~~", mwnatSI$eSImod, "*mwnat",
                        "mwequ ~~", mwequSI$eSImod, "*mwequ",
                        "cons ~~",  consSI$eSImod, "*cons",
                        "relig ~~", religSI$eSImod, "*relig",
                        "vax ~~",   vaxSI$eSImod,    "*vax"
)
fitVaxEvoCamfull <- sem(modelVaxEvoCAMfull, compositedarwin, std.lv=TRUE, estimator="ML")
summary(fitVaxEvoCamfull,standardized=TRUE, fit.measures=TRUE)


modelVaxEvoCAM <- c("vaxFac   ~ fmFac +  mwequFac + consFac
                      evoFac   ~ fmFac +  religFac + mwequFac + mwevoFac  + consFac 
                      camFac   ~ fmFac +  religFac
                      
                      evoFac ~~ 0*camFac
                 
                         mwnatFac =~ mwnat
                         fmFac    =~ fm
                         evoFac   =~ evo
                         camFac   =~ cam
                         mwevoFac =~ mwevo
                         mwequFac =~ mwequ
                         religFac =~ relig
                         vaxFac   =~ vax
                         consFac  =~ bc*cons
                         bc > 0
                         
                         fm ~~",     fmSI$eSImod,    "*fm",
                    "evo ~~ ",  evoSI$eSImod,   "*evo",
                    "cam ~~",   camSI$eSImod,   "*cam",
                    "mwevo ~~", mwevoSI$eSImod, "*mwevo",
                    "mwnat ~~", mwnatSI$eSImod, "*mwnat",
                    "mwequ ~~", mwequSI$eSImod, "*mwequ",
                    "cons ~~",  consSI$eSImod, "*cons",
                    "relig ~~", religSI$eSImod, "*relig",
                    "vax ~~",   vaxSI$eSImod,    "*vax"
)

fitVaxEvoCam <- sem(modelVaxEvoCAM, compositedarwin, std.lv=TRUE, estimator="ML")
summary(fitVaxEvoCam,standardized=TRUE, fit.measures=TRUE)
modificationIndices(fitVaxEvoCam, sort. = TRUE, maximum.number = 6)

aovfullvfin <- anova(fitVaxEvoCamfull,fitVaxEvoCam)
aovfullvfindf<- aovfullvfin[["Df"]][2]- aovfullvfin[["Df"]][1]

finalmodelgof <- fitMeasures(fitVaxEvoCam)
modificationIndices(fitVaxEvoCam, sort. = TRUE, maximum.number = 4)
lavInspect(fitVaxEvoCam, what = "cor.lv")
eigen( inspect(fitVaxEvoCam, "cov.lv") )$values 

pdf(file=paste(figdir,"finalSEM.pdf",sep="/"),height=8,width=8) 
semPaths(fitVaxEvoCam, what="std", title =FALSE, curvePivot = TRUE,residuals=FALSE, intercepts=FALSE,
         sizeLat=12,sizeLat2=12,
         curvature=5, nCharNodes=0, edge.label.cex=1.2,mar=c(3,10,3,3),
         structural=TRUE, layout="tree",rotation=2, 
         label.scale.equal=TRUE,label.prop=0.95,
         nodeLabels=c("M&W\nnaturally \n different","Free Market","Evolution","Reject\nCAM",
                      "M&W\nevolved \ndifferently","M&W\nthe same","Religiosity","Vaccinations","Conservatism"))
dev.off()



#raw4fa <- darwin1.5 %>% select(contains("Q")) %>% select(contains("."))

#Factor analysis of the data using various ways of extracting general factor
#see e.g. https://www.personality-project.org/revelle/presentations/general.pdf 
# https://psyarxiv.com/2y3w9/ 
#darwinfa <- fa(r = raw4fa, nfactors = 4, rotate="bifactor")
#darwinfa
#use of omega for general factor is recommended
#darwinfa3 <- omega(raw4fa,nfactors=6)
#get factor scores for g, regress all items onto that and convert scores to residuals
#gfacscores <- darwinfa3$scores[,"g"]
#darwin1.5 <- darwin1.5 %>% mutate_at(vars(QC.1:QG.5),regong)

#rescale all variables 0-1
#rescale7   <- darwin2 %>% select(contains("Q")) %>% select(contains(".")) %>% lapply(., FUN=function(x) scales::rescale(x,to=c(0,1))) %>% as.data.frame  #mutate_all(funs(. -1))  %>% mutate_all(funs(. /6))
#rescale100 <- darwin2 %>% select(contains("QB_")) %>% mutate_all(funs(. /100))
darwin2norescale <- darwin2  #keep a copy for political extreme analysis at the very end that is not scaled 0-1
#darwin2[,names(rescale7)] <- rescale7
#darwin2[,names(rescale100)] <- rescale100

# #examine effects of affirmation bias -- so split on affirmation bias not politics
# ggpanels2 <- list()
# k <- 0
# for (qqs in c(.25,.5)) {
#   topcons <- meanIndicators %>% filter(affirmbias > quantile(affirmbias,1-qqs))  
#   toplibs <- meanIndicators %>% filter(affirmbias < quantile(affirmbias,qqs)) 
#   k<-k+1
#   ggpanels2[[k]] <- plotexploration(toplibs, toplibs$MWevomean, toplibs$MWequmean,toplibs$MWnatmean, 
#                                     "Men and women evolved differently",
#                                     "Men and women are equal",
#                                     "Nat Diff","YlGnBu",1,
#                                     paste("Affirmation bias (bottom ",qqs*100,"%)",sep=""))
#   k<-k+1
#   ggpanels2[[k]] <- plotexploration(topcons, topcons$MWevomean, topcons$MWequmean, topcons$MWnatmean, 
#                                     "Men and women evolved differently",
#                                     "Men and women are equal",
#                                     "Nat Diff","YlGnBu",1,
#                                     paste("Affirmation bias (top ",qqs*100,"%)",sep=""))
# }
# x11(width=12,height=8)
# gridplt2 <- do.call("grid.arrange", c(ggpanels2, nrow=2))
# ggsave("affirmbiasgrid.pdf", plot=gridplt2, width=12, height=8)
# 
# 
# 
# # check out based on top and bottom political quartile
# topconsptr <- consmean > quantile(consmean,.75)  
# toplibsptr <- consmean < quantile(consmean,.25) 
# 
# #pdf(file="polgenderextremehisto.pdf",height=6,width=6) 
# par(mfrow=c(2,2))
# 
# tt<-darwin2norescale[toplibsptr,] %>% select("QI.4") %>% 
#   hist(las=1,xlim=c(0,1),xlab="Average score",main=clusterlabels["J"],col="light gray")
# 
# 
# tt<-darwin2norescale[toplibsptr,] %>% select(num_range(paste("Q","J",".",sep=""),1:5)) %>% melt %>% select(value) %>%
#                hist(las=1,xlim=c(0,1),xlab="Average score",main=clusterlabels["J"],col="light gray")
# 
# dev.off()
