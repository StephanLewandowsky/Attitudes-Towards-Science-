#' ---
#' title: "Analysis of Qualtrics Data L'sky - Woike - Oberauer"
#' ---
#analyze qualtrics darwinian evolution representative sample 
rm(list=ls())
library(lattice)
library(stargazer)
library(tidyverse)
library(lme4)
library(lavaan)
library(semPlot)
library(semTools)
library(RColorBrewer)
library(gridExtra)
library(stringr)
library(Hmisc)
library(reshape2)
library(psych)
library(scales)
library(summarytools) #contains descr()

##--- correct affirmation bias or not?
correctaffirmbias <- 1

#define working directories
texdir   <- "../_outputTex"
figdir   <- "../_figures"
inputdir <- "../_input"

##--- Begin by reading data and variable names -----------------------------------------------------------------------
if (rev(str_split(getwd(),"/")[[1]])[1] == "_paper") {
  setwd("../_sourceR")
}
source("darwinFuncs.r")
darwin <- read.csv(paste(inputdir,"full1200downloadedAnonymized.csv",sep="/"),header=TRUE,row.names=NULL) 
#row 2&3 (with verbose names) manually deleted from Q file. 

#read variable names for table of raw responses
vn <- read.csv(paste(inputdir,"variableNames.csv",sep="/"),header=TRUE,row.names=NULL,stringsAsFactors = FALSE) 

darwin1.5 <- darwin %>% filter(QB_15 == 20) %>%       #must choose 20 from the AFQ slider (this eliminates NA softlaunch subject)
                    filter(QAC == 4) %>%              #table is not an animal
                    #filter(Q_TotalDuration<1801) %>%    #30 minutes should be ample
                    select(contains("Q"),contains("CRT")) %>% 
                    select(-contains("Q_Tot"),-contains("QZ"),-contains("QX"),-contains("QY"),
                           -contains("Click"),-contains("Submit")) %>%        #drop the Qs that ain't qs.
                    select(-QB_15)                    #drop AFQ

# first fix the Qualtrics-induced scale problems, so SD=1 and SA=7
darwin1.5 <- darwin1.5 %>%  mutate_at(c(paste("QC.",c(1:5),sep=""),
                                   paste("QE.",c(1:5),sep=""),
                                   paste("QF.",c(1:5),sep=""),
                                   paste("QH.",c(1:5),sep=""),
                                   paste("QI.",c(1:5),sep=""),
                                   paste("QJ.",c(1:5),sep=""),
                                   paste("QG.",c(1:5),sep="")),fixscore,mm=14) %>%
                            mutate_at(paste("QD.",c(1:5),sep=""), fixscore,mm=28) 

##--- identify keyhitters before reverse-scoring ----------------------------------------------------------------------
neutral <- 0 #if set to zero, any sequence of identical keys is eliminated. If set to 4, only non-neutral responses are dropped
keyhitters <- NULL
for (cluster in c("C","E","F","D","G")) { #eliminate "H","I","J", which have no reverse scoring
  keyhitters <- cbind(keyhitters,
                      darwin1.5 %>% select(num_range(paste("Q",cluster,".",sep=""),1:5)) %>% 
                                    apply(.,1,FUN=function(x) ifelse((var(x)==0 & mean(x)!=neutral),1,0)) 
                      )
  }
table(rowSums(keyhitters))
#include "50" for sliders in identifying keyhitters
keyhitters <- cbind(keyhitters, darwin1.5 %>% select(contains("QB")) %>% 
                      apply(.,1,FUN=function(x) ifelse((var(x)==0 & mean(x)==50),1,0))) 
table(rowSums(keyhitters))
#darwin1.5 <- filter(darwin1.5,!(rowSums(keyhitters)>0)) #filter(darwin1.5, Qage<50) #

#demographics
males <- table(darwin1.5$Qgender)["1"]
females <- table(darwin1.5$Qgender)["2"]
mage <- round(mean(darwin1.5$Qage),1)
mdage <- round(median(darwin1.5$Qage),1)
minage <- min(darwin1.5$Qage)
maxage <- max(darwin1.5$Qage)
histogram(~Qage|Qgender,data=darwin1.5)
agebygender <- darwin1.5 %>% group_by(Qgender) %>% summarise(mean(Qage)) %>% unlist %>% as.numeric %>% round(.,1)

#*------------------------  some summaries before reverse-scoring -----------------------------------------------------
slider.raw <- darwin1.5 %>% select(contains("QB_"))
#histograms of slider raw responses
label4plot <- vn$shortname[str_detect(vn$qvarname,fixed("QB_"))][1:14]
s4p<-gather(slider.raw,factor_key = TRUE)
s4p$key <- factor(s4p$key,labels=label4plot)
slidp <- ggplot(s4p, aes(value)) + 
                geom_histogram(bins = 10) + 
                xlab("Slider value") + ylab("Frequency") +
                facet_wrap(~key, scales = 'free_x',labeller=label_value)
ggsave(paste(figdir,"slidersSocDarw.pdf",sep="/"), plot=slidp, width=6.5, height=9)

#add table of all responses
itemResppercent <- darwin1.5 %>% select(contains("Q")) %>% select(contains(".")) %>% lapply(table) %>% lapply(as.numeric) %>%
  lapply(FUN=function(x) c(x,round(x/sum(x)*100))) 
#generate latex code for insertion into document
t4l <- NULL
for (i in 1:length(itemResppercent)) {
  mychar <-paste(vn$shortname[vn$qvarname==names(itemResppercent)[i]], interleave(itemResppercent[[i]]), "\\","\\", sep="")
  t4l<-rbind(t4l,mychar,deparse.level = 0)
}
write.table(t4l[1:5,],file=paste(texdir,"_t.freemarket.tex",sep="/"),quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(t4l[6:10,],file=paste(texdir,"_t.evo.tex",sep="/"),quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(t4l[11:15,],file=paste(texdir,"_t.cam.tex",sep="/"),quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(t4l[16:20,],file=paste(texdir,"_t.mwevo.tex",sep="/"),quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(t4l[21:25,],file=paste(texdir,"_t.mwnat.tex",sep="/"),quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(t4l[26:30,],file=paste(texdir,"_t.mwequ.tex",sep="/"),quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(t4l[31:35,],file=paste(texdir,"_t.relig.tex",sep="/"),quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(t4l[36:40,],file=paste(texdir,"_t.vax.tex",sep="/"),quote=FALSE,col.names=FALSE,row.names=FALSE)


#* --------------------- examine affirmation bias before reverse scoring -----------------------
if (correctaffirmbias) {
  #look at all item clusters that contain reverse-scored items, identify affirmation bias that way 
  #by just averaging all responses. This can be used to adjust responses of those item clusters
  revvars <- vn$qvarname[vn$reverse=="R"]
  revmean <- darwin1.5 %>% select(revvars) %>% apply(1,mean) 
  posmean <- darwin1.5 %>% select(matches("Q[C-G].[1-5]")) %>% select(-revvars) %>% apply(1,mean) 
  
  darwin1.5$affirmbias <- darwin1.5 %>% select(matches("Q[C-G].[1-5]")) %>% apply(1,mean)
  #use affirmation bias scores to adjust responses to the item clusters that do not contain reverse-scored items
  darwin1.5 %>% select(matches("Q[H-J].[1-5]")) %>% apply(1,mean) %>% plot(darwin1.5$affirmbias,.) 
  darwin1.5 <- darwin1.5 %>% mutate_at(vars(matches("Q[H-J].[1-5]")), regong ,darwin1.5$affirmbias)
  
} else {
  darwin1.5$affirmbias <- 0
}

#* --------------------- reverse score such that polarity is: -----------------------------------------------------
# QC Free market endorsement  "QC.2","QC.4","QC.5"
# QE accept evolution         "QE.4","QE.5"
# QF **reject** CAM           "QF.1","QF.3"  <--- note this used to point to acceptance
# QH men & women different evolved
# QI men & women different nature
# QJ men & women same
# QD religion                 "QD.4","QD.5"
# QG vaccinations             "QG.2","QG.4"
darwin2 <- darwin1.5 %>% mutate_at(vn$qvarname[vn$reverse=="R"],  revscore,mm=7)  %>% 
                         mutate_at(c("QB_1","QB_2","QB_3","QB_4"),revscore,mm=99) %>% #now reverse score the sliders 
                         mutate_at(vars(contains("QF")),revscore,mm=7) #flip CAM to point to rejection


#explore structure of sliders after reverse scoring
darwin2 %>% select(paste("QB_",as.character(1:4),sep="")) %>% cor(use="complete.obs")
darwin2 %>% select(paste("QB_",as.character(5:14),sep="")) %>% cor(use="complete.obs")
slider <- darwin2 %>% select(contains("QB_"))
cor(slider,use="complete.obs")
sliderpc <- prcomp(slider)
plot(sliderpc)
print(sliderpc)

# compute pairwise correlations within each cluster
for (cluster in c("C","E","F","H","I","J","D","G")) {
  darwin2 %>% select(num_range(paste("Q",cluster,".",sep=""),1:5)) %>% cor(use="complete.obs") %>% print
}
#summarize each cluster
for (cluster in c("C","E","F","H","I","J","D","G")) {
  darwin2 %>% select(num_range(paste("Q",cluster,".",sep=""),1:5)) %>% summary %>% print
}
# construct histogram for summary statistics
pdf(file=paste(figdir,"histoSocDarw.pdf",sep="/") ,height=9,width=6.5) 
par(mfrow=c(4,2))
clusterlabels <- c("Free market", "Evolution", "Rejection of CAM", "Men/women evolved differently", "Men/women naturally different",
                   "Men/women are the same", "Religiosity", "Vaccinations")
names(clusterlabels) <- c("C","E","F","H","I","J","D","G")
for (cluster in names(clusterlabels)) {
    darwin2 %>% select(num_range(paste("Q",cluster,".",sep=""),1:5)) %>% rowMeans %>% 
                hist(las=1,xlim=c(1,7),xlab="Average score",main=clusterlabels[cluster],col="light gray")
}
dev.off()
#construct histogram based on political quartiles at the end (after necessary intermediate variables have been computed)

#rescale sliders to 1-7 approximately to reduce variance differential with other items
rescale100 <- darwin2 %>% select(contains("QB_")) %>% mutate_all(funs(. /(100/6))) %>% mutate_all(funs(. +1))
darwin2[,names(rescale100)] <- rescale100

#*--------- compute measurement models for all constructs so we can use single-indicator models later --------------
freemarketvars <- paste("QC.",c(1:5),sep="")
fmmod <- c("fm =~ ",paste(freemarketvars,collapse=" + "),"QC.1 ~~ QC.3")
invisible(fmgof<-fitMM (fmmod,darwin2))

evovars <- paste("QE.",c(1:5),sep="")
evomod <- c("evo =~ ",paste(evovars,collapse=" + "),"QE.4 ~~ QE.5")
invisible(evogof<-fitMM(evomod,darwin2))

camvars <- paste("QF.",c(1:5),sep="")
cammod <- c("cam =~ ",paste(camvars,collapse=" + "),"QF.1 ~~ QF.3")
invisible(camgof<-fitMM(cammod,darwin2))

mwevovars <- paste("QH.",c(1:5),sep="") 
mwevomod <- c("mwevo =~ ",paste(mwevovars,collapse=" + "),"QH.2 ~~ QH.3")
invisible(mwevogof<-fitMM(mwevomod,darwin2))

mwnaturevars <- paste("QI.",c(1:5),sep="")
mwnaturemod <- c("mwnature =~ ",paste(mwnaturevars,collapse=" + "),"QI.2 ~~ QI.3")
invisible(mwnatgof<-fitMM(mwnaturemod,darwin2))

mwequalvars <- paste("QJ.",c(1:5),sep="")
mwequalmod <- c("mwequal =~ ",paste(mwequalvars,collapse=" + "),"QJ.2 ~~ QJ.5")
invisible(mwequgof<-fitMM(mwequalmod,darwin2))

religvars <- paste("QD.",c(1:5),sep="")
religmod <- c("relig =~ ",paste(religvars,collapse=" + "),"QD.2 ~~ QD.3")
invisible(religgof <- fitMM(religmod,darwin2))

vaxvars <- paste("QG.",c(1:5),sep="")
vaxmod <- c("vax =~ ",paste(vaxvars,collapse=" + "),"QG.2 ~~ QG.4")
invisible(vaxgof<-fitMM(vaxmod,darwin2))

## Deal with Conservatism scale 
slidervars <- names(slider) [-c(1:4)]
consertismmodel <- c("consertism =~ ",paste(slidervars,collapse=" + "),"QB_9 ~~ QB_10", "QB_6 ~~ QB_14" )
consertismgof <- fitMM(consertismmodel,darwin2)

sliderlibvars <- names(slider)[1:4]
liblismmodel <- c("liblism =~ ",paste(sliderlibvars,collapse=" + "))
liblismgof <- fitMM(liblismmodel,darwin2)

#two-factor model shows no correlation (or negative) between liberal and conservative items
conslib.two <- c(consertismmodel, "\n", 
                 liblismmodel, "QB_7 ~~  QB_1")
conslibfit <- sem(conslib.two,darwin2)
summary(conslibfit, standardized=TRUE, fit.measures=TRUE)
conslibfitgof <- fitmeasures(conslibfit)
modificationIndices(conslibfit, sort. = TRUE, maximum.number = 4)
semPaths(conslibfit, "std", title =FALSE, curvePivot = TRUE)


#* ---------- compute single-indicators models -------------------------------------------------------------------
fmSI      <- singleindmodel(freemarketvars,list(c("QC.1","QC.3")),darwin2)
evoSI     <- singleindmodel(evovars, list(c("QE.4","QE.5")),darwin2)
camSI     <- singleindmodel(camvars, list(c("QF.1","QF.3")),darwin2)
mwevoSI   <- singleindmodel(mwevovars, list(c("QH.2","QH.3")),darwin2) 
mwnatSI   <- singleindmodel(mwnaturevars, list(c("QI.2","QI.3")),darwin2)
mwequSI   <- singleindmodel(mwequalvars, list(c("QJ.2","QJ.5")),darwin2)
religSI   <- singleindmodel(religvars, list(c("QD.2","QD.3")),darwin2)
vaxSI     <- singleindmodel(vaxvars, list(c("QG.2","QG.4")),darwin2)
consSI    <- singleindmodel(slidervars, list(c("QB_9","QB_10"), c("QB_6","QB_14")),darwin2)

compositedarwin <- data.frame ( #compute composite scores for the SI models
                      fm =   apply(darwin2[,freemarketvars], 1,mean),
                      evo =  apply(darwin2[,evovars], 1,mean),
                      cam =  apply(darwin2[,camvars], 1,mean),
                      mwevo =   apply(darwin2[,mwevovars], 1,mean),
                      mwnat =   apply(darwin2[,mwnaturevars], 1,mean),
                      mwequ = apply(darwin2[,mwequalvars], 1,mean),
                      relig = apply(darwin2[,religvars], 1,mean),
                      vax = apply(darwin2[,vaxvars], 1,mean),
                      cons = apply(darwin2[,slidervars], 1,mean),
                      gender = darwin2$Qgender,
                      age = darwin2$Qage
                    )
cor(compositedarwin)
histogram(~cam|gender,data=compositedarwin)
aggregate(cam~gender,data=compositedarwin,FUN=mean)

#*------------- correlation structure among latent constructs -----------------------------------------------------
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

# and now the full model ...  
modelCorrel <- c("
                  fmFac    =~ fm
                 evoFac   =~ evo
                 camFac   =~ cam
                 mwevoFac =~ mwevo
                 mwnatFac =~ mwnat
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
                 "relig ~~", religSI$eSImod, "*relig",
                 "cons ~~",  consSI$eSImod,  "*cons",
                 "vax ~~",   vaxSI$eSImod,   "*vax"
)
 fitCorrel <- sem(modelCorrel, compositedarwin, std.lv=TRUE, estimator="ML") # http://lavaan.ugent.be/tutorial/groups.html
 summary(fitCorrel,standardized=TRUE, fit.measures=TRUE)
 x11()
 semPaths(fitCorrel, what="std", title =FALSE, curvePivot = TRUE,residuals=FALSE, structural=TRUE, layout="circle")
 lavInspect(fitCorrel, what = "cor.lv")
 eigen( inspect(fitCorrel, "cov.lv") )$values 
 
 #now compute p-values
 pvals2tailed <- pnorm(abs(inspect(fitCorrel,what="est")$psi/inspect(fitCorrel,what="se")$psi),lower.tail = FALSE)*2
 colnames(pvals2tailed) <- rownames(pvals2tailed) <- c(clusterlabels,"Conservatism")
 pvals2tailed[upper.tri(pvals2tailed)] <- 0
 
  #get matrix printed
 lvcormat <- lavInspect(fitCorrel, what = "cor.lv")
 lvcormat[upper.tri(lvcormat,diag=TRUE)] <- NA
 colnames(lvcormat) <- rownames(lvcormat) <- c(clusterlabels,"Conservatism")
 cormat <- stargazer(lvcormat, title="Correlations among all latent variables") 
 for (i in 1:9) {
   lvcp<-i+11
   cormat[lvcp] <- str_replace_all(cormat[lvcp],fixed("$-$"),"-")
   sigs<-unique(as.numeric(round(lvcormat[i,which(pvals2tailed[i,]>.05)],3)))
   if (length(sigs)>0) {
     for (j in 1:length(sigs)){
       subss<-substr(paste("$",as.character(sigs[j]),"0000",sep=""),1,ifelse(sigs[j]<0,7,6))
       w2r <- (str_locate_all(cormat[lvcp],fixed(subss)))[[1]]
       if (dim(w2r)[1]>0) {
         for (k in 1:dim(w2r)[1]) {
           str_sub(cormat[lvcp],w2r[k,1],w2r[k,2])<-paste(subss,"ns",sep="") #"^{*}"
         }#print(c(substr(cormat[lvcp],x[1],x[2]),subss)))
       }
     }
   }
 }
 write.table(cormat[13:20],file=paste(texdir,"_t.lvcor.tex",sep="/"),quote=FALSE,col.names=FALSE,row.names=FALSE) 
 
 #now look at effects of gender: unconstrained first
 fitCorrelFM <- sem(modelCorrel, compositedarwin, 
                                      std.lv=TRUE, estimator="ML",
                                      group="gender")
 summary(fitCorrelFM,standardized=TRUE, fit.measures=TRUE)
 fitCorrelFM2 <- sem(modelCorrel, compositedarwin, 
                                      std.lv=TRUE, estimator="ML",
                                      group="gender",
                                      group.equal = c("lv.covariances"))
 summary(fitCorrelFM2,standardized=TRUE, fit.measures=TRUE)
 aovresult <- anova(fitCorrelFM,fitCorrelFM2)
 equcovgof <- fitMeasures(fitCorrelFM2)

 
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
 

 modelVaxEvoCAM <- c("  vaxFac   ~ fmFac +  mwequFac + consFac 
                        evoFac   ~ fmFac +  religFac + mwequFac + mwevoFac  + consFac 
                        camFac   ~ fmFac +  religFac 

                         mwevoFac ~~ 0*consFac
                         evoFac   ~~ 0*camFac

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
 fitVaxEvoCam <- sem(modelVaxEvoCAM, compositedarwin, std.lv=TRUE, estimator="ML")
 summary(fitVaxEvoCam,standardized=TRUE, fit.measures=TRUE)
 
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
 
 
 
## ---- Deal with CRT --------------------------------------------------------------------------------
 crt1 <- darwin2$CRT.1_1 %>% as.character %>% str_to_lower %>% str_replace(.,"cents","") %>% str_replace(.,"$","") %>% as.numeric
 crt1 <- sapply(crt1,FUN=function(x) ifelse(x>=5,x/100,x)) #convert from cents to dollars where applicable
 
 crt2 <- darwin2$CRT.2_1 %>% as.character %>% str_to_lower %>% str_replace(.,"minutes","") %>%  as.numeric
 
 crt3 <- darwin2$CRT.3_1 %>% as.character %>% str_to_lower %>% str_replace(.,"days","") %>% as.numeric
 
 darwin2$crtScore <- apply(cbind(crt1==.05,crt2==5,crt3==47),1,sum,na.rm=TRUE)
 crtnums <- table(darwin2$crtScore)
 crtperc <- round(crtnums/sum(crtnums)*100,2)
 crtmn   <- mean(darwin2$crtScore) 
 
 
## ---- Create mean-indicators for all clusters for correlation with CRT -------
meanIndicators <- NULL
for (cluster in c("C","E","F","H","I","J","D","G")) {
  meanIndicators <- cbind(meanIndicators,
                      darwin2 %>% select(num_range(paste("Q",cluster,".",sep=""),1:5)) %>% 
                        apply(.,1,mean) 
  )
}

# include all of conservatism items, including those with liberal polarity
consmean <- darwin2 %>% select(num_range("QB_",1:14)) %>% apply(.,1,mean) 
meanIndicators <- as.data.frame(cbind(meanIndicators,consmean,darwin2$crtScore))
# overall correlations, to parallel those at latent-variable level
compcormat <- rcorr(as.matrix(meanIndicators), type="pearson")$r
compcormat[upper.tri(compcormat,diag=TRUE)] <- NA
pcomposite <- rcorr(as.matrix(meanIndicators), type="pearson")$P
pcomposite[upper.tri(pcomposite,diag=TRUE)] <- NA

#get matrix printed
colnames(compcormat) <- rownames(compcormat) <- c(clusterlabels,"Conservatism","CRT")
cormat4p <- stargazer(compcormat, title="Correlations among composite indicators") 
for (i in 1:10) {
  lvcp<-i+11
  cormat4p[lvcp] <- str_replace_all(cormat4p[lvcp],fixed("$-$"),"-")
  sigs<-unique(as.numeric(round(compcormat[i,which(pcomposite[i,]>.05)],3)))
  if (length(sigs)>0) {
    for (j in 1:length(sigs)){
      subss<-substr(paste("$",as.character(sigs[j]),"0000",sep=""),1,ifelse(sigs[j]<0,7,6))
      w2r <- (str_locate_all(cormat4p[lvcp],fixed(subss)))[[1]]
      if (dim(w2r)[1]>0) {
        for (k in 1:dim(w2r)[1]) {
          str_sub(cormat4p[lvcp],w2r[k,1],w2r[k,2])<-paste(subss,"ns",sep="") #"^{*}"
        }#print(c(substr(cormat4p[lvcp],x[1],x[2]),subss)))
      }
    }
  }
}
cormat4p[12:21] <- str_replace_all(cormat4p[12:21],fixed("& $$ \\"),"\\")
write.table(cormat4p[13:20],file=paste(texdir,"_t.compcor.tex",sep="/"),quote=FALSE,col.names=FALSE,row.names=FALSE) 
write.table(cormat4p[21],file=paste(texdir,"_t.crtcor.tex",sep="/"),quote=FALSE,col.names=FALSE,row.names=FALSE) 

#update mean indicators with demographics and rename variables to something less decorative
meanIndicators <- as.data.frame(cbind(meanIndicators,darwin2$Qgender,darwin2$Qage,darwin2$affirmbias))
names(meanIndicators) <- c("FMmean","Evomean","CAMmean","MWevomean","MWnatmean","MWequmean",
                           "Religmean","Vaxmean","consmean","crtscore","gender","age","affirmbias")
# median split on conservatism 
meanIndicators$medsplcons <- ifelse(meanIndicators$consmean < median(meanIndicators$consmean), 0,1)
histogram(~MWequmean|medsplcons,data=meanIndicators)
meanIndicators %>% group_by(medsplcons) %>% summarise(mean(MWnatmean))
meanIndicators %>% filter(medsplcons==0) %>% select(starts_with("MW")) %>% cor
meanIndicators %>% filter(medsplcons==1) %>% select(starts_with("MW")) %>% cor

#chasing conservative hypocrisy
ggpanels <- list()
lowerab <- 0   #if these are set to 0 and 1, then all obs are used. Otherwise eliminate extreme affirmation bias folks
upperab <- 1
k <- 0
for (qqs in c(.25,.5)) {
  topcons <- meanIndicators %>% filter(affirmbias >= quantile(affirmbias,lowerab) & affirmbias <= quantile(affirmbias,upperab))  %>% 
    filter(consmean > quantile(consmean,1-qqs))  
  toplibs <- meanIndicators %>% filter(affirmbias >= quantile(affirmbias,lowerab) & affirmbias <= quantile(affirmbias,upperab))  %>% 
    filter(consmean < quantile(consmean,qqs)) 
  k<-k+1
  ggpanels[[k]] <- plotexploration(toplibs, toplibs$MWevomean, toplibs$MWnatmean, toplibs$Evomean,
                  "Men and women evolved differently",
                  "Men and women are naturally different",
                  "Evolution","YlOrBr",-1,
                  paste("Liberals (top ",qqs*100,"%)",sep=""))
  k<-k+1
  ggpanels[[k]] <- plotexploration(topcons, topcons$MWevomean, topcons$MWnatmean, topcons$Evomean,
                                   "Men and women evolved differently",
                                   "Men and women are naturally different",
                                   "Evolution", "YlOrBr",-1,
                                   paste("Conservatives (top ",qqs*100,"%)",sep=""))
}
x11(width=12,height=8)
gridplt <- do.call("grid.arrange", c(ggpanels, nrow=2))
if (lowerab==0 & upperab==1) {
  ggfn <- paste(figdir,"conslibevogrid.pdf",sep="/")  
} else {
  ggfn <- paste(figdir,paste("conslibevogrid",lowerab,upperab,".pdf",sep=""),sep="/")  
}
ggsave(ggfn, plot=gridplt, width=12, height=8)



#chasing liberal hypocrisy
ggpanels2 <- list()
lowerab <- .0   #if these are set to 0 and 1, then all obs are used. Otherwise eliminate extreme affirmation bias folks
upperab <- 1.
k <- 0
for (qqs in c(.25,.5)) {
  topcons <- meanIndicators %>% filter(affirmbias >= quantile(affirmbias,lowerab) & affirmbias <= quantile(affirmbias,upperab))  %>% 
                                filter(consmean > quantile(consmean,1-qqs))  
  toplibs <- meanIndicators %>% filter(affirmbias >= quantile(affirmbias,lowerab) & affirmbias <= quantile(affirmbias,upperab))  %>% 
                                filter(consmean < quantile(consmean,qqs)) 
  k<-k+1
  ggpanels2[[k]] <- plotexploration(toplibs, toplibs$MWevomean,toplibs$MWnatmean,  toplibs$MWequmean,
                                    "Men and women evolved differently",  
                                   "Men and women naturally different",
                                   "MW Same","YlGnBu",1,
                                   paste("Liberals (top ",qqs*100,"%)",sep=""))
  k<-k+1
  ggpanels2[[k]] <- plotexploration(topcons, topcons$MWevomean,topcons$MWnatmean,  topcons$MWequmean, 
                                    "Men and women evolved differently",
                                   "Men and women naturally different",
                                   "MW Same","YlGnBu",1,
                                   paste("Conservatives (top ",qqs*100,"%)",sep=""))
}
x11(width=12,height=8)
gridplt2 <- do.call("grid.arrange", c(ggpanels2, nrow=2))
if (lowerab==0 & upperab==1) {
    ggfn <- paste(figdir,"conslibequalitygrid.pdf",sep="/") 
} else {
    ggfn <- paste(figdir,paste("conslibequalitygrid",lowerab,upperab,".pdf",sep=""),sep="/") 
}
ggsave(ggfn, plot=gridplt2, width=12, height=8)



