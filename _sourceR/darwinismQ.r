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
library(DiagrammeR)
library(webshot)
library(htmltools)

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
  revmean <- NA   #so LaTex compilation does not crash
  posmean <- NA
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


## ---- Deal with CRT --------------------------------------------------------------------------------
crt1 <- darwin2$CRT.1_1 %>% as.character %>% str_to_lower %>% str_replace(.,"cents","") %>% str_replace(.,"$","") %>% as.numeric
crt1 <- sapply(crt1,FUN=function(x) ifelse(x>=5,x/100,x)) #convert from cents to dollars where applicable

crt2 <- darwin2$CRT.2_1 %>% as.character %>% str_to_lower %>% str_replace(.,"minutes","") %>% str_replace(.,"min","") %>%  as.numeric

crt3 <- darwin2$CRT.3_1 %>% as.character %>% str_to_lower %>% str_replace(.,"days","") %>% as.numeric

darwin2$crtScore <- apply(cbind(crt1==.05,crt2==5,crt3==47),1,sum,na.rm=TRUE)
crtnums <- table(darwin2$crtScore)
crtperc <- round(crtnums/sum(crtnums)*100,2)
crtmn   <- mean(darwin2$crtScore) 

#check if CRT can be added to SEM by using measurement model with ordered variables
tmp.data <- data.frame(crt1==.05,crt2==5,crt3==47)
names(tmp.data) <- c("crt1", "crt2", "crt3") 
compositedarwin <- data.frame(compositedarwin,tmp.data)

crtmodel <- c("crtFac =~ a*crt1 + a*crt2 + a*crt3")
crtmodfit <- sem(crtmodel, data=compositedarwin, missing="pairwise", ordered = c("crt1", "crt2", "crt3"), std.lv=TRUE) 
summary(crtmodfit, standardized=TRUE, fit.measures=TRUE)
crtmodgof <- fitMeasures(crtmodfit)



#*------------- correlation structure among latent constructs -----------------------------------------------------
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
                 crtFac =~ a*crt1 + a*crt2 + a*crt3
                 
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
colnames(pvals2tailed) <- rownames(pvals2tailed) <- c(clusterlabels,"Conservatism","CRT")
pvals2tailed[upper.tri(pvals2tailed)] <- 0

#get matrix printed
lvcormat <- lavInspect(fitCorrel, what = "cor.lv")
lvcormat[upper.tri(lvcormat,diag=TRUE)] <- NA
colnames(lvcormat) <- rownames(lvcormat) <- c(clusterlabels,"Conservatism","CRT")
cormat <- stargazer(lvcormat, title="Correlations among all latent variables") 
for (i in 1:10) {
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
write.table(cormat[13:21],file=paste(texdir,"_t.lvcor.tex",sep="/"),quote=FALSE,col.names=FALSE,row.names=FALSE) 


#*----------------------------- now look at effects of gender: unconstrained first -------------------------
fitCorrelFM <- sem(modelCorrel, compositedarwin, 
                   std.lv=TRUE, estimator="ML",
                   group="gender")
summary(fitCorrelFM,standardized=TRUE, fit.measures=TRUE)
#now constrain all
fitCorrelFM2 <- sem(modelCorrel, compositedarwin, 
                    std.lv=TRUE, estimator="ML",
                    group="gender",
                    group.equal = c("lv.covariances"))
summary(fitCorrelFM2,standardized=TRUE, fit.measures=TRUE)

#and finally partial constraints embodied in model
modelCorrelcons <- c("
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
                  crtFac =~ a*crt1 + a*crt2 + a*crt3
                 
                 mwevoFac ~~ c(v1,v1)*mwnatFac
                 mwevoFac ~~ c(v2,v2)*mwequFac
                 mwnatFac ~~ c(v3,v3)*mwequFac
                 
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
fitCorrelFM3 <- sem(modelCorrelcons, compositedarwin, 
                    std.lv=TRUE, estimator="ML",
                    group="gender")
summary(fitCorrelFM3,standardized=TRUE, fit.measures=TRUE)

aovresult <- anova(fitCorrelFM,fitCorrelFM3,fitCorrelFM2)
equcovgof <- fitMeasures(fitCorrelFM2)
equpartcovgof <- fitMeasures(fitCorrelFM3)



#------Beast route to modeling: first zero in on politics alone, then add other components.
polmodelVaxEvoCAM <- c("vaxFac   ~ fmFac + 0*allPolFac
                         evoFac   ~ religFac + 0*allPolFac 
                         camFac   ~ allPolFac 
                        
                         allPolFac =~ fmFac + consFac + religFac
                        
                         evoFac ~~ 0*camFac

                         fmFac    =~ fm
                         evoFac   =~ evo
                         camFac   =~ cam
                         religFac =~ relig
                         vaxFac   =~ vax
                         consFac  =~ cons
                         
                         fm ~~",     fmSI$eSImod,    "*fm",
                       "evo ~~ ",  evoSI$eSImod,   "*evo",
                       "cam ~~",   camSI$eSImod,   "*cam",
                       
                       "cons ~~",  consSI$eSImod, "*cons",
                       "relig ~~", religSI$eSImod, "*relig",
                       "vax ~~",   vaxSI$eSImod,    "*vax"
)
fitPolVaxEvoCam <- sem(polmodelVaxEvoCAM, compositedarwin, std.lv=TRUE, estimator="ML")
summary(fitPolVaxEvoCam,standardized=TRUE, fit.measures=TRUE) 
finalpolmodelgof <- fitMeasures(fitPolVaxEvoCam)

fitPolVaxEvoCamfree <- sem(str_replace_all(polmodelVaxEvoCAM,"0\\*a","a"), compositedarwin, std.lv=TRUE, estimator="ML")
summary(fitPolVaxEvoCamfree,standardized=TRUE, fit.measures=TRUE) 

aovpol <- anova(fitPolVaxEvoCamfree,fitPolVaxEvoCam)
aovpolfindf<- aovpol[["Df"]][2]- aovpol[["Df"]][1]


# Extract to-be-graphed elements from lavaan object: based on https://rpubs.com/tjmahr/sem_diagrammer
t4modP <- "
 digraph {

graph [layout = neato,
       overlap = scale, splines=true]

node [shape = circle ]

allPolFac  [pos = '0,0!', label = ' All \\n politics ', shape = circle, fontsize=14 width=1.,fixedsize=true]
fmFac  [pos = '-4,1.2!', label = ' Free \\n market ', shape = circle, fontsize=14 width=1.,fixedsize=true]
evoFac  [pos = '4,-1.2!', label = ' Evolution ', shape = circle, fontsize=14 width=1.,fixedsize=true]
camFac  [pos = '4,0!', label = ' Reject \\n CAM ', shape = circle, fontsize=14 width=1.,fixedsize=true]

C [pos='5.2,0!' label = 'C', style = invis width=.01 fixedsize=true]
D [pos='5.2,0.8!' label = 'D', style = invis width=.01 fixedsize=true]

religFac  [pos = '-4,-1.2!', label = ' Religiosity ', shape = circle, fontsize=14 width=1.,fixedsize=true]
vaxFac  [pos = '4,1.2!', label = ' Vax ', shape = circle, fontsize=14 width=1.,fixedsize=true]
consFac  [pos = '-4,0!', label = ' Conser- \\n vatism ', shape = circle, fontsize=14 width=1.,fixedsize=true]
"
pathsP <- fitPolVaxEvoCam %>%
  parameterestimates(.,standardized=TRUE) %>%
  select(lhs, op, rhs, std.all)

# Latent variables are left-hand side of "=~" lines
latentP <- pathsP %>%
  filter(op == "=~") %>%
  select(nodes = lhs) %>%
  distinct 
# Edges will be labeled by the parameter estimates
all_pathsP <- pathsP %>%
  filter(op != "~1") %>%
  mutate(label = round(std.all, 2)) %>%
  select(-std.all)

# Factor loadings are the paths in the "=~" lines, drop all but second-order factor
loadingsP <- all_pathsP %>%
  filter(op == "=~") %>% filter(lhs == "allPolFac") %>%
  rename(to = rhs, from = lhs) %>%
  select(from, to, label)
t4modP <- c(t4modP,apply(loadingsP,1,FUN=function(x) (paste(x["from"]," -> ",x["to"], " [label = '", x["label"], "'  fontsize=18]"))))

regressionsP <- all_pathsP %>%
  filter(op == "~") %>%
  rename(to = lhs, from = rhs) %>% filter(label!=0) %>%
  select(from, to, label)
t4modP<-c(t4modP,apply(regressionsP,1,FUN=function(x) (paste(x["from"]," -> ",x["to"], " [label = '", x["label"], "'  fontsize=18]"))))

#extract covariances 
covarsP <- all_pathsP %>%
  filter(op == "~~") %>%
  rename(to = rhs, from = lhs) %>%
  select(from, to, label) %>% filter(from != to) %>% filter(label!=0)
t4modP<-c(t4modP, paste("C:se -> ", covarsP[1,"from"],":e [splines=curved  style = 'dashed' label = '", covarsP[1,"label"], "'  fontsize=18] "))
t4modP<-c(t4modP, paste("C:ne -> ", covarsP[1,"to"],":e [splines=curved  style = 'dashed'] "))

t4modP<-c(t4modP, paste("D:se -> ", covarsP[2,"from"],":e [splines=curved  style = 'dashed'] "))
t4modP<-c(t4modP, paste("D:ne -> ", covarsP[2,"to"],":e [splines=curved  style = 'dashed' label = '", covarsP[2,"label"], "'  fontsize=18] "))
t4modP <- c(t4modP, "}")

writeLines(t4modP,con="finalpolSEM.txt") 
polhtml <- grViz("finalpolSEM.txt")
html_print(add_mathjax(polhtml)) %>% webshot(file = paste(figdir,"/finalpolSEM.pdf",sep="/"),cliprect = "viewport",zoom=.7)




#---- now add CRT to get super model -----------------------------------------------------
supermodel <- c("   vaxFac   ~ fmFac +  mwequFac + crtFac 
                     evoFac   ~ religFac + mwevoFac  + mwnatFac + crtFac
                     camFac   ~ allPolFac + crtFac
                            
                     
                       evoFac ~~ 0* camFac
                        mwevoFac ~~ 0* allPolFac
                        crtFac ~~ 0*mwnatFac
                        crtFac ~~ 0*mwevoFac
               
                      crtFac =~ a*crt1 + a*crt2 + a*crt3
                      allPolFac =~ fmFac + consFac + religFac
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
supermodelfit <- sem(supermodel, data=compositedarwin, missing="pairwise", ordered = c("crt1", "crt2", "crt3"), std.lv=TRUE) 
summary(supermodelfit, standardized=TRUE, fit.measures=TRUE)
supermodelgof <- fitMeasures(supermodelfit)
modificationIndices(supermodelfit, sort. = TRUE, maximum.number = 6)


# Extract to-be-graphed elements from lavaan object: based on https://rpubs.com/tjmahr/sem_diagrammer
t4mod <- "
 digraph {

graph [layout = neato,
       overlap = true,
       outputorder = edgesfirst]

node [shape = circle ]

crtFac     [pos = '0,3!', label = ' CRT ', shape = circle, fontsize=14 width=1.,fixedsize=true]
allPolFac  [pos = '0,0!', label = ' All \\n politics ', shape = circle, fontsize=14 width=1.,fixedsize=true]
mwnatFac  [pos = '-0.5,-2.5!', label = ' M&W \\n NatDiff ', shape = circle, fontsize=14 width=1.,fixedsize=true]
fmFac  [pos = '-4,1.2!', label = ' Free \\n market ', shape = circle, fontsize=14 width=1.,fixedsize=true]
evoFac  [pos = '4,-1.2!', label = ' Evolution ', shape = circle, fontsize=14 width=1.,fixedsize=true]
camFac  [pos = '4,0!', label = ' Reject \\n CAM ', shape = circle, fontsize=14 width=1.,fixedsize=true]
mwevoFac  [pos = '1,-3!', label = ' M&W \\n Evo ', shape = circle, fontsize=14 width=1.,fixedsize=true]
mwequFac  [pos = '-2,-2!', label = ' M&W \\n Equal ', shape = circle, fontsize=14 width=1.,fixedsize=true]
religFac  [pos = '-4,-1.2!', label = ' Religiosity ', shape = circle, fontsize=14 width=1.,fixedsize=true]
vaxFac  [pos = '4,1.2!', label = ' Vax ', shape = circle, fontsize=14 width=1.,fixedsize=true]
consFac  [pos = '-4,0!', label = ' Conser- \\n vatism ', shape = circle, fontsize=14 width=1.,fixedsize=true]
"
paths <- supermodelfit %>%
  parameterestimates(.,standardized=TRUE) %>%
  select(lhs, op, rhs, std.all)

# Latent variables are left-hand side of "=~" lines
latent <- paths %>%
  filter(op == "=~") %>%
  select(nodes = lhs) %>%
  distinct 
# Edges will be labeled by the parameter estimates
all_paths <- paths %>%
  filter(op != "~1") %>%
  mutate(label = round(std.all, 2)) %>%
  select(-std.all)

# Factor loadings are the paths in the "=~" lines, drop all but second-order factor
loadings <- all_paths %>%
  filter(op == "=~") %>% filter(lhs == "allPolFac") %>%
  rename(to = rhs, from = lhs) %>%
  select(from, to, label)
t4mod <- c(t4mod,apply(loadings,1,FUN=function(x) (paste(x["from"]," -> ",x["to"], " [label = '", x["label"], "'  fontsize=18]"))))

regressions <- all_paths %>%
  filter(op == "~") %>%
  rename(to = lhs, from = rhs) %>%
  select(from, to, label)

t4mod<-c(t4mod,apply(regressions,1,FUN=function(x) (paste(x["from"]," -> ",x["to"], " [label = '", x["label"], "'  fontsize=18]"))))
t4mod <- c(t4mod, "}")
writeLines(t4mod,con="supermodelViz.txt") 
superhtml <- grViz("supermodelViz.txt")
html_print(add_mathjax(superhtml)) %>% webshot(file = paste(figdir,"/superSEM.pdf",sep="/"),cliprect = "viewport",zoom=0.5)


#extract covariances for separate table
covars <- all_paths %>%
  filter(op == "~~") %>%
  rename(to = rhs, from = lhs) %>%
  select(from, to, label) %>% filter(from != to) %>% filter(label!=0)



##------ use politics to predict gender attitudes ------------------------------
polmodelgender <- c("mwnatFac   ~  allPolFac 
                         mwevoFac   ~ allPolFac + religFac
                         mwequFac   ~ allPolFac 
                         
                         allPolFac =~ fmFac + consFac + religFac
                         
                      mwnatFac =~ mwnat
                      mwevoFac =~ mwevo
                      mwequFac =~ mwequ  
                         
                         fmFac    =~ fm
                         religFac =~ relig
                         consFac  =~ cons
                         
                         fm ~~",     fmSI$eSImod,    "*fm",
                       "cons ~~",  consSI$eSImod, "*cons",
                       "relig ~~", religSI$eSImod, "*relig",
                       
                       "mwevo ~~", mwevoSI$eSImod, "*mwevo",
                       "mwnat ~~", mwnatSI$eSImod, "*mwnat",
                       "mwequ ~~", mwequSI$eSImod, "*mwequ"
                       )
fitpolmodelgender <- sem(polmodelgender, compositedarwin, std.lv=TRUE, estimator="ML")
summary(fitpolmodelgender,standardized=TRUE, fit.measures=TRUE) 
gofpolmodelgender <- fitMeasures(fitpolmodelgender)
modificationIndices(fitpolmodelgender, sort. = TRUE, maximum.number = 6)

# Extract to-be-graphed elements from lavaan object: based on https://rpubs.com/tjmahr/sem_diagrammer
t4modG <- "
 digraph {

graph [layout = neato,
       overlap = scale, splines=true]

node [shape = circle ]

allPolFac  [pos = '0,0!', label = ' All \\n politics ', shape = circle, fontsize=14 width=1.,fixedsize=true]
fmFac  [pos = '-4,1.2!', label = ' Free \\n market ', shape = circle, fontsize=14 width=1.,fixedsize=true]
mwevoFac  [pos = '4,-1.2!', label = ' M&W \n Evo ', shape = circle, fontsize=14 width=1.,fixedsize=true]
mwnatFac  [pos = '4,0!', label = ' M&W \\n NatDiff ', shape = circle, fontsize=14 width=1.,fixedsize=true]

C [pos='4.8,-0.8!' label = 'C', style = invis width=.01 fixedsize=true]
D [pos='4.8,0.8!' label = 'D', style = invis width=.01 fixedsize=true]
E [pos='5.8,0.!' label = 'D', style = invis width=.01 fixedsize=true]

religFac  [pos = '-4,-1.2!', label = ' Religiosity ', shape = circle, fontsize=14 width=1.,fixedsize=true]
mwequFac  [pos = '4,1.2!', label = ' M&W \n Equal ', shape = circle, fontsize=14 width=1.,fixedsize=true]
consFac  [pos = '-4,0!', label = ' Conser- \\n vatism ', shape = circle, fontsize=14 width=1.,fixedsize=true]
"
pathsG <- fitpolmodelgender %>%
  parameterestimates(.,standardized=TRUE) %>%
  select(lhs, op, rhs, std.all)

# Latent variables are left-hand side of "=~" lines
latentG <- pathsG %>%
  filter(op == "=~") %>%
  select(nodes = lhs) %>%
  distinct 
# Edges will be labeled by the parameter estimates
all_pathsG <- pathsG %>%
  filter(op != "~1") %>%
  mutate(label = round(std.all, 2)) %>%
  select(-std.all)

# Factor loadings are the paths in the "=~" lines, drop all but second-order factor
loadingsG <- all_pathsG %>%
  filter(op == "=~") %>% filter(lhs == "allPolFac") %>%
  rename(to = rhs, from = lhs) %>%
  select(from, to, label)
t4modG <- c(t4modG,apply(loadingsG,1,FUN=function(x) (paste(x["from"]," -> ",x["to"], " [label = '", x["label"], "'  fontsize=18]"))))

regressionsG <- all_pathsG %>%
  filter(op == "~") %>%
  rename(to = lhs, from = rhs) %>% filter(label!=0) %>%
  select(from, to, label)
t4modG<-c(t4modG,apply(regressionsG,1,FUN=function(x) (paste(x["from"]," -> ",x["to"], " [label = '", x["label"], "'  fontsize=18]"))))

#extract covariances 
covarsG <- all_pathsG %>%
  filter(op == "~~") %>%
  rename(to = rhs, from = lhs) %>%
  select(from, to, label) %>% filter(from != to) %>% filter(label!=0)

t4modG<-c(t4modG, paste("C:se -> ", covarsG[1,"from"],":se [splines=curved  style = 'dashed' label = '", covarsG[1,"label"], "'  fontsize=18] "))
t4modG<-c(t4modG, paste("C:ne -> ", covarsG[1,"to"],":e [splines=curved  style = 'dashed'] "))

t4modG<-c(t4modG, paste("D:se -> ", covarsG[2,"from"],":ne [splines=curved  style = 'dashed'  label = '", covarsG[2,"label"], "'  fontsize=18] "))
t4modG<-c(t4modG, paste("D:ne -> ", covarsG[2,"to"],":e [splines=curved  style = 'dashed'] "))

t4modG<-c(t4modG, paste("E:se -> ", covarsG[3,"from"],":e [splines=curved  style = 'dashed'  label = '", covarsG[3,"label"], "'  fontsize=18] "))
t4modG<-c(t4modG, paste("E:ne -> ", covarsG[3,"to"],":e [splines=curved  style = 'dashed'] "))

t4modG <- c(t4modG, "}")

writeLines(t4modG,con="finalgenderSEM.txt") 
polgenderhtml <- grViz("finalgenderSEM.txt")
html_print(add_mathjax(polgenderhtml)) %>% webshot(file = paste(figdir,"/finalgenderSEM.pdf",sep="/"),cliprect = "viewport",zoom=.7)



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

#now combine both figures as one using only the top/bottom 50%
x11(width=12,height=8)
ggpanels3 <- list(ggpanels[[3]],ggpanels[[4]],ggpanels2[[3]],ggpanels2[[4]])
gridplt3 <- do.call("grid.arrange", c(ggpanels3, nrow=2))
if (lowerab==0 & upperab==1) {
  ggfn <- paste(figdir,"conslibequalitygridAll.pdf",sep="/") 
} else {
  ggfn <- paste(figdir,paste("conslibequalitygridAll",lowerab,upperab,".pdf",sep=""),sep="/") 
}
ggsave(ggfn, plot=gridplt3, width=12, height=8)



