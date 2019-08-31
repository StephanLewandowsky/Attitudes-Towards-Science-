# extract items for table in paper
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

mklatexline<- function(rowof2) {
   paste("\\emph{",rowof2["shortname"],"} & ", rowof2["longname"],"\\","\\",sep="")
}

wd <- "C:/Users/Lewan/Documents/Research Projects/Jan Woike social darwinism and evolution/Q full study/data from Q"
vn <- read.csv(paste(wd,"variableNames.csv",sep="/"),header=TRUE,row.names=NULL,stringsAsFactors = FALSE) 

fortab <- str_detect(vn$qvarname,fixed(".")) & !str_detect(vn$qvarname,"CRT")
vnfortab <- subset(vn,fortab) %>% select(c(shortname,longname))

ltx <- as.character(apply(vnfortab,1,mklatexline))
write.table(ltx,file="_t.items.tex",quote=FALSE,col.names=FALSE,row.names=FALSE)
