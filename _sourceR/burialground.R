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
