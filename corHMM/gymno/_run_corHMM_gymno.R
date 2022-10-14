##############################################################################################
##############################################################################################
##############################################################################################
### corHMM

# This script was used to model trait evolution using the corHMM package, for the Gymnophiona
# dataset of Liedtke et al. 2022.
##############################################################################################
##############################################################################################
##############################################################################################

##############################################################################################
# set working directory and load required libraries

setwd("~/Documents/amphibian_diversity_project/2021/_SI_scripts/corHMM/gymno/")

library(ape)
library(geiger)
library(phytools)
library(corHMM)
library(tidyverse)
library(parallel)
library(MuMIn)
library(qgraph)



##############################################################################################
# load MCC tree and data
phy<-read.tree("gymno.tre")
dat<-read_csv("gymno.csv")

# check tree and data match
name.check(phy, data.names = dat$species_asw)

# sort data to match tree
dat<-as.data.frame(dat[,c("species_asw","rep_mode1")])

# summarize reproductive mode
summary(as.factor(dat$rep_mode1))

# where:
# D = Direct development
# S = Semi-terrestrial
# V = Viviparity/Live-bearing (L in manuscript)

##############################################################################################
# call on auxiliary script to build transition matrices

pdf("corHMM_gymno_q.pdf")
source("_corHMM_models_gymno.R")
dev.off()

##############################################################################################
## code unknowns as uncertain

dat$rep_mode1[is.na(dat$rep_mode1)]<-"D&S&V"
summary(as.factor(dat$rep_mode1))

##############################################################################################
## append shared arguments to models

q<-lapply(q, FUN=append,
          list(phy=phy,
               data=dat,
               node.states="none",
               get.tip.states=F,
               nstarts = 25,
               n.cores = 25))## warning!! high number of cores


## remove the 'legend' term to allow do.call to run

for(i in 1:length(q)){
  q[[i]]$legend<-NULL
}

# check that no ancestral states or tip states are being inferred and check rootstate settings

sapply(X=q, `[[`, "root.p")
sapply(X=q, `[[`, "node.states")
sapply(X=q, `[[`, "get.tip.states")


### force roots to be fixed as semi-terrestrial

for(i in 1:length(q)){
  q[[i]]$root.p<-c(0,1,0)
}


##############################################################################################
# RUN corHMM

## write function to run corHMM and save output
run.corHMM<-function(i) {
  tmp<-do.call(corHMM, args=q[[i]])
  if(!dir.exists("corHMM_output")) dir.create("corHMM_output")
  saveRDS(tmp,file = paste0("corHMM_output/",i,"_",names(q)[i], ".rds"))
  return(tmp)
}

#test
#run.corHMM(i = 1)

# loop through all models in parallel

corHMM_fit<-list()
corHMM_fit<-mclapply(FUN=run.corHMM,
           X=1:length(q2),
           mc.cores = 4)  # warning!!! each model already uses 25 cores!

names(corHMM_fit)<-names(q)


########
save.image("corHMM_gymno.RData")
saveRDS(corHMM_fit, "corHMM_fit_gymno.rds")
########



##############################################################################################
# summarize models


fit_sum<-data.frame(models=names(corHMM_fit),
                    loglik=sapply(corHMM_fit[sapply(corHMM_fit, length)>1], `[[`, "loglik")[names(corHMM_fit)],
                    AIC=sapply(corHMM_fit[sapply(corHMM_fit, length)>1], `[[`, "AIC")[names(corHMM_fit)],
                    AICc=sapply(corHMM_fit[sapply(corHMM_fit, length)>1], `[[`, "AICc")[names(corHMM_fit)])

# add Akaike weights
AcW<-Weights(fit_sum$AICc[!is.na(fit_sum$AICc)])
names(AcW)<-fit_sum$models[!is.na(fit_sum$AICc)]
fit_sum$AcW<-AcW[as.character(fit_sum$models)]

Aw<-Weights(fit_sum$AIC[!is.na(fit_sum$AIC)])
names(Aw)<-fit_sum$models[!is.na(fit_sum$AIC)]
fit_sum$Aw<-Aw[as.character(fit_sum$models)]

fit_sum<-fit_sum[order(fit_sum$AcW, decreasing = T, na.last = T),]

head(fit_sum)


# export summary
write.csv(fit_sum, "gymno_fit_summary.csv", row.names = T)


########
save.image("corHMM_gymno.RData")
########
