################################################################################
################################################################################
#### RUN MUHISSE
################################################################################
################################################################################

# this script runs muhisse on the anura dataset


# following vignette from here: https://cran.r-project.org/web/packages/hisse/vignettes/hisse-vignette.pdf

library(parallel)
library(ape)
library(hisse)
library(diversitree)
library(geiger)
library(tidyverse)
library(phytools)

setwd("~/Documents/amphibian_diversity_project/2021/_SI_scripts/muhisse/anura/")
pdf("muhisse_anura.pdf", paper="a4")

################################################################################
## LOAD DATA AND MODELS
################################################################################

phy<-read.tree("./anura.tre")
dat.all<-read.csv("./anura.csv", row.names = 1) ## all data
dat<-dat.all[phy$tip.label,]
name.check(phy, dat) ## should be OK

# load custom functions used by the scripts
source("../aux_scripts/muhisse_functions.R")

# load the model designs to be tested
source("../aux_scripts/muhisse_models_2022.R")

################################################################################
## MAKE MUHISSE DATASET
################################################################################

states<-data.frame(sp=rownames(dat), larva=dat$SSE_larva, adult=dat$SSE_adult)
head(states)

table(states$larva, useNA = "ifany")
table(states$adult, useNA = "ifany")
any(is.na(states)) ## warning!! should not contain NA's
table(states[,c("larva","adult")])


## do a final plot to check things look right
plot(phy, show.tip.label=F)
tiplabels(pch=15,col=ifelse(states$larva==1, "mediumvioletred","gold"), cex=0.5)
tiplabels(pch=15,col=ifelse(states$adult==1, "mediumvioletred","gold"), cex=0.5, adj = 3)

# calculate sampling fraction
f<-c("00"=length(which(dat$SSE_larva==0 & dat$SSE_adult==0))/length(which(dat.all$SSE_larva==0 & dat.all$SSE_adult==0)),
     "01"=length(which(dat$SSE_larva==0 & dat$SSE_adult==1))/length(which(dat.all$SSE_larva==0 & dat.all$SSE_adult==1)),
     "10"=length(which(dat$SSE_larva==1 & dat$SSE_adult==0))/length(which(dat.all$SSE_larva==1 & dat.all$SSE_adult==0)),
     "11"=length(which(dat$SSE_larva==1 & dat$SSE_adult==1))/length(which(dat.all$SSE_larva==1 & dat.all$SSE_adult==1))
)
f
f[is.na(f)]<-1 # update any NAs for states that don't exist (e.g. aquatic eggs and terrestrial larvae)
f

################################################################################
## append run-specific arguments to models
################################################################################

models<-lapply(models, FUN=append,
               list(f=f,
                    phy=phy,
                    data=states,
                    root.type="madfitz",
                    root.p=NULL))


### set all 00 and all 10 states in the rate vectors and transition matrix to 0

for(i in 1:length(models)){
  models[[i]]$trans.rate<-trans.correct(models[[i]]$trans.rate, index=1)
  models[[i]]$turnover<-drop.rate(models[[i]]$turnover, index=1)
  models[[i]]$eps<-drop.rate(models[[i]]$eps, index=1)
}

for(i in 1:length(models)){
  models[[i]]$trans.rate<-trans.correct(models[[i]]$trans.rate, index=3)
  models[[i]]$turnover<-drop.rate(models[[i]]$turnover, index=3)
  models[[i]]$eps<-drop.rate(models[[i]]$eps, index=3)
}


### set correct root state

for(i in 1:length(models)){
  models[[i]]$root.p<-set_root(n.hidden = (length(models[[i]]$turnover)/4), root = 4 )
}


## add names to models for StartingValueTry function

for(i in 1:length(models)){
  models[[i]]$name<-as.character(models.meta$model.name)[i]
}

################################################################################
# PARALLELIZE AND RUN
################################################################################

## set number of cores to use (max 20) and a seed for mclapply
#ncores=ifelse(length(models)>20, 20, length(models))
ncores=3 # remember that each StartingValueTry will take 25 cores
set.seed(12345)


## set up function to run either MuHiSSE with one set of starting values (for models with no hidden state; MuSSE) or trying 20 different starting values (for models with starting values).
run.muhisse<-function(i) do.call(StartingValueTry, args=models[[i]])

## run...

mclapply(FUN=run.muhisse,
         X=1:length(models),
         mc.cores = ncores)

save.image("muhisse_anura.RData")


################################################################################
# PROCESS RUNS
################################################################################

file.names<-list.files("./starting_tries/")
file.path<-paste0("./starting_tries/", file.names)

all.fit<-list()
for(i in 1:length(file.names)){
  load(file.path[i])
  all.fit[[i]]<-tmp
}

starting.tries<-data.frame(file.no=1:length(all.fit),
                           name=gsub(file.names, pattern="_iter.*", replacement = ""),
                           aic=unlist(lapply(all.fit, `[[`,"AIC")),
                           loglik=unlist(lapply(all.fit, `[[`,"loglik")))


# pull out best scoring models
max.loglik<-starting.tries %>%
  group_by(name) %>%
  slice(which.max(loglik)) %>%
  arrange(match(name,as.character(models.meta$model.name)))

max.loglik

muhisse.fit<-all.fit[max.loglik$file.no]

save.image("muhisse_anura.RData")

#############################################################################################
# INSPECT MODEL PERFORMANCE

model.performance<-data.frame(model=1:length(models),
                              model_type=models.meta$model.type,
                              model_name=models.meta$model.name,
                              lnL=sapply(muhisse.fit, '[[', 1),
                              AIC=sapply(muhisse.fit, '[[', 2),
                              AICc=sapply(muhisse.fit, '[[', 3))

model.performance$dAIC<-model.performance$AIC-min(model.performance$AIC)
model.performance$AW<-exp(-0.5*model.performance$dAIC)/sum(exp(-0.5*model.performance$dAIC))

(model.performance<-model.performance[order(model.performance$dAIC),])

write.table(model.performance, "muhisse_fit_anura.txt", quote=F, sep="\t", row.names = F)

####
save.image("muhisse_anura.RData")
####

################################################################################
#  PERFORM CHARACTER RECONSTRUCTION FOR ALL MODELS
################################################################################

#load("muhisse_anura.RData")

## extract only the "best" models for each model category/type and reconstruct ancestral states and div rates

best.models<-model.performance$model[!duplicated(model.performance$model_type)]
top.models<-model.performance$model[model.performance$AW>0.01]
rec.ids<-sort(unique(c(top.models, best.models)))
rec.ids

margin.rec<-list()
for(i in 1:length(rec.ids)){
  id<-rec.ids[i]
  n.hidden<-nrow(muhisse.fit[[id]]$trans.matrix)/4
  margin.rec[[i]]<-MarginReconMuHiSSE(phy=muhisse.fit[[id]]$phy,
                                       data=muhisse.fit[[id]]$data,
                                       f=muhisse.fit[[id]]$f,
                                       hidden.states = n.hidden,
                                       pars=muhisse.fit[[id]]$solution,
                                       aic=muhisse.fit[[id]]$AIC, ### NOTE: newer versions use capitalized AIC
                                       root.type="user",
                                       root.p=models[[id]]$root.p,
                                       n.cores = 10)
  print(paste0("...done with reconstruction ", i, " of " ,length(rec.ids)))
}

names(margin.rec)<-as.character(model.performance$model_name)[match(rec.ids, model.performance$model)]

####
save.image("muhisse_anura.RData")
####

for(i in 1:length(margin.rec)){
  plot.muhisse.states(margin.rec[[i]], rate.param="net.div", show.tip.label=FALSE)
  legend("topright", bty="n", legend = names(margin.rec)[[i]])
}

####
save.image("muhisse_anura.RData")
####

################################################################################
####  PLOTS WITH utilhisse FUNCTIONS
################################################################################


library(gghisse)
library(gridExtra)
library(tidyverse)

## extract plotting data for all the reconstructed models

plot.dat<-list()
for(i in 1:length(margin.rec)){
  plot.dat[[i]]<-m_process_recon(margin.rec[[i]])
}

names(plot.dat)<-names(margin.rec)


## get model of interest [best model]

mod<-muhisse.fit[[model.performance$model[1]]]


## get model parameters
muhisse.param<-m_collect_rates(model_fit=mod,
                               hidden_traits=mod$hidden.states,
                               character_states=c("nLnA","LnA","nLA","LA"))

muhisse.param


# model average top x models

models.to.average<-as.character(model.performance$model_name) [match(model.performance$model,x=top.models)]
plot.dat$top.avg<-m_process_recon(margin.rec[models.to.average])


# diversification rates for binary states under best models for each category


for(i in 1:length(plot.dat)){
  params=c("net.div","speciation","extinction")
  s<-list()
  r<-list()

  for(j in params){
    s[[j]]<-m_scatterplot(plot.dat[[i]],
                          parameter = j,
                          states_names = c("impossibility","Direct", "Paedomorphic","Biphasic"),
                          colors = c("#2D2E83","#951B81"))


    r[[j]]<-m_ridgelines(plot.dat[[i]],
                         parameter = j,
                         states_names = c("impossibility","Direct", "Paedomorphic","Biphasic"),
                         line_colors = c("#2D2E83","#951B81"),
                         fill_colors = c(NA,NA))

  }

  grid.arrange(nrow = 2, ncol=3, grobs=c(s,r), top=names(plot.dat)[i])

  map.traits<-m_trait_recon(plot.dat[[i]],
                            show_tip_labels = FALSE,
                            states_of_first_character = c("no larva", "larva"),
                            states_of_second_character = c("no adult","adult"),
                            colors=c("gray","#2D2E83", "gold","#951B81"),
                            cutoff = 0.5,
                            tree_layout = "circular",
                            tree_direction = "right",
                            time_axis_ticks = 4,
                            open_angle = 10) +
    theme(legend.position="bottom")

  map.rates<-m_rate_recon(plot.dat[[i]],
                          show_tip_labels = FALSE,
                          parameter = "net.div",
                          discrete = FALSE,
                          breaks = seq(0, 1, 0.2),
                          colors = c("mediumvioletred","gold"),
                          plot_as_waiting_time = FALSE,
                          tree_layout = "circular",
                          tree_direction = "right",
                          time_axis_ticks = 4,
                          open_angle = 10) +
    theme(legend.position="bottom")

  grid.arrange(nrow = 1, ncol=2, map.traits,map.rates, top=names(plot.dat)[i])


}




####
save.image("muhisse_anura.RData")
####


################################################################################
### best muhisse model stats
################################################################################

names(plot.dat)
plot.dat$LhID_4$tip_rates %>%
  select(-id) %>%
  mutate(state=factor(max.col(select(.,contains("state"))), labels = c("direct","biphasic"))) %>%
  group_by(state) %>%
  summarise_all(c(mean,sd)) %>%
  print(width=400)

### best model stats

##  model stats

best.mod.stats<-plot.dat$LhID_4$tip_rates %>%
  select(-id) %>%
  mutate(state=factor(max.col(select(.,contains("state"))), labels = c("direct","biphasic"))) %>%
  group_by(state) %>%
  summarise_all(c(mean,sd)) %>%
  print(width=400)

write.csv(best.mod.stats, "best_mod_tiprates_anura.csv")

################################################################################
### Plot tip rates on tree
################################################################################

par(mfrow=c(1,4))
par(mar=c(4,1,4,1))
rep.cols=c("00"="grey", "10"="gold", "01"="#2D2E83", "11"="#951B81")

plot(phy, show.tip.label = F, main="LhD tip rates")
plot(x=plot.dat$LhD_4$tip_rates$net.div, y=1:Ntip(phy),
     pch=16,
     col=rep.cols[apply(MARGIN = 1, X=states[,2:3], paste0, collapse="")],
     bty="n",
     xlab = "r",ylab=NA, yaxt="n",
     main="Net Diversification")
plot(x=plot.dat$LhD_4$tip_rates$turnover, y=1:Ntip(phy),
     pch=16,
     col=rep.cols[apply(MARGIN = 1, X=states[,2:3], paste0, collapse="")],
     bty="n",
     xlab = "t",ylab=NA, yaxt="n",
     main="turnover")
plot(x=plot.dat$LhD_4$tip_rates$extinct.frac, y=1:Ntip(phy),
     pch=16,
     col=rep.cols[apply(MARGIN = 1, X=states[,2:3], paste0, collapse="")],
     bty="n",
     xlab = "e",ylab=NA, yaxt="n",
     main="ext.fr")

dev.off()

####
save.image("muhisse_anura.RData")
####

