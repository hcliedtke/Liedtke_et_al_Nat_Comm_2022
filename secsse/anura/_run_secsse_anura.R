##############################################################################################
##############################################################################################
##############################################################################################
### SecSSE

# This script was used to model trait evolution using the secsse package, for the Anura
# dataset of Liedtke et al. 2022.
##############################################################################################
##############################################################################################
##############################################################################################

##############################################################################################
# set working directory and load required libraries

setwd("~/Documents/amphibian_diversity_project/2021/_SI_scripts/secsse/anura/")

library(pacman)
p_load(ape,secsse,DDD,doMC,tidyverse,parallel,qgraph)

##############################################################################################
## LOAD DATA AND MODELS

phy_full<-read.tree("../amphibia.tre")
dat_full<-read_csv("../rep_modes.csv")

#############################
## Prune to only anura

# adjust data
dat_anura<-dat_full %>%
  filter(order=="Anura") %>%
  mutate(secsse_state=case_when(
    rep_mode=="A"~1,
    rep_mode=="S"~2,
    rep_mode=="T"~3,
    rep_mode=="D"~4,
    rep_mode=="L"~5))


table(dat_anura$secsse_state,useNA = "always") ### compare
table(dat_anura$rep_mode, useNA="always")

# pair tree and data
phy<-drop.tip(phy_full, as.character(dat_full$species_asw[dat_full$order!="Anura"]))
phy

dat<-dat_anura %>%
  filter(species_asw %in% phy$tip.label) %>%
  arrange(match(species_asw, phy$tip.label))
all(dat$species_asw==phy$tip.label)

# make secsse states
states<-pull(dat, secsse_state)
states
table(states,useNA = "always")

##############################################################################################
## number of states and sampling fraction

# number of observed states
n_states<-length(unique(na.omit(states)))


## sampling fraction
f=table(dat$secsse_state)/table(dat_anura$secsse_state)
f

#############################################################
# LOAD MODELS

## loads auxiliary functions
source("../aux_scripts/secsse_functions.R")
# loads the models/hypotheses to be tested
source("../aux_scripts/secsse_base_models.R")


length(models) # should be 4 models plus 4 fixed mu0 models
names(models)


#############################################################
### MODIFY MODELS - fix all transitions to and from viviparity and from terrestrial to be symmetrical (with up to 3 hidden states)

for(i in 1:length(models)){
  
  models[[i]]$idparslist$Q[3,1]<-models[[i]]$idparslist$Q[1,3]
  models[[i]]$idparslist$Q[3,2]<-models[[i]]$idparslist$Q[2,3]
  models[[i]]$idparslist$Q[3,4]<-models[[i]]$idparslist$Q[4,3]
  models[[i]]$idparslist$Q[3,5]<-models[[i]]$idparslist$Q[5,3]
  
  models[[i]]$idparslist$Q[5,1]<-models[[i]]$idparslist$Q[1,5]
  models[[i]]$idparslist$Q[5,2]<-models[[i]]$idparslist$Q[2,5]
  models[[i]]$idparslist$Q[5,3]<-models[[i]]$idparslist$Q[3,5]
  models[[i]]$idparslist$Q[5,4]<-models[[i]]$idparslist$Q[4,5]
  
  
  if(nrow(models[[i]]$idparslist$Q)>n_states) {
    
    models[[i]]$idparslist$Q[8,6]<-models[[i]]$idparslist$Q[6,8]
    models[[i]]$idparslist$Q[8,7]<-models[[i]]$idparslist$Q[7,8]
    models[[i]]$idparslist$Q[8,9]<-models[[i]]$idparslist$Q[9,8]
    models[[i]]$idparslist$Q[8,10]<-models[[i]]$idparslist$Q[10,8]
    
    models[[i]]$idparslist$Q[10,6]<-models[[i]]$idparslist$Q[6,10]
    models[[i]]$idparslist$Q[10,7]<-models[[i]]$idparslist$Q[7,10]
    models[[i]]$idparslist$Q[10,8]<-models[[i]]$idparslist$Q[8,10]
    models[[i]]$idparslist$Q[10,9]<-models[[i]]$idparslist$Q[9,10]
    
  }
  if(nrow(models[[i]]$idparslist$Q)>n_states*2) {
    
    models[[i]]$idparslist$Q[13,11]<-models[[i]]$idparslist$Q[11,13]
    models[[i]]$idparslist$Q[13,12]<-models[[i]]$idparslist$Q[12,13]
    models[[i]]$idparslist$Q[13,14]<-models[[i]]$idparslist$Q[14,13]
    models[[i]]$idparslist$Q[13,15]<-models[[i]]$idparslist$Q[15,13]
    
    models[[i]]$idparslist$Q[15,11]<-models[[i]]$idparslist$Q[11,15]
    models[[i]]$idparslist$Q[15,12]<-models[[i]]$idparslist$Q[12,15]
    models[[i]]$idparslist$Q[15,13]<-models[[i]]$idparslist$Q[13,15]
    models[[i]]$idparslist$Q[15,14]<-models[[i]]$idparslist$Q[14,15]
    
  }
  
  models[[i]]$idparslist$Q<-reindex_q(models[[i]]$idparslist$Q)
  
}

#############################################################
# MODIFY MODELS - set specific transition rates to 0 to leave only the best resulting Q matrix from the corHMM analysis

### prepare a mask that can be overlaid the Q matrix to set specific transitions to 0
mask<-matrix(1, nrow=5, ncol=5)
diag(mask)<-NA
mask[1,3]<-0
mask[3,1]<-0
mask[1,4]<-0
mask[4,1]<-0
mask[2,5]<-0
mask[5,2]<-0
mask[3,4]<-0
mask[4,3]<-0
mask[3,5]<-0
mask[5,3]<-0
mask[4,5]<-0
mask[5,4]<-0

# visualize
#qgraph(mask, directed=TRUE,layout="circular")

# mask and reindex Q matrices

for(i in 1:length(models)){
  models[[i]]$idparslist$Q<-mask_q(q=models[[i]]$idparslist$Q, mask=mask, n_states = n_states)
}

######################################
### update number of parameters to optimize

for(i in 1:length(models)){
  
  idparsopt = sort(unique(na.omit(unname(unlist(models[[i]]$idparslist)))))
  if(any(idparsopt %in% 0)){
    idparsopt = idparsopt[idparsopt!=0]
    idparsfix = c(0)
    parsfix = c(0)
  } else {
    idparsfix = c()
    parsfix = c()
  }
  
  ## append
  models[[i]]<-append(models[[i]],
                      list(
                        idparsopt=idparsopt,
                        idparsfix=idparsfix,
                        parsfix=parsfix
                      )
  )
}

###############################################################################################
## set different starting parameters for models. 
## replicate models three times, implementing different starting values for lambda, mu and Q.

## try1: optimized
## try2: optimized/2
## try3: optimizes*2

startingpoint <- bd_ML(brts = ape::branching.times(phy))
intGuessLambda <- startingpoint$lambda0
intGuessMu <- startingpoint$mu0

intGuessLambdas<-c(intGuessLambda, intGuessLambda/2, intGuessLambda*2)
intGuessMus <- c(intGuessMu, intGuessMu/2, intGuessMu*2)
initTrans <- intGuessLambdas/n_states

nrep=length(intGuessLambdas)

# replicate models and rename
replicate_models<-rep(models,each=nrep)
names(replicate_models)<-paste0(rep(names(models), each=nrep), rep(paste0("_try", 1:nrep),nrep))
names(replicate_models)

# replace models object
models<-replicate_models


# add starting values 
iterator=rep(1:nrep, length(models)/nrep)
for(i in 1:length(models)){

  ## prep:
  initparsopt = c(
    rep(intGuessLambdas[iterator[i]],length(unique(models[[i]]$idparslist$lambdas))),
    ifelse(unique(models[[i]]$idparslist$mus)!=0,
      rep(intGuessMus[iterator[i]],length(unique(models[[i]]$idparslist$mus))),
      0
    ),
    rep(initTrans[iterator[i]],length(which(unique(c(models[[i]]$idparslist$Q))>0)))
  )
  initparsopt<-initparsopt[initparsopt!=0]

  ## append
  models[[i]]$initparsopt<-initparsopt

}

## check idparsopt is the same length as initparsopt

all(sapply(models, FUN = function(x) length(x$initparsopt)==length(x$idparsopt)))


# one final check:
lapply(models, `[[`, "idparslist")


###############################################################################################
### make a function to run secsse in a parallelized loop

## add remaining arguments for secsse_ml() function
models<-lapply(models, FUN=append,
               list(phy = phy,
                    traits = states,
                    cond="proper_cond",
                    root_state_weight = c(1,0,0,0,0),
                    sampling_fraction=f,
                    num_cycles = 3, ### WARNING!!!! set this to something higher!!! - takes a very long time
                    run_parallel=T)
               )

## export lists of starting parameters and model specifications
saveRDS(models,"secsse_models_anura.rds")

# make run function
run_secsse<-function(i) {

  if(!dir.exists("secsse_out")) dir.create("secsse_out")

  try(
    saveRDS(do.call(secsse_ml, args=models[[i]]),
            file = paste0("secsse_out/",names(models)[i],"_out.rds")))
}

###############################################################################################
# RUN MODELS
###############################################################################################

# set number of cores
n_cores=ifelse(length(models)>30, 30, length(models))
n_cores

# test run single model
#do.call(secsse_ml, args = models[[1]])

# run in parallel
mclapply(FUN=run_secsse,
         X=1:length(models),
         mc.cores = n_cores)


###############################################################################################
## Evaluate models 
###############################################################################################

## load secsse output files
files<-list.files(path = "secsse_out", pattern=".*_out.rds")
files
secsse_out<-lapply(paste0("secsse_out/",files), readRDS)
names(secsse_out)<-str_remove_all(files, pattern="_out.rds")

## check all models match
names(secsse_out) %in% names(models)
names(models) %in% names(secsse_out)

## make sure that all have converged
all(sapply(secsse_out, `[[`, "conv") == 0) # should all be 0 (i.e. TRUE)

## tabulate likelihoods
secsse_fit_itr<-data.frame(model_type=str_remove(names(secsse_out), pattern="_try\\d+$"),
                           model_itr = as.numeric(replace_na(str_extract(names(secsse_out), pattern="\\d+$"),"0")),
                           ml = sapply(secsse_out, `[[`, "ML"))

secsse_fit_itr

## keep only models with highest likelihood per try
secsse_fit <- secsse_fit_itr %>%
  rownames_to_column("model_name") %>%
  group_by(model_type) %>%
  filter(ml == max(ml)) %>%
  distinct(model_type, .keep_all=TRUE)

secsse_fit

### filter out mu0 models

secsse_fit <- secsse_fit %>%
  filter(!str_detect(model_name, "mu0"))

secsse_fit
## calculate AIC scores and Aw

### functions for manual calculations:
calc_AIC<-function(loglik, np) { -2 * loglik + 2 * np}
calc_AW <- function(dAIC) { exp(-0.5*dAIC)/sum(exp(-0.5*dAIC)) }

secsse_fit$AIC<-NA
for(i in 1:nrow(secsse_fit)) {
  secsse_fit$AIC[i]<-calc_AIC(loglik = secsse_fit$ml[[i]],
                              np = length(models[[secsse_fit$model_name[[i]]]]$idparsopt))
}

secsse_fit$dAIC<-secsse_fit$AIC-min(secsse_fit$AIC)
secsse_fit$Aw<-calc_AW(dAIC = secsse_fit$dAIC)

secsse_fit<-secsse_fit %>%
  arrange(desc(Aw))

secsse_fit

write_csv(secsse_fit,"secsse_fit_summary.csv")

###############################################################################################
## plot model rate estimates

#### tabulate lambda, mu and div rate for best model
secsse_fit$model_name
best<-secsse_fit$model_name[1]
best
best_model<-secsse_out[[best]]


best_params<-best_model$MLpars[[1]] %>%
  enframe(name = "state", value = "lambda") %>%
  bind_cols(
    best_model$MLpars[[2]] %>%
      enframe(name = "state_mu", value = "mu") 
  ) %>% #### BE CAREFUL HERE THAT THE COLUMNS BIND IN THE CORRECT ORDER
  select(-state_mu) %>%
  mutate(div=lambda-mu) %>%
  pivot_longer(-state) %>%
  mutate(state=factor(state, levels=unique(state)),
         hidden=str_remove(state, pattern="\\d+"), # if hidden states exist
         observed=str_extract(state, pattern="\\d+")) %>%
  mutate(observed = c("1"="A",
                      "2"="S",
                      "3"="D",
                      "4"="L",
                      "5"="P")[as.character(observed)])
best_params

#export
write_csv(best_params, "best_model_div_rates.csv")
#### plot observed and hidden side by side (IF APPLICABALE)

## colour scheme
rep.cols1<-c("A"="#24b9e9",
             "S"="#009E73",
             "T"="#BF9C00",
             "D"="#55FF7E",
             "L"="#f6776f",
             "P"="gold",
             "none"="#ebebeb")
rep.cols1

best_params %>%
  ggplot(aes(x=observed, y=value, fill=observed, group=hidden)) +
  geom_bar(stat = "identity", position="dodge", color="black") +
  scale_fill_manual(values=rep.cols1) +
  facet_wrap(~name) +
  ylab("Rate") +
  xlab("") + 
  theme(legend.position = "none")

# plot as averages across hidden and observed states if applicable
best_params %>%
  group_by(name, observed) %>%
  summarise(value=mean(value)) %>%
  ggplot(aes(x=observed, y=value, fill=observed)) +
  geom_bar(stat = "identity", position="dodge") +
  scale_fill_manual(values=rep.cols1) +
  facet_wrap(~name, scales = "free_y") +
  coord_flip() +
  theme(legend.position = "none")
