##############################################################################################
##############################################################################################
##############################################################################################
### corHMM: ancestral states [join estimates]

# This script was used to reconstruct ancestral states using joint estimates, for the Anura
# dataset of Liedtke et al. 2022.
##############################################################################################
##############################################################################################
##############################################################################################

### set working directory and load libraries
setwd("~/Documents/amphibian_diversity_project/2021/_SI_scripts/corHMM/anura/anc_states_joint/")

library(ape)
library(geiger)
library(phytools)
library(corHMM)
library(tidyverse)
library(qgraph)

### load custom functions
source("../../aux_scripts/find_transitions.R")

##############################################################################################
# load corHMM results from previous script [_run_corHMM_caudata.R]
corHMM_fit<-readRDS("../corHMM_fit_anura.rds")

fit_sum<-read_csv("../anura_fit_summary.csv")
head(fit_sum)

#############################################################################################
#### defined the best fitting model (based on AIC) and extract model results

best_fit<-fit_sum[min(fit_sum$AIC)==fit_sum$AIC,]$models
best_fit

best_fit<-corHMM_fit[[best_fit]]
#############################################################################################
#### function for preparing parameters in the right order for ancRECON

#- index: index matrix
#- q: parameter matrix

p<-function(index, q) {
  sapply(1:max(index, na.rm = TRUE),
         function(x) na.omit(c(q))[na.omit(c(index) == x)][1])
}

##############################################################################################
# reconstruct ancestral states

## extract parameters needed:
best_fit_params <- p(index=best_fit$index.mat, q=best_fit$solution)
best_fit_params

best_fit_anc<-ancRECON(phy = best_fit$phy,
                      data = best_fit$data,
                      p = best_fit_params,
                      method = "joint",
                      get.tip.states = TRUE,
                      rate.cat = best_fit$rate.cat,
                      rate.mat = best_fit$index.mat,
                      root.p = best_fit$root.p)

## add states to corHMM object
best_fit$states<-best_fit_anc$lik.anc.states


###### save RData
save.image("anura_anc_joint.RData")
#load("anura_anc_joint.RData")


##############################################################################################
# painted phylogeny with ancestral states

### define tip states coding
tip_states<-best_fit$data$rep_mode1
tip_states[grepl(tip_states, pattern="\\&")]<-NA # set all unknowns to NA

## make look-up for node state coding
states.codes<-rep(c("A","D","S","T","V"),2) # double because of the best model has hidden states
names(states.codes)<-as.character(1:10)
states.codes

# translate numeric code to our lettering
node_states_best<-states.codes[as.character(best_fit$states)]


# paint tree
phy_painted_best<-find_transitions(phy=best_fit$phy,
                                   tip_states = tip_states,
                                   node_states = node_states_best,
                                   simmap = T,
                                   stem_prop = 1)

##############################################################################################
# Plot tree


### define colour scheme

rep.cols1<-c("A"="#24b9e9",
             "S"="#009E73",
             "T"="#BF9C00",
             "D"="#55FF7E",
             "V"="#f6776f",
             "P"="gold",
             "none"="#ebebeb")
rep.cols1



pdf("anura_best_joint.pdf",width=25, height=75)
par(mfrow=c(1,1))
par(mar=c(1,1,1,1))
plotSimmap(phy_painted_best, fsize=0.1, lw=1, colors = rep.cols1, offset = 0.2)
nodelabels(pch=19, col=rep.cols1[states.codes[best_fit$states]], cex=0.5)
tiplabels(pch=15, col=rep.cols1[tip_states], cex=0.2, adj = 0.7)
legend("bottomleft", bty="n", pch=19, col=rep.cols1[unique(tip_states)],
       legend=unique(tip_states),
       title="Life-history modes")
dev.off()






##############################################################################################
# Plot the number of transitions and the transition rates



### count the number of transitions mapped on the tree:
(best_tr<-countSimmap(phy_painted_best)$Tr)


# make plotting function
#### see full list of graphic arguments here: https://cran.r-project.org/web/packages/qgraph/qgraph.pdf
tr_plots<-function(tr_counts, states=c("A","S","T","D","V"), title) {
  qgraph(title=title,
         tr_counts[states,states],
         #nodes:
         labels=states,
         shape="rectangle",
         vsize=10,
         node.height=0.5,
         label.cex=1.25,
         label.color=rep.cols1[states],
         # edges:
         edge.labels=T,
         edge.label.cex=2,
         #minimum=0.001, # threshold for edges to include
         edge.width=1.5,
         posCol="black",
         fade=T,
         colFactor=0.2,
         # themes and layout:
         #theme="TeamFortress",
         layout="circular",
         directed=TRUE)
}

# plotting function for transition rates
rate_plots<-function(tr_rates,states=c("A","S","T","D","V"), state_codes, title) {
  
  ordered_rates<-tr_rates[match(states, state_codes), match(states, state_codes)]
  
  qgraph(title=title,
         ordered_rates,
         #nodes:
         labels=states,
         shape="rectangle",
         vsize=10,
         node.height=0.5,
         label.cex=1.25,
         label.color=rep.cols1[states],
         # edges:
         edge.labels=ordered_rates,
         edge.label.cex=2,
         #minimum=0.001, # threshold for edges to include
         edge.width=1.5,
         posCol="black",
         fade=T,
         colFactor=0.2,
         # themes and layout:
         #theme="TeamFortress",
         layout="circular",
         directed=TRUE)
}



pdf("anura_transitions_joint.pdf",width=11, height=8)
par(mar=c(1,2,1,1))
par(mfrow=c(1,3))
tr_plots(best_tr, title="Number of transitions")
rate_plots(round(best_fit$solution[1:5,1:5],4), state_codes=states.codes, title="Transition rates (observed)")
rate_plots(round(best_fit$solution[6:10,6:10],4), state_codes=states.codes, title="Transition rates (hidden")
dev.off()

save.image("anura_anc_joint.RData")

