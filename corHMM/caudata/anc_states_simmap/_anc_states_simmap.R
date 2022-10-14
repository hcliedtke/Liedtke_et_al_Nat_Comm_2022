##############################################################################################
##############################################################################################
##############################################################################################
### corHMM: ancestral states [stochastic character mapping]

# This script was used to reconstruct ancestral states using stochasitc character mapping, for the Caudata
# dataset of Liedtke et al. 2022.
##############################################################################################
##############################################################################################
##############################################################################################

### set working directory and load libraries

setwd("~/Documents/amphibian_diversity_project/2021/_SI_scripts/corHMM/caudata/anc_states_simmap/")


library(ape)
library(geiger)
library(phytools)
library(corHMM)
library(tidyverse)
library(qgraph)
library(phytools)


### load custom functions

source("../../aux_scripts/find_transitions.R")


##############################################################################################
# load corHMM results from previous script [_run_corHMM_caudata.R]
corHMM_fit<-readRDS("../corHMM_fit_caudata.rds")

fit_sum<-read_csv("../caudata_fit_summary.csv")
head(fit_sum)

#############################################################################################
#### defined the best fitting model (based on AIC) and extract model results

best_fit<-fit_sum[min(fit_sum$AIC)==fit_sum$AIC,]$models
best_fit

best_fit<-corHMM_fit[[best_fit]]
#############################################################################################

best_fit_sim<-makeSimmap(tree = best_fit$phy,
                         data = best_fit$data,
                         model = best_fit$solution,
                         rate.cat = best_fit$rate.cat,
                         root.p= best_fit$root.p,
                         nSim=999, # this may take a while, see below for loading pre-baked run
                         nCores=8)


class(best_fit_sim)<-c("multiSimmap","multiPhylo")

## describe simmap
best_fit_simmap<-describe.simmap(best_fit_sim)

# save
saveRDS(best_fit_sim, "anc_states_sim.rds")
saveRDS(best_fit_simmap, "anc_states_simmap.rds")

# load saved run 
best_fit_sim<-readRDS("anc_states_sim.rds")
best_fit_simmap<-readRDS("anc_states_simmap.rds")


##############################################################################################
##### plot
##############################################################################################


### define colour scheme

rep.cols1<-c("A"="#24b9e9",
             "S"="#009E73",
             "T"="#BF9C00",
             "D"="#55FF7E",
             "V"="#f6776f",
             "P"="gold",
             "none"="#ebebeb")
rep.cols1

### make painted map for final tree

# define states
states<-c("A","D","P","S","V")
names(states)<-as.character(1:5)

## extract node states with highest probability
node_states_max<-states[as.character(apply(X = best_fit_simmap$ace, M=1, FUN = which.max))]

# define tip states
tip_states<-best_fit$data$rep_mode1
tip_states[grepl(tip_states, pattern="\\&")]<-NA # set all unknowns to NA

phy_painted<-find_transitions(phy=best_fit$phy,
                              tip_states = tip_states,
                              node_states = node_states_max,
                              simmap = T,
                              stem_prop = 1)


pdf("caudata_best_stochmap.pdf", paper="a4", width=8, height=11)
par(mfrow=c(1,1))
par(mar=c(1,1,1,1))
plotSimmap(phy_painted, fsize=0.1, lw=1, colors = rep.cols1)
tiplabels(pch=15, col=rep.cols1[tip_states], cex=0.2, adj = 1.2)
nodelabels(pie=best_fit_simmap$ace, piecol= rep.cols1[states], cex = 0.25)
dev.off()


##############################################################################################
# Plot the number of transitions

### count the number of transitions mapped on the tree:

tr_tbl<-function(x) {
  tr<-x %>%
    as.data.frame() %>%
    select(-N) %>%
    summarise_all(mean) %>%
    pivot_longer(everything(), names_to="trans", values_to="n") %>%
    mutate(from=str_remove(trans, ",\\w+"),
           to=str_remove(trans, "\\w+,")) %>%
    select(from, to, n) %>%
    filter(n>0) %>%
    mutate(from=states[from],
           to = states[to])
  
  return(tr)
}


(best_tr<-tr_tbl(x=best_fit_simmap$count) )

# make plotting function
#### see full list of graphic arguments here: https://cran.r-project.org/web/packages/qgraph/qgraph.pdf
tr_plots<-function(tr_counts, states, title) {
  
  
  qgraph(title=title,
         tr_counts,
         #nodes:
         labels=states[c(1,4,2,3,5)],
         shape="rectangle",
         vsize=10,
         node.height=0.5,
         label.cex=1.25,
         label.color=rep.cols1[states],
         # edges:
         edge.labels=T,
         edge.label.cex=2,
         #minimum=0.01, # threshold for edges to include
         edge.width=1.5,
         posCol="black",
         fade=T,
         colFactor=0.2,
         # themes and layout:
         #theme="TeamFortress",
         layout="circular",
         directed=TRUE)
}

pdf("caudata_transitions_stochmap.pdf",width=11, height=8)
par(mar=c(1,2,1,1))
par(mfrow=c(1,2))
tr_plots(tr_counts = best_tr, states=states, title="Mean Stochastic transitions")
dev.off()


