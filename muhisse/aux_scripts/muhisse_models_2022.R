################################################################################
################################################################################
#### SET UP BiSSE/HISSE/MuHiSSE MODELS
################################################################################
################################################################################

# following vignette from here: https://cran.r-project.org/web/packages/hisse/vignettes/hisse-vignette.pdf

library(hisse)
library(diversitree)


################################################################################
# SET UP MODEL LIST

models<-list()

################################################################################
### CD models: constant diversification rate 
################################################################################


#turnover and extinction fraction different for each state, full Q matrix
models[[1]]<-list(hidden.states=FALSE,
                  turnover=c(1,1,1,1),
                  eps=c(1,1,1,1),
                  trans.rate= TransMatMakerMuHiSSE(hidden.traits = 0,include.diagonals = T))

################################################################################
### LhD models: diversification rates vary across life history traits, with and without hidden states
################################################################################

#LhD_1: No hidden states
models[[2]]<-list(hidden.states=FALSE,
                  turnover=c(1,2,3,4),
                  eps=c(1,2,3,4),
                  trans.rate= TransMatMakerMuHiSSE(hidden.traits = 0, include.diagonals = T))

#LhD_2: Observed + 1 hidden state (2 rate categories)
models[[3]]<-list(hidden.states=TRUE,
                  turnover=c(1:8),
                  eps=c(1:8),
                  trans.rate=TransMatMakerMuHiSSE(hidden.traits = 1,include.diagonals = T))

#LhD_3: Observed + 2 hidden state (3 rate categories)
models[[4]]<-list(hidden.states=TRUE,
                  turnover=c(1:12),
                  eps=c(1:12),
                  trans.rate=TransMatMakerMuHiSSE(hidden.traits = 2,include.diagonals = T))

#LhD_4: Observed + 3 hidden state (4 rate categories)
models[[5]]<-list(hidden.states=TRUE,
                  turnover=c(1:16),
                  eps=c(1:16),
                  trans.rate=TransMatMakerMuHiSSE(hidden.traits = 3,include.diagonals = T))



################################################################################
### LhID:  diversification rates vary across hidden states (i.e. indepdentent form life history modes) ["CID" models for the hisse nomenclature]
################################################################################

#LhID_2: Observed + 1 hidden state (2 rate categories) [like CID2]
models[[6]]<-list(hidden.states=TRUE,
                  turnover=rep(1:2, each=4),
                  eps=rep(1:2, each=4),
                  trans.rate=TransMatMakerMuHiSSE(hidden.traits = 1, make.null = T, include.diagonals = T))


#LhID_3: Observed + 2 hidden state (3 rate categories) [like CID3]
models[[7]]<-list(hidden.states=TRUE,
                  turnover=rep(1:3, each = 4),
                  eps=rep(1:3, each = 4),
                  trans.rate=TransMatMakerMuHiSSE(hidden.traits = 2, make.null = T, include.diagonals = T))


#LhID_4: Observed + 3 hidden state (4 rate categories) [like CID4]
models[[8]]<-list(hidden.states=TRUE,
                  turnover=rep(1:4, each = 4),
                  eps=rep(1:4, each = 4),
                  trans.rate=TransMatMakerMuHiSSE(hidden.traits = 3, make.null = T, include.diagonals = T))



################################################################################
### define model names and types
################################################################################

models.meta<-data.frame(model.number=1:length(models),
                        model.name=c("CD",
                                     "LhD_1","LhD_2","LhD_3","LhD_4",
                                     "LhID_2","LhID_3","LhID_4"
                                     ),
                        model.type=c("CD", rep("LhD",4), rep("LhID",3)))

names(models)<-models.meta$model.name
models
