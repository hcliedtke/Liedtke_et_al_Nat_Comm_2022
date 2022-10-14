################################################################################
################################################################################
### MuHiSSE functions
################################################################################
################################################################################

# this script defines custom functions that are used in the muhisse work flow

################################################################################
### try multiple starting values
################################################################################

StartingValueTry <-
  function(phy=phy,
           data=data,
           f=f,
           trans.rate=trans.rate,
           root.type = root.type,
           root.p = root.p,
           turnover=turnover,
           eps=eps,
           hidden.states=hidden.states,
           name=name) {
    freqs <- table(apply(data[, 2:3], 1,
                         function(x)
                           switch(
                             paste0(x, collapse = ""),
                             "00" = 1,
                             "01" = 2,
                             "10" = 3,
                             "11" = 4
                           )))
    
    ### ADDED FOLLOWING CODE TO SORT OUT PROBLEMS IF ANY STATES ARE MISSING
    if(length(freqs)<4){
      freqs<-c(freqs, rep(0, 4-length(freqs)))
      names(freqs)[names(freqs)==""]<-as.character(1:4)[!1:4 %in% as.numeric(names(freqs))]
      freqs<-freqs[order(names(freqs))]
    }
    ####
    
    samp.freq.tree <- Ntip(phy) / sum(freqs / f)
    init.pars <-
      hisse:::starting.point.generator(phy, 4, samp.freq.tree, yule = FALSE)

    turnover <- turnover
    eps <- eps
    NewStarting <- function(iteration) {
      turn.start <- exp(rnorm(4, log(init.pars[1] + init.pars[5]), 1))
      eps.start <- runif(1, 0, 1)
      trans.start <- exp(rnorm(12, log(init.pars[9])))
      starting.vals <- c(turn.start, rep(eps.start, 4), trans.start)
      print(starting.vals)
      tmp <-
        MuHiSSE(
          phy,
          data,
          f = f,
          turnover = turnover,
          eps = eps,
          trans.rate = trans.rate,
          hidden.states = hidden.states,
          starting.vals = starting.vals,
          root.type = root.type
        )
      if(!dir.exists("starting_tries")) dir.create("starting_tries")
      save(tmp,
           file = paste0("starting_tries/",name, "_iter", iteration, ".Rsave"))
    }
    mclapply(1:10, NewStarting, mc.cores = 10) ## set the number of tries and cores to use
}



################################################################################
### Dropping parameters from models
################################################################################


## e.g.: state 01 (aquatic egg but no larvae) does not exist. For this reason. all rates will be set to 0 and transitions to and from will also be set to 0. to aid with this, at the end of script, all transition matrices in the models() list will be revised/corrected:

### function to correct Q matrix
trans.correct<-function(x, index){
  # set indexed transitions to 0
  i=seq(from=index, to=nrow(x), 4) ## multiple if hidden states
  x[,i]<-0
  x[i,]<-0

  # break down into sub matrices and update inner matrices
  last.id<-0
  for(k in seq(1, nrow(x), 4)){

    # extract submatrix
    sub.mat<-x[k:(k+3),k:(k+3)]
    # extract parameter ids
    params<-sub.mat[sub.mat>0 & !is.na(sub.mat)]

    # renumber parameter ids
    params[1]<-ifelse(params[1]<=last.id, 1, last.id+1)
    for(j in 2:length(params)){
      params[j]<-ifelse(params[j-1]==params[j] | params[j-1]==params[j]+1, params[j], params[j-1]+1)
    }
    sub.mat[sub.mat>0 & !is.na(sub.mat)]<-params

    # keep tabs on the last id
    last.id<-c(max(sub.mat, na.rm = T))

    #re-insert sub matrix
    x[k:(k+3),k:(k+3)]<-sub.mat
  }

  # update all values outside the diagonal sub-matrices
  x[x>last.id & !is.na(x)]<-last.id+1


  # make sure diagonals NA and return
  diag(x)<-NA
  return(x)
}


################################################################################
### function to correct speciation and extinction  rate vectors
################################################################################


drop.rate<-function(x,index){
  i=index
  # set indexed trait to 0
  x[seq(1, length(x), 4)+(i-1)]<-0
  # re-renumber trait ids if necessary
  params<-x[x>0 & !is.na(x)]
  params[1]<-1 # fix first parameter to be 1
  if(length(params)>0){
    for(k in 2:length(params)){
      params[k]<-ifelse(params[k-1]==params[k] | params[k-1]==params[k]+1, params[k], params[k-1]+1)
    }
  }
  x[x>0 & !is.na(x)]<-params
  return(x)
}


################################################################################
# SET UP ROOT PRIORS AND ROOT TYPE
################################################################################

# setting the root prior for states 1 to 4 (A,T,D,V) for musse and muhisse models.
# we can say with some confidence that we expect the root state to have been A (i.e. 1):
# this function will generate appropriate root priors, depending on the number of hidden states

set_root<-function(n.hidden, root) {
  tmp=c(0,0,0,0)
  p=1/n.hidden
  tmp[root]=p
  rep(tmp,n.hidden)
}
