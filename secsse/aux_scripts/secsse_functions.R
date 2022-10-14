################################################################################
################################################################################
################################################################################
################################################################################
### Function to build transition rate matrix with more flexibility

### to add still:
## way to have each hidden transition to be different

secsse_q<-function(n_states,n_hidden, starting_index=1,no_dual=TRUE, hidden_transitions_fixed=TRUE,CID=FALSE){
  
  ### build basic matrix
  q<-matrix(1:(n_states*n_hidden)^2,
            ncol=n_states*n_hidden,
            nrow=n_states*n_hidden)
  
  ### name rows and columns
  rownames(q)<-paste0(rep(1:n_states,n_hidden),
                      rep(LETTERS[1:n_hidden], each=n_states))
  colnames(q)<-rownames(q)
  
  ### set diagonal to NA, toggle on/off dual transitions, fix transitions between hidden states to be the same
  for(i in 1:nrow(q)){
    string.dist<-c(adist(colnames(q), rownames(q)[i]))
    num.dist<-as.numeric(str_extract(colnames(q), pattern ="\\d+"))/as.numeric(str_extract(rownames(q), pattern ="\\d+"))[i]
    if(hidden_transitions_fixed) { q[i,num.dist  == 1] <-999 }
    if(no_dual) { q[i,string.dist  == 2]<-0 }
    q[i,string.dist  == 0] <- NA 
  }
  
  ### make CID 
  if(CID) {
    for(k in seq(1, nrow(q), n_states)[-1]){
      q[k:(k+n_states-1),k:(k+n_states-1)]<-q[1:n_states,1:n_states]
    }
  }
  
  ## re-index rowwise, with the exception of 999, with always goes last
  x = c(t(q))
  unique_x = unique(x[!is.na(x) & x!=0])
  unique_x = c(unique_x[unique_x!=999],999) # position 999 on the end
  y = seq(from=starting_index, to=starting_index+length(unique_x)-1)
  names(y)<-unique_x
  y["0"]<-0
  q[]<-unname(y[as.character(c(q))])
  
  ### return matrix
  return(q)
}


################################################################################
################################################################################
################################################################################
################################################################################
## re-index q

reindex_q<-function(q, start_index=min(q[q>0], na.rm = T)) {
  x = c(t(q))
  unique_x = unique(x[!is.na(x) & x!=0])
  y = seq(from=start_index, to=start_index+length(unique_x)-1)
  names(y)<-unique_x
  y["0"]<-0
  q[]<-unname(y[as.character(c(q))])
  return(q)
}

##################################################################################
##################################################################################
##################################################################################
##################################################################################
## mask and reindex q
### this function takes a binary matrix indicating which transitions should be set to 0 and then masks a transition matrix to set these transitions to 0

mask_q<-function(q, mask, n_states){
  
  for(i in 1:(nrow(q)/nrow(mask))){
    coordinates<-(1:n_states)+(n_states*(i-1))
    q[coordinates,coordinates][mask==0]<-0
  }
  return(reindex_q(q))
}  


##################################################################################
##################################################################################
##################################################################################
##################################################################################
## function to generate random starting tries for lambda, mu and transition rates. Takes an object "params" (which is a list that contains the idparslist), and initiation values for lambda and mu.
## starting tries are sampled from a normal distribution of the log starting try. 

gen_starting_values <-
  function(params, initLamb=intGuessLambda, initMu=intGuessMu) {
    
    n_lambda = max(params$idparslist$lambdas)
    n_mu = max(params$idparslist$mu)-n_lambda
    n_trans = max(params$idparslist$Q, na.rm = T)-n_lambda-n_mu
    
    lambda_start = exp(rnorm(n_lambda, log(initLamb)))
    mu_start = ifelse(initMu==0, 0, exp(rnorm(n_mu, log(initMu))))
    trans_start = exp(rnorm(n_trans, log((initLamb/length(params$idparslist$lambdas)))))
    
    return(c(lambda_start, mu_start, trans_start))
  }

##################################################################################
##################################################################################
##################################################################################
##################################################################################
### looping function to run seccse_lm with different starting tries (in parallel using mclapply)
## takes the list of models to run, the model names, the number of starting tries and the numbers of cores to use. 

secsse_multi_starting_tries<-function(model,model_name, n_starts, n_cores){
  
  if(!dir.exists("starting_tries")) dir.create("starting_tries")
  
  mclapply(FUN=function(x){
    params<-model
    params$initparsopt<-gen_starting_values(model)
    tmp<-do.call(secsse_ml,args=params)
    saveRDS(tmp,file = paste0("starting_tries/",model_name,"_try",j,"_out.rds"))},
    X=1:n_starts,
    mc.cores=n_cores)
}