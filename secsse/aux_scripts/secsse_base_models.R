################################################################################
################################################################################
#### SET UP SecSSE MODELS
################################################################################
################################################################################


### NOTE:
## 1. THESE MODELS HAVE TRANSITIONS TO AND FROM LIVE BEARING (VIVIPARITY) AND TERRESTRIAL SET TO BE SYMMETRICAL
## 2. All models are run once with full exitinction rate estimates and again with extinction rate parameters are set to be equal.

################################################################################
# SET UP MODEL LIST
models<-list()

################################################################################
### "MuSSE" MODELS - no hidden states
################################################################################

### "MuSSE"-style models where rates vary with trait states, but no hidden states

#musse FULL: turnover and extinction fraction different for each state, full Q matrix
models[["secsse_1_full"]]<-list(
  idparslist=id_paramPos(states,num_concealed_states=1),
  num_concealed_states=1
)
models[["secsse_1_full"]]$idparslist

######
#musse FULL, with mu FIXED to 0

models[["secsse_1_full_mu0"]]<-list(
  idparslist=id_paramPos(states,num_concealed_states=1),
  num_concealed_states=1
)

models[["secsse_1_full_mu0"]]$idparslist$mus[]<-rep(0, length(models[["secsse_1_full_mu0"]]$idparslist$mus))
models[["secsse_1_full_mu0"]]$idparslist$Q<-secsse_q(n_states=n_states, n_hidden=1, starting_index = max(unlist(models[["secsse_1_full_mu0"]]$idparslist[1:2]))+1)

models[["secsse_1_full_mu0"]]


################################################################################
### "CD" MODELS - constant diversification
################################################################################

### "MuSSE"-style models (no hidden states) where rates are fixed to be the same across trait states (i.e. constant diversfication)

#musse FIXED: turnover and extinction fraction fixed for each state, full Q matrix
models[["secsse_1_fixed"]]<-list(
  idparslist=id_paramPos(states, num_concealed_states=1),
  num_concealed_states=1
)

models[["secsse_1_fixed"]]$idparslist$lambdas[]<-rep(1, n_states)
models[["secsse_1_fixed"]]$idparslist$mus[]<-rep(2, n_states)
models[["secsse_1_fixed"]]$idparslist$Q<-secsse_q(n_states=n_states, n_hidden=1, starting_index = max(unlist(models[["secsse_1_fixed"]]$idparslist[1:2]))+1)


# musse FIXED with mu fixed to 0

models[["secsse_1_fixed_mu0"]]<-list(
  idparslist=id_paramPos(states, num_concealed_states=1),
  num_concealed_states=1
)
models[["secsse_1_fixed_mu0"]]$idparslist$lambdas[]<-rep(1, n_states)
models[["secsse_1_fixed_mu0"]]$idparslist$mus[]<-rep(0, length(models[["secsse_1_fixed_mu0"]]$idparslist$mus))
models[["secsse_1_fixed_mu0"]]$idparslist$Q<-secsse_q(n_states=n_states, n_hidden=1, starting_index = max(unlist(models[["secsse_1_fixed_mu0"]]$idparslist[1:2]))+1)


###############################################################################
### "MuHiSSE" MODELS 
################################################################################

### "MuHiSSE"-style models where rates vary with trait states and hidden states

#muhisse FULL 1 hidden state: turnover and extinction fraction different for each state, full Q matrix
models[["secsse_2_full"]]<-list(
  idparslist=id_paramPos(states, num_concealed_states=2),
  num_concealed_states=2
)

models[["secsse_2_full"]]$idparslist$Q<-secsse_q(n_states=n_states, n_hidden=2, starting_index = max(unlist(models[["secsse_2_full"]]$idparslist[1:2]))+1)

#####

#muhisse FULL 1 with fixed mu
models[["secsse_2_full_mu0"]]<-list(
  idparslist=id_paramPos(states, num_concealed_states=2),
  num_concealed_states=2
)

models[["secsse_2_full_mu0"]]$idparslist$mus[]<-rep(0, length(models[["secsse_2_full_mu0"]]$idparslist$mus))
models[["secsse_2_full_mu0"]]$idparslist$Q<-secsse_q(n_states=n_states, n_hidden=2, starting_index = max(unlist(models[["secsse_2_full_mu0"]]$idparslist[1:2]))+1)



###############################################################################
### "CID" MODELS 
################################################################################

### "CID"-style models [character independent diversification] where rates vary only across hidden states

#muhisse CID2: turnover and extinction fraction different observed and hidden trait, full Q matrix, 1 hidden state

models[["secsse_2_cid"]]<-list(
  idparslist=id_paramPos(states, num_concealed_states=2),
  num_concealed_states=2
)

models[["secsse_2_cid"]]$idparslist$lambdas[]<-rep(c(1,2), each=n_states) # fix speciation rates
models[["secsse_2_cid"]]$idparslist$mus[]<-rep(c(3,4), each=n_states) # fix extinction rates
models[["secsse_2_cid"]]$idparslist$Q<-secsse_q(n_states=n_states, n_hidden=2, starting_index = max(unlist(models[["secsse_2_cid"]]$idparslist[1:2]))+1,
                                   CID = TRUE)

####

#muhisse CID2 with mu fixed

models[["secsse_2_cid_mu0"]]<-list(
  idparslist=id_paramPos(states, num_concealed_states=2),
  num_concealed_states=2
)

models[["secsse_2_cid_mu0"]]$idparslist$lambdas[]<-rep(c(1,2), each=n_states) # fix speciation rates
models[["secsse_2_cid_mu0"]]$idparslist$mus[]<-rep(0, n_states*2) # fix extinction rates
models[["secsse_2_cid_mu0"]]$idparslist$Q<-secsse_q(n_states=n_states, n_hidden=2, starting_index = max(unlist(models[["secsse_2_cid_mu0"]]$idparslist[1:2]))+1,
                                   CID = TRUE)
