##############################################################################################
##############################################################################################
##############################################################################################
### corHMM model builder

# This is an auxiliary script for `_run_corHMM_anura.R`. It builds and plots the models
# model names correspond to descriptions in manuscript.
##############################################################################################
##############################################################################################
##############################################################################################



##############################################################################################
library(corHMM)
library(qgraph)


# Set up rate matrices

q<-list()


#################################################################
#################################################################
#################################################################
# NULL models



## ER
q$er<-getStateMat4Dat(dat, model="ER")
q$er$model<-"ER"
q$er$rate.cat<-1
q$er$root.p="maddfitz"
q$er

## SYM
q$sym<-getStateMat4Dat(dat, model="SYM")
q$sym$model<-"SYM"
q$sym$rate.cat<-1
q$sym$root.p="maddfitz"
q$sym

## ARD
q$ard<-getStateMat4Dat(dat, model="ARD")
q$ard$model<-"ARD"
q$ard$rate.cat<-1
q$ard$root.p="maddfitz"
q$ard


#################################################################
#################################################################
#################################################################
# Radial models


## radial model ARD
q$radial_ard<-getStateMat4Dat(dat, model="ER")
q$radial_ard$rate.mat[q$radial_ard$rate.mat>0]<-0
q$radial_ard$rate.mat[1,2:5]<-c(1,2,3,4)
q$radial_ard$rate.mat[2:5,1]<-c(5,6,7,8)

q$radial_ard$model<-"ARD"
q$radial_ard$rate.cat<-1
q$radial_ard$root.p=c(1,0,0,0,0)
q$radial_ard

## radial model ER
q$radial_er<-q$radial_ard
q$radial_er$rate.mat[q$radial_er$rate.mat>0]<-1
q$radial_er
q$radial_er$model<-"ER"

## radial SYM
q$radial_sym<-getStateMat4Dat(dat, model="ER")
q$radial_sym$rate.mat[q$radial_sym$rate.mat>0]<-0
q$radial_sym$rate.mat[1,2:5]<-c(1,2,3,4)
q$radial_sym$rate.mat[2:5,1]<-c(1,2,3,4)

q$radial_sym$model<-"SYM"
q$radial_sym$rate.cat<-1
q$radial_sym$root.p=c(1,0,0,0,0)
q$radial_sym


#################################################################
#################################################################
#################################################################
# Classic terrestrial serial models


## serial ARD

q$serial_ard<-getStateMat4Dat(dat, model="ER")
q$serial_ard$rate.mat[q$serial_ard$rate.mat>0]<-0
q$serial_ard$rate.mat[1,3]<-1
q$serial_ard$rate.mat[3,4]<-2
q$serial_ard$rate.mat[4,2]<-3
q$serial_ard$rate.mat[2,5]<-4

q$serial_ard$rate.mat[3,1]<-5
q$serial_ard$rate.mat[4,3]<-6
q$serial_ard$rate.mat[2,4]<-7
q$serial_ard$rate.mat[5,2]<-8

q$serial_ard$model<-"ARD"
q$serial_ard$rate.cat<-1
q$serial_ard$root.p=c(1,0,0,0,0)
q$serial_ard

## serial  ER
q$serial_er<-q$serial_ard
q$serial_er$rate.mat[q$serial_er$rate.mat>0]<-1

q$serial_er$model<-"ER"
q$serial_er$rate.cat<-1
q$serial_er$root.p=c(1,0,0,0,0)
q$serial_er

## serial  SYM
q$serial_sym<-q$serial_ard
q$serial_sym$rate.mat[3,1]<-1
q$serial_sym$rate.mat[4,3]<-2
q$serial_sym$rate.mat[2,4]<-3
q$serial_sym$rate.mat[5,2]<-4

q$serial_sym$model<-"SYM"
q$serial_sym$rate.cat<-1
q$serial_sym$root.p=c(1,0,0,0,0)
q$serial_sym


#################################################################
#################################################################
#################################################################
# Complex terrestrial


## terrestrial v1 ARD

q$terrestrial_v1_ard<-getStateMat4Dat(dat, model="ER")
q$terrestrial_v1_ard$rate.mat[q$terrestrial_v1_ard$rate.mat>0]<-0
q$terrestrial_v1_ard$rate.mat[1,3]<-1
q$terrestrial_v1_ard$rate.mat[2,3]<-2
q$terrestrial_v1_ard$rate.mat[4,3]<-3

q$terrestrial_v1_ard$rate.mat[3,1]<-4
q$terrestrial_v1_ard$rate.mat[3,2]<-5
q$terrestrial_v1_ard$rate.mat[3,4]<-6

q$terrestrial_v1_ard$rate.mat[1,5]<-7
q$terrestrial_v1_ard$rate.mat[5,1]<-8

q$terrestrial_v1_ard$model<-"ARD"
q$terrestrial_v1_ard$rate.cat<-1
q$terrestrial_v1_ard$root.p=c(1,0,0,0,0)
q$terrestrial_v1_ard

## terrestrial v1 SYM

q$terrestrial_v1_sym<-getStateMat4Dat(dat, model="ER")
q$terrestrial_v1_sym$rate.mat[q$terrestrial_v1_sym$rate.mat>0]<-0
q$terrestrial_v1_sym$rate.mat[1,3]<-1
q$terrestrial_v1_sym$rate.mat[2,3]<-2
q$terrestrial_v1_sym$rate.mat[4,3]<-3
q$terrestrial_v1_sym$rate.mat[1,5]<-4


q$terrestrial_v1_sym$rate.mat[3,1]<-1
q$terrestrial_v1_sym$rate.mat[3,2]<-2
q$terrestrial_v1_sym$rate.mat[3,4]<-3
q$terrestrial_v1_sym$rate.mat[5,1]<-4

q$terrestrial_v1_sym$model<-"SYM"
q$terrestrial_v1_sym$rate.cat<-1
q$terrestrial_v1_sym$root.p=c(1,0,0,0,0)
q$terrestrial_v1_sym

## terrestrial v1 ER

q$terrestrial_v1_er<-getStateMat4Dat(dat, model="ER")
q$terrestrial_v1_er$rate.mat[q$terrestrial_v1_er$rate.mat>0]<-0
q$terrestrial_v1_er$rate.mat[1,3]<-1
q$terrestrial_v1_er$rate.mat[2,3]<-1
q$terrestrial_v1_er$rate.mat[4,3]<-1
q$terrestrial_v1_er$rate.mat[1,5]<-1


q$terrestrial_v1_er$rate.mat[3,1]<-1
q$terrestrial_v1_er$rate.mat[3,2]<-1
q$terrestrial_v1_er$rate.mat[3,4]<-1
q$terrestrial_v1_er$rate.mat[5,1]<-1


q$terrestrial_v1_er$model<-"ER"
q$terrestrial_v1_er$rate.cat<-1
q$terrestrial_v1_er$root.p=c(1,0,0,0,0)
q$terrestrial_v1_er



## terrestrial v2 ARD

q$terrestrial_v2_ard<-getStateMat4Dat(dat, model="ER")
q$terrestrial_v2_ard$rate.mat[q$terrestrial_v2_ard$rate.mat>0]<-0
q$terrestrial_v2_ard$rate.mat[1,3]<-1
q$terrestrial_v2_ard$rate.mat[2,3]<-2
q$terrestrial_v2_ard$rate.mat[4,3]<-3

q$terrestrial_v2_ard$rate.mat[3,1]<-4
q$terrestrial_v2_ard$rate.mat[3,2]<-5
q$terrestrial_v2_ard$rate.mat[3,4]<-6

q$terrestrial_v2_ard$rate.mat[3,5]<-7
q$terrestrial_v2_ard$rate.mat[5,3]<-8


q$terrestrial_v2_ard$model<-"ARD"
q$terrestrial_v2_ard$rate.cat<-1
q$terrestrial_v2_ard$root.p=c(1,0,0,0,0)
q$terrestrial_v2_ard

## terrestrial v2 SYM

q$terrestrial_v2_sym<-getStateMat4Dat(dat, model="ER")
q$terrestrial_v2_sym$rate.mat[q$terrestrial_v2_sym$rate.mat>0]<-0
q$terrestrial_v2_sym$rate.mat[1,3]<-1
q$terrestrial_v2_sym$rate.mat[2,3]<-2
q$terrestrial_v2_sym$rate.mat[4,3]<-3
q$terrestrial_v2_sym$rate.mat[5,3]<-4

q$terrestrial_v2_sym$rate.mat[3,1]<-1
q$terrestrial_v2_sym$rate.mat[3,2]<-2
q$terrestrial_v2_sym$rate.mat[3,4]<-3
q$terrestrial_v2_sym$rate.mat[3,5]<-4

q$terrestrial_v2_sym$model<-"SYM"
q$terrestrial_v2_sym$rate.cat<-1
q$terrestrial_v2_sym$root.p=c(1,0,0,0,0)
q$terrestrial_v2_sym

## terrestrial v2 ER

q$terrestrial_v2_er<-getStateMat4Dat(dat, model="ER")
q$terrestrial_v2_er$rate.mat[q$terrestrial_v2_er$rate.mat>0]<-0
q$terrestrial_v2_er$rate.mat[1,3]<-1
q$terrestrial_v2_er$rate.mat[2,3]<-1
q$terrestrial_v2_er$rate.mat[4,3]<-1
q$terrestrial_v2_er$rate.mat[5,3]<-1

q$terrestrial_v2_er$rate.mat[3,1]<-1
q$terrestrial_v2_er$rate.mat[3,2]<-1
q$terrestrial_v2_er$rate.mat[3,4]<-1
q$terrestrial_v2_er$rate.mat[3,5]<-1

q$terrestrial_v2_er$model<-"ER"
q$terrestrial_v2_er$rate.cat<-1
q$terrestrial_v2_er$root.p=c(1,0,0,0,0)
q$terrestrial_v2_er


# terrestrial v3

q$terrestrial_v3_ard<-getStateMat4Dat(dat, model="ER")
q$terrestrial_v3_ard$rate.mat[q$terrestrial_v3_ard$rate.mat>0]<-0
q$terrestrial_v3_ard$rate.mat[1,3]<-1
q$terrestrial_v3_ard$rate.mat[2,3]<-2
q$terrestrial_v3_ard$rate.mat[4,3]<-3
q$terrestrial_v3_ard$rate.mat[5,2]<-4

q$terrestrial_v3_ard$rate.mat[3,1]<-5
q$terrestrial_v3_ard$rate.mat[3,2]<-6
q$terrestrial_v3_ard$rate.mat[3,4]<-7
q$terrestrial_v3_ard$rate.mat[2,5]<-8

q$terrestrial_v3_ard$model<-"ARD"
q$terrestrial_v3_ard$rate.cat<-1
q$terrestrial_v3_ard$root.p=c(1,0,0,0,0)
q$terrestrial_v3_ard

## terrestrial v3 SYM

q$terrestrial_v3_sym<-getStateMat4Dat(dat, model="ER")
q$terrestrial_v3_sym$rate.mat[q$terrestrial_v3_sym$rate.mat>0]<-0
q$terrestrial_v3_sym$rate.mat[1,3]<-1
q$terrestrial_v3_sym$rate.mat[2,3]<-2
q$terrestrial_v3_sym$rate.mat[4,3]<-3
q$terrestrial_v3_sym$rate.mat[5,2]<-4

q$terrestrial_v3_sym$rate.mat[3,1]<-1
q$terrestrial_v3_sym$rate.mat[3,2]<-2
q$terrestrial_v3_sym$rate.mat[3,4]<-3
q$terrestrial_v3_sym$rate.mat[2,5]<-4

q$terrestrial_v3_sym$model<-"SYM"
q$terrestrial_v3_sym$rate.cat<-1
q$terrestrial_v3_sym$root.p=c(1,0,0,0,0)
q$terrestrial_v3_sym

## Terrestrial v3 ER

q$terrestrial_v3_er<-getStateMat4Dat(dat, model="ER")
q$terrestrial_v3_er$rate.mat[q$terrestrial_v3_er$rate.mat>0]<-0
q$terrestrial_v3_er$rate.mat[1,3]<-1
q$terrestrial_v3_er$rate.mat[2,3]<-1
q$terrestrial_v3_er$rate.mat[4,3]<-1
q$terrestrial_v3_er$rate.mat[5,2]<-1

q$terrestrial_v3_er$rate.mat[3,1]<-1
q$terrestrial_v3_er$rate.mat[3,2]<-1
q$terrestrial_v3_er$rate.mat[3,4]<-1
q$terrestrial_v3_er$rate.mat[2,5]<-1

q$terrestrial_v3_er$model<-"ER"
q$terrestrial_v3_er$rate.cat<-1
q$terrestrial_v3_er$root.p=c(1,0,0,0,0)
q$terrestrial_v3_er




#################################################################
#################################################################
#################################################################
#####################
# set all transition related to Viviparity and Terrestrial to be symmetric


## function for re_indexing to remove gaps in the index:

re_index<-function(x) {
  y<-1:length(unique(x))
  names(y)<-unique(x)
  return(as.numeric(y[as.character(x)]))
}


for(i in 1:length(q)){
  q_temp<-q[[i]]$rate.mat
  q_temp[5,]<- q_temp[,5]
  q_temp[4,]<- q_temp[,4]
  
  # re index
    q_temp[q_temp>0 & !is.na(q_temp)]<-re_index(q_temp[q_temp>0 & !is.na(q_temp)])
  
  # replace
    q[[i]]$rate.mat<-q_temp
}

#################################################################
#################################################################
#################################################################
#####################
# make hidden states equivalents


RateClassMat <- getRateCatMat(2) #
RateClassMat <- equateStateMatPars(RateClassMat, c(1, 2)) # fix transitions between rate classes to be the same
RateClassMat

hmm_q<-list()

for(i in 1:length(q)){
  hmm_q[[i]]<-q[[i]]
  hmm_q[[i]]$rate.mat<-getFullMat(list(q[[i]]$rate.mat, q[[i]]$rate.mat), RateClassMat)
  hmm_q[[i]]$rate.cat<-2
  names(hmm_q)[i]<-paste0("HMM_",names(q)[i])
}


## combine models

q<-c(q, hmm_q)

#########################
##########################


# plot all models

par(mfrow=c(6,3))
plot.order<-c(1,3,4,2,5)
  
  
for(i in 1:length(q)){
  
  if(q[[i]]$rate.cat==1){
    
    qgraph(q[[i]]$rate.mat[plot.order, plot.order],
           title=names(q)[i],
           #nodes:
           labels=q[[i]]$legend[plot.order],
           shape="rectangle",
           vsize=30,
           node.height=0.5,
           #label.cex=3,
           label.color=c("#24b9e9", "#009E73", "#BF9C00","#55FF7E", "#f6776f"),
           # edges:
           fade=F,
           esize=1.5,
           edge.labels=T,
           edge.label.cex=3,
           asize=10,
           posCol="black",
           layout="circular",
           directed=TRUE)
  }
  
  
}

par(mfrow=c(1,1))
