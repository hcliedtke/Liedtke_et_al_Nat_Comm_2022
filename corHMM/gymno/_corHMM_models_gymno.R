##############################################################################################
##############################################################################################
##############################################################################################
### corHMM model builder

# This is an auxiliary script for `_run_corHMM_gymno.R`. It builds and plots the models
# model names correspond to descriptions in manuscript.
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
library(corHMM)
library(qgraph)


# Set up rate matrices

q<-list()

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

## serial ARD
q$serial_ard<-q$ard
q$serial_ard$rate.mat[2,3]<-0
q$serial_ard$rate.mat[3,2]<-0
q$serial_ard$rate.mat[1,1:3]<-c(0,3,4)
q$serial_ard

## serial ER
q$serial_er<-q$serial_ard
q$serial_er$rate.mat[q$serial_er$rate.mat>0]<-1
q$serial_er$model="ER"
q$serial_er

## serial SYM
q$serial_sym<-q$serial_ard
q$serial_sym$rate.mat[1,1:3]<-q$serial_sym$rate.mat[1:3,1]
q$serial_sym$model="SYM"
q$serial_sym


## terrestrial_v2 ARD
q$terrestrial_v2_ard<-q$ard
q$terrestrial_v2_ard$rate.mat[1:3,1]<-c(0,1,0)
q$terrestrial_v2_ard$rate.mat[1:3,2]<-c(2,0,3)
q$terrestrial_v2_ard$rate.mat[1:3,3]<-c(0,4,0)
q$terrestrial_v2_ard

## terrestrial_v2 ER
q$terrestrial_v2_er<-q$terrestrial_v2_ard
q$terrestrial_v2_er$rate.mat[q$terrestrial_v2_er$rate.mat>0]<-1
q$terrestrial_v2_er
q$terrestrial_v2_er$model<-"ER"

## terrestrial_v2 SYM
q$terrestrial_v2_sym<-q$terrestrial_v2_ard
q$terrestrial_v2_sym$rate.mat[1,1:3]<-c(0,1,0)
q$terrestrial_v2_sym$rate.mat[2,1:3]<-c(1,0,2)
q$terrestrial_v2_sym$rate.mat[3,1:3]<-c(0,2,0)
q$terrestrial_v2_sym
q$terrestrial_v2_sym$model<-"SYM"


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

q<-append(q, hmm_q)



#########################
##########################

# plot all models


par(mfrow=c(3,3))
plot.order<-c(2,1,3)
  
  
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
           label.color=c("#009E73","#55FF7E", "#f6776f"),
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
