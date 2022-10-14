##############################################################################################
##############################################################################################
##############################################################################################
### Find transitions function

# This function facilitates finding branches where the nodes at each end are in a different state
# This function returns either a table, or (simmap=T) a painted tree (simmap format)
##############################################################################################
##############################################################################################
##############################################################################################

## library

library(ape)
library(phytools)

##

## start of function
find_transitions<-function(phy, tip_states, node_states, simmap=T, stem_prop=0.25) {
  
  # make data frame of all branches with start and end states
  branch_states<-data.frame(start_node=phy$edge[,1],
                            end_node=phy$edge[,2],
                            start_state=c(tip_states,node_states)[phy$edge[,1]],
                            end_state=c(tip_states,node_states)[phy$edge[,2]])
  
  # identify branches with shifts
  branch_states$shifts<-branch_states$start_state!=branch_states$end_state
  
  ## return simmap 
  if(simmap) {
    # reorder nodes to paint only downstream from root
    branch_states<-branch_states[order(branch_states$start_node),]
    # delete tip.states as this causes problems with paintSubTree()
    tree<-phy
    tree$tip.states<-NULL
    # apply 'base coat' of root state
    painted<-paintSubTree(tree = tree,
                          node = 1+length(phy$tip.label),
                          state = node_states[1],
                          anc.state = node_states[1],
                          stem = F)
    
    # loop through shifts to paint downstream clades
    for(i in which(branch_states$shifts)){
      painted<-paintSubTree(tree = painted,
                            node = branch_states$end_node[i],
                            state = branch_states$end_state[i],
                            stem = stem_prop)
    }
    
  }
  ifelse(simmap, return(painted), return(branch_states))
}