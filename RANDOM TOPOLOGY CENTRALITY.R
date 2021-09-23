### RANDOM TOPOLOGY CENTRALITY 

library("igraph") 
library("network")  
library("intergraph")  

# Number of nodes
A <- 16
# Number of links
B <- 105

Randon<-erdos.renyi.game(A, B, type = c("gnm"), directed = F,  loops = F) # Network creation 

RandonNetwork<- asNetwork (Randon) # Transforming the original network into an Igraph object 

plot(RandonNetwork) 

centr_degree(Randon, mode = c("total"), loops = TRUE, normalized = TRUE) # Calculating centrality 
