
# This is a general code, which allows to simulate the population dynamics of community networks and metacommunity networks.
# It allows all types of fixed ecological interactions to be included in the networks, as well as two types of conditional interactions (i/e., linear and quadratic).


## LINK DOWN PATH

# Here the case of a single community is presented, with all interactions of the conditional quadratic type 



#--------------------------------------------------- 



## BASIC DEFINITIONS 



# Interaction functions 

# linear corresponds to 1 

# cuadratic corresponds to 2 

# fixed negative corresponds to 3 

# fixed positive corresponds to 4 

# null corresponds to 5 



# Number of species, comunities and iterations (time) 

n=16 # Number of species, must be greater than 2 

c=1 # Number of comunities, must be greater than 1 

t=120 # Number of iterations 

tamaño=n*c # General matrix dimension 



# Migration rates 

MigraMin=0.00 # Max migration rate 

MigraMax=0.00 # Min migration rate 



# Dummy parameters 

tie=1 # dummy 

cero<-0 # dummy 

AlfaFijoN<--1 # negative interaction value 

AlfaFijoP<-1 # positive interaction value 





#--------------------------------------------------- 



## FIXED MATRICES 



Ks <- diag(round (1/ runif(tamaño, min = 40, max = 150),2)) # Matrix with parameter k 

Rs<- diag(round(runif(tamaño, min = 0.5, max = 0.5),2)) # Matrix with parameter r 

Bs<- round(runif(tamaño, min = 1.5, max = 4)) # parameter vector for interaction function 

Cs <- round(runif(tamaño, min = 0.2, max = 0.4),2) # parameter vector for interaction function 

Es<- round(runif(tamaño, min = 1.5, max = 4)) # parameter vector for interaction function 

Fs <- round(runif(tamaño, min = 0.2, max = 0.4),2) # parameter vector for interaction function 



#--------------------------------------------------- 



## SPECIES DENSITIES MATRIX 



Ns<- matrix(NA, nrow =tamaño, ncol = t + 1) # Densities matrix 

Ns[,1]<- round(runif(tamaño, min = 1, max = 10)) # assignment of initial densities 



#--------------------------------------------------- 



## INTERACTION NETWORK (TOPOLOGY) 



MatInt<- array(NA,dim=c(tamaño, n, c)) # This array allows to create the Interaction Matrix subcolumns  - Each subcolumn is in a different dimension of the array – Inside the subcolumn lays a matrix that specifies the interaction between species 

Mat <- array(NA,dim=c(n, n, c)) # This step is also necesary fot setting the Interaction Matrix structure 

VectorInteraciones<-c(0,2,5,5,5,2,2,5,5,5,5,5,5,5,5,5,0,0,0,0,2,0,5,0,5,2,5,5,5,5,5,5,5,2,0,5,5,5,2,5,5,5,2,5,2,2,5,5,5,2,5,0,5,2,2,5,5,2,2,5,2,5,5,5,5,2,5,5,0,2,2,5,2,2,5,2,2,2,5,5,2,2,5,2,2,0,2,2,2,2,2,2,5,2,2,2,2,5,2,2,2,2,0,2,2,2,2,2,5,5,2,2,5,5,5,5,5,2,2,0,5,5,2,5,0,5,5,5,5,5,5,5,0,2,2,5,0,2,0,2,5,5,5,5,5,2,5,0,0,2,0,5,2,0,0,2,2,5,0,2,5,5,2,2,2,2,2,2,2,2,0,2,2,2,2,2,5,5,5,5,0,2,2,5,2,2,0,0,5,2,2,2,5,5,2,2,2,5,5,2,5,5,0,5,0,2,2,2,5,2,5,5,5,2,5,5,5,2,2,2,2,0,5,5,5,5,5,5,5,2,2,5,5,2,2,2,2,5,0,5,5,5,5,5,5,2,2,5,5,2,0,2,2,5,5,0) # This vector represents the interactions between species and translates to the Ecological Network Topology of the Conuco 



Inter <- t(matrix(VectorInteraciones, nrow =n, ncol =n)) # Interaction Matrix creation 

NoInter<- matrix(0, nrow =n, ncol = n) # Null matrix that will fill all the other subcolumns 

AA <- Inter # Dummy variable to identify each interaction sub-matrix 

BB<- NoInter # Dummy variable to identify each null sub-matrix 

for (col in 1:c) { 
  
  for (f in 1:c){ 
    
    BB->Mat[,,f]  
    
    AA ->Mat[,,col]} 
  
  lista <- list(Mat) 
  
  unidas<- do.call(rbind, lista) 
  
  columna<- t(matrix(unidas, nrow =n, ncol =tamaño)) 
  
  MatInt[,,col]<- columna} # This for loop creates the Interaction Matrix as an array and combines different submatrices in the right order inside each dimension of the bigger array 

listaInt <- list(MatInt) # Transforming the array into a list 

unidasInt<- do.call(rbind, listaInt) # this allows to join the subcolumns 

MatrizInteraccion<- matrix(unidasInt, nrow =tamaño, ncol =tamaño) # Transforming the list into a matrix 



#-------------------------------------------------- 



## REMOVING EDGES 



# This is the part where the number of edges can be changed starting from the original topology - this must be done manually and step by step, introducing the amount of edges you want the network to have initally



numenlac <- 105 # This is the number of edges the the initial network is going to have – It must be a number between 105 and 15 



for (contadorenlaces in 1:10000) 
  
{ 
  
  q<-round(runif(1, max = 256, min = 1), 0) 
  
  
  
  w <- MatrizInteraccion[[q]] 
  
  
  
  if (w == 2) MatrizInteraccion[[q]] <- 0 else MatrizInteraccion[[q]] <- MatrizInteraccion[[q]] 
  
  ww <- sum (MatrizInteraccion == 2) 
  
  if (ww == numenlac) break 
   
} 



sum (MatrizInteraccion == 2) # Confirmation of the number of edges 



#--------------------------------------------------- 



## FUNCTIONS FOR DENSITY-DEPENDENT CONDITIONAL INTERACTIONS 



FuncionLineal<- function(b, c, N) {  
  
  Alf<-b-c*N 
  
  Alf} # Linear function 

FuncionCuadratica<- function(e, f, N) {  
  
  Alf<- ((e*N)-(N^2))/(1+(f*N^2)) 
  
  Alf} # Quadratic function 



#--------------------------------------------------- 



## ALPHA-ITERATIVE MATRIX 



As <- array(NA,dim=c(tamaño, tamaño, t + 1)) # This array is to save the alphas 

MatrizAlfa<- function(MatrizInteraccion, As, Bs, Cs, Es, Fs, Ns) { for (i in 1:tamaño) for (j in 1:tamaño)  
  
  if (MatrizInteraccion [i,j]==0) As [j,i, tie] <- 0 else 
    
    if (MatrizInteraccion [i,j]==1) As [j,i, tie] <- FuncionLineal (Bs[i],Cs[i],Ns[j,tie]) else 
      
      if (MatrizInteraccion [i,j]==2) As [j,i, tie] <- FuncionCuadratica (Es[i],Fs[i],Ns[j, tie]) else 
        
        if (MatrizInteraccion [i,j]==3) As [j,i, tie] <- AlfaFijoN else 
          
          if (MatrizInteraccion [i,j]==4) As [j,i, tie] <- AlfaFijoP else 
            
            if (MatrizInteraccion [i,j]==5) As [j,i, tie] <- cero 
            
            diag(As [,,tie]) <- c(-1) 
            
            As} # This function gives the alpha that determines the interaction between species dictated by the Interaction Topology  

#As<- MatrizAlfa(MatrizInteraccion, As, Bs, Cs, Es, Fs, Ns) # This is confirmation  

#As 



#--------------------------------------------------- 



## MIGRATION DYNAMICS 



NsAntes<- rep(NA, times=tamaño) 

DensiAntes<- function(Ns, Rs, Ks, As) { Dummy1<-Ns[,tie]+Rs%*%Ns[,tie]+Rs%*%Ks%*% diag(Ns[,tie])%*%As[,,tie] %*% Ns[,tie]; 

for (D1 in 1:tamaño) { 
  
  ifelse(Dummy1[D1]<0, NsAntes[D1]<-0, NsAntes[D1]<-Dummy1[D1]); 
  
  if (Dummy1[1]==0) NsAntes<-rep(0, times=tamaño) 
  
}; 

NsAntes} # Determines densities before migration 



#--------------------------------------------------- 



## MIGRATION MATRIX 



MatMig<- array(NA,dim=c(tamaño, n, c)) # This array works in a similar way as the one used for the Interaction Matrix. It crates subcolumns that hold matrices which represent migrations 

Mat <- array(NA,dim=c(n, n, c)) # This creates a set a of subcolumns that are in different dimensions 

for (col in 1:c) { 
  
  AAA<- diag(round(runif(n, min = MigraMin, max = MigraMax),2)); 
  
  Aneg<- AAA*(-1) 
  
  BBB<-AAA/(c-1); 
  
  for (f in 1:c){ 
    
    BBB->Mat[,,f]  
    
    Aneg ->Mat[,,col]} 
  
  lista <- list(Mat) 
  
  unidas<- do.call(rbind, lista) 
  
  columna<- t(matrix(unidas, nrow =n, ncol =tamaño)) 
  
  MatMig[,,col]<- columna} # This loop creates the entire Migration Matrix linking all the subcolumns in a certain order 

listaMig <- list(MatMig) # Transforming the array into a list 

unidasMig<- do.call(rbind, listaMig) # Joins all the list components 

MatrizMigracion<- matrix(unidasMig, nrow =tamaño, ncol =tamaño) # Transforming the list into a matrix 



#--------------------------------------------------- 



## ITERATIVE DYNAMIC 



Npremigra<- matrix(NA, nrow =tamaño, ncol = t + 1) # Dummy matrix to save the densities previous to the migrations 



for (tie in 1:t) {  
  
  As<- MatrizAlfa(MatrizInteraccion, As, Bs, Cs, Es, Fs, Ns); 
  
  Npremigra[,tie]<- DensiAntes(Ns, Rs, Ks, As); 
  
  Dummy2<- Npremigra[,tie]+ MatrizMigracion %*% Npremigra[,tie]; 
  
  for (D2 in 1:tamaño) { 
    
    ifelse(Dummy2[D2]<0, Ns[D2,tie+1]<-0, Ns[D2,tie+1]<-Dummy2[D2]); 
    
    if (Dummy2[1]==0) Ns[,tie+1]<-rep(0, times=tamaño) 
    
  }} # This loop simulates metacomunity dynamics 



#--------------------------------------------------- 



## POPULATION DYNAMICS PLOT 



# Densities 

# Axes and names 

matplot(0:t+1, t(Ns), xlab = "Tiempo", ylab = "Densidades", pch = 1: tamaño, type = "l", col = 1:tamaño, main="Densidades")#graficando 

# Legend 

Poblaciones <- rep("Pobla", times=tamaño) 

Numeros<-c(1:tamaño) 

#legend("topleft", paste(Poblaciones, Numeros), lty = 1: tamaño, col = 1:tamaño,title = "Poblaciones") 



#--------------------------------------------------- 



## ANALYSIS 



#--------------------------------------------------- 



library("network")  

library("intergraph")  

library("igraph") 

vecrespob <- rep(NA, c+3) # Vector that saves population related results 

vecresred <- rep(NA, 12) # Vector that saves network analysis results  

vecresint <- rep(NA, 2) # Vector that saves results related to interaction between species  



#--------------------------------------------------- 



# First community population dynamics plot 



matplot(0:t+1, t(Ns[1:16,0:t+1]), xlab = "Tiempo", ylab = "Densidades", pch = 1: tamaño, type = "l", col = 1:tamaño, main="Densidades de la Primera Comunidad") # This plots the population dynamics for the first community 

# Legend 

Poblaciones <- rep("Pobla", times=tamaño) 

Numeros<-c(1:n) 

#legend("topleft", paste(Poblaciones, Numeros), lty = 1: tamaño, col = 1:tamaño,title = "Poblaciones") 



#--------------------------------------------------- 



layout(matrix(1:2, nrow = 1))  # Dividing plot area in two parts 



#--------------------------------------------------- 



# POPULATION 



# Number of surviving species  



generacion<-t 

VivosT<-0 

VivosMax<-0 

VivosMin<-16 

for (comuni in 1:c){ 
  
  Vivos<-0 
  
  for (sobre in ((((n*comuni)-n)+1):(n*comuni))){  
    
    ifelse (Ns[sobre, generacion]>0, Vivos<- Vivos+1, Vivos<- Vivos)}  
  
  print (Vivos) 
  
  VivosT <- VivosT+Vivos 
  
  MediaVivos <-VivosT/c 
  
  if (Vivos> VivosMax) VivosMax <-Vivos 
  
  if (Vivos < VivosMin) VivosMin <-Vivos 
  
  vecrespob [[comuni]] <- Vivos 
  
} # This loop prints the number of survivors each community has – That means it only saves in Vivos the number of survivors for the last community – Also saves the number of surviving species (VivosT) and the community with maximum and minimum richness (VivosMax) 



# Survival percentage 



MEDIAsobrevivencia <-round (MediaVivos/n*100,0) 

# Percentage of survivors 

MEDIAsobrevivencia 

vecrespob [[c+1]] <- MEDIAsobrevivencia 



MAXIMOsobrevivencia <-round (VivosMax/n*100,0) 

# Max percentage of survivors 

MAXIMOsobrevivencia 

vecrespob [[c+2]] <- MAXIMOsobrevivencia 



MINIMOsobrevivencia <-round (VivosMin/n*100,0) 

# Min percentage of survivors 

MINIMOsobrevivencia 

vecrespob [[c+3]] <- MINIMOsobrevivencia 



#--------------------------------------------------- 



# NETWORK ANALYSIS 



# Network extraction 



# Complete network 

RedTotal<- function(As,Red){ 
  
  for (TotalX in ((((n*cont)-n)+1):(n*cont))) {for (TotalY in ((((n*cont)-n)+1):(n*cont))) {  
    
    if (As [TotalX, TotalY,tie]==0) Red[TotalX, TotalY, cont] <-0 else 
      
      Red[TotalX, TotalY, cont] <-1}}  
  
  diag(Red [,,cont]) <- c(0) 
  
  Vector1<- rep(NA, times=tamaño) 
  
  for (tama in 1:tamaño) {ifelse (Ns[tama, tie] >0, Vector1[tama]<-1, Vector1[tama]<-0)} 
  
  MatrizDiago<- matrix(0, nrow = tamaño, ncol = tamaño) 
  
  MatrizDiago<- diag(Vector1) 
  
  Red[,,cont] <-MatrizDiago%*%Red[,,cont] %*% MatrizDiago 
  
  Red[,,cont]} # This function extracts the network from the As matrix without the unlinked nodes 



#--------------------------------------------------- 



# Extracting the Initial Networks 



tie<-1 # This allows to take the initial networks before the simulation process 

RedIni<- array(0,dim=c(tamaño, tamaño, c)) 

for (cont in 1:c){ 
  
  RedIni[,,cont]<- RedTotal(As,RedIni)} # This loop creates the networks 



REDinicial<- network(RedIni[1:n,1:n,1]) # Transforming the object into a network 

plot(REDinicial, main="Red Inicial") # Plotting the graph 



# Extracting the Final Network 



tie<-t  

RedFin<- array(0,dim=c(tamaño, tamaño, c)) 

for (cont in 1:c){ 
  
  RedFin[,,cont]<- RedTotal(As,RedFin)} 



REDfinal<- network(RedFin[1:n,1:n,1]) 

plot(REDfinal, main="Red Final")  



#--------------------------------------------------- 



# Basic network analysis 



# Initial Network  



REDinicialIgraph <- asIgraph(REDinicial) # Transforming the original network into an Igraph object 



vecresred [[1]] <- transitivity(REDinicialIgraph , type="global")  # Transitivity calculation 



vecresred [[3]] <- reciprocity(REDinicialIgraph) # Reciprocity calculation 



vecresred [[5]] <- edge_density(REDinicialIgraph, loops=F) # Density calculation 



vecresred [[7]] <- mean_distance(REDinicialIgraph, directed=T) # Distance calculation 



vecresred [[9]] <- diameter(REDinicialIgraph, directed=F, weights=NA) # Diameter calculation 



centr_degree(REDinicialIgraph, mode = c("total"), loops = TRUE, normalized = TRUE) 



vecresred [[11]] <- centr_degree(REDinicialIgraph)$centralization 



plot(REDinicialIgraph, edge.color= "dark red", vertex.color="gray40",vertex.label.cex = 0.8 , vertex.label.color="black", vertex.label.font = 3, vertex.label.dist = 2, layout=layout.circle, main="Red Inicial", edge.arrow.size=0.2) # Customizing the graph 



#Red final 



REDfinalIgraph<- asIgraph(REDfinal)  



vecresred [[2]] <- transitivity(REDfinalIgraph, type="global")   



vecresred [[4]] <- reciprocity(REDfinalIgraph) 



vecresred [[6]] <- edge_density(REDfinalIgraph, loops=F)  



vecresred [[8]] <- mean_distance(REDfinalIgraph, directed=T) 



vecresred [[10]] <- diameter(REDfinalIgraph, directed=F, weights=NA)  



centr_degree(REDfinalIgraph, mode = c("total"), loops = TRUE, normalized = TRUE) 



vecresred [[12]] <- centr_degree(REDfinalIgraph)$centralization 



plot(REDfinalIgraph, edge.color= "dark red", vertex.color="gray40",vertex.label.cex = 0.8 , vertex.label.color="black", vertex.label.font = 3, vertex.label.dist = 2, layout=layout.circle , main="Red Final", edge.arrow.size=0.2)  



#--------------------------------------------------- 



## INTRACCIONES 



# Last 6 iterations interaction values 



# Alpha values in the last 6 iterations 

a<-round (As[1:n,1:n,t-5],2) 

#a 

b<-round (As[1:n,1:n,t-4],2) 

#b 

c<-round (As[1:n,1:n,t-3],2) 

#c 

d<-round (As[1:n,1:n,t-2],2) 

#d 

e<-round (As[1:n,1:n,t-1],2) 

#e 

f<-round (As[1:n,1:n,t],2) 

#f 



# Alphas average 



Promedio <- round((a+b+c+d+e+f)/6,2) 

Promedio  



# Alphas variance 

Varianza <- round((a-Promedio)^2 + (b-Promedio)^2  +(c-Promedio)^2  +(d-Promedio)^2 + (e-Promedio)^2  + (f-Promedio)^2 , 2) 



Varianza 



#--------------------------------------------------- 



# Percentage of positive interactions 



# Number of interactions 

NumeroInteracciones <- 0 

for (enX in 1:n) {for (enY in 1:n) {  
  
  if (Promedio [enX, enY] >0) NumeroInteracciones<- NumeroInteracciones +1 
  
  if (Promedio [enX, enY] <0) NumeroInteracciones<- NumeroInteracciones +1 
  
}} 



NumeroInteracciones 



# Number of positive interactions 

NumeroInteraccionesP <- 0 

for (enX in 1:n) {for (enY in 1:n) {  
  
  if (Promedio [enX, enY] >0) NumeroInteraccionesP<- NumeroInteraccionesP +1 
    
}} 



NumeroInteraccionesP 



# Percentage of positive interactions in the community 



vecresint [[1]] <- round((NumeroInteraccionesP*100)/NumeroInteracciones,2) 



#--------------------------------------------------- 



# Percentage of varying interactions 



NumeroInteraccionesOscilantes <- 0 

for (enX in 1:n) {for (enY in 1:n) {  
  
  if (Varianza [enX, enY] >0) NumeroInteraccionesOscilantes<- NumeroInteraccionesOscilantes +1 
   
}} 



NumeroInteraccionesOscilantes	 



vecresint [[2]] <- round((NumeroInteraccionesOscilantes	*100)/NumeroInteracciones,2) 



# This percentage indicates the number of interactions that vary during the last 6 iterations 



#--------------------------------------------------- 



# Results vectors 



vecrespob  

vecresred  

vecresint 



sum (MatrizInteraccion == 2) # Confirmation for the number of interactions (edges) in the Ecological Network 





