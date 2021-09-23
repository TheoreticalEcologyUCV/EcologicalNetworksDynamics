# EcologicalNetworksDynamics

EcologicalNetworksDynamics is a set of programs developed in R to simulate the dynamics of ecological networks, which in turn respond to the intrinsic population, community, metapopulation and/or metacommunity dynamics. 

The population dynamics are based on Generalized Lotka-Volterra Equations; the population interactions include the possibility of conditional outcomes; the spatial component assumes migrations between populations and/or between communities. The networks can be generated randomly, or based on actual data or natural (environmental, social, etc.) input. The program includes the possibility to generate networks with topologies that represent different types of ecosystems and agro-ecosystems; theoretical transitions between them can be simulated by modifying some network's central features.

The development of this program started in 2015; more components and options have been added as more theoretical explorations have been confronted. 

Results using these programs, or parts of them, have been published in: Griffon D (2015, Doctoral Thesis); Griffon D and Hernandez MJ (2014, 2019); Griffon D and Rodríguez G (2017); Ramírez D (2020, Bachelor Thesis); Griffon D, Hernandez MJ and Ramírez D (2021).

Details:
Each run of any of these programs is considered a single simulation or experiment. E. g.: If you want the data for 1000 simulations for a specific set of conditions, then you need to run the script (with the changes of the selected parameters) 1000 times.

To change the number of edges of the network check the sections "ADDING EDGES", "REMOVING EDGES" and the function erdos.renyi.game for (LINK UP PATH, LINK DOWN PATH AND RANDOM NETWORK CENTRALITY, respectively)

For additional help contact: griffondiego@gmail.com or davidearamirezp@gmail.com
