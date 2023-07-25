# Ensemble inference of unobserved infections in networks using partial observations

Codes for Zhang, R., Tai, J., Pei, S. Ensemble inference of unobserved infections in networks using partial observations. PLOS Computational Biology (2023). Codes were programmed by Jilei Tai (taijilei@mail.dlut.edu.cn).


## 1. Outbreak simulation

Input the network structure (hamster_1.txt),the transmission rates of individuals (beta_.txt), the time length of the simulation, and the initial number of seeds. Run simulation.cpp to generate the simulated outbreak file (State.txt) and observation record (diagno.txt) as outputs.

### Function: 

simulation.cpp - Simulation of an outbreak in a real-world network and observation of the states of partial nodes.

## 2. Ensemble inference

Input the network structure (hamster_1.txt), the transmission rates of individuals (beta_.txt), observation record (diagno.txt), and give the appropriate range of parameters (e.g. the number of ensemble members, the upper and lower bounds of susceptible population, an appropriate range of transmission rates ,the time length of the simulation). Run infer1.cpp or infer2.cpp to get infection probability data (infer.txt) as the output.

### Function: 

infer1.cpp - the original ensemble inference algorithm where the likelihood for each node is computed separately.

infer2.cpp - the expedited ensemble inference algorithm where the likelihood for all nodes are computed simultaneously.

## 3. Modified DMP

Input the network structure (hamster_1.txt), the transmission rates of individuals (beta_.txt), observation record (diagno.txt), and give the appropriate range of parameters (e.g. the number of ensemble members, the upper and lower bounds of susceptible population, an appropriate range of transmission rates ,the time length of the simulation). Run DMP.cpp to get infection probability data (DMP.txt) as the output.

### Function:

DMP.cpp - the modified DMP algorithm to estimate infection probability.
