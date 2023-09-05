# Heterogeneous-Equilibrium
The code of two functions written in Julia programming language for calculating the equilibrium composition of multicomponent heterogeneous thermodynamic systems and several examples are presented. 
To run the examples JuMP and Ipopt packages should be installed.

Description of the arguments:
m - number of elements;
k - number of substances; 
np - number of solutions (mixtures); 
nc - number of pure condensed subatances;
g - array of Gibbs energies of substances at temperature T and pressure P;
jj
A - formula matrix, A[j,i] represents the number of atoms of the ith element in the jth species;
b - array of amounts of chemical elements that form the thermodynamic system;
ion = 1 if ions are present in the gas phase, else ion = 0.
