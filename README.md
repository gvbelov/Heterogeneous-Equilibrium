# Heterogeneous-Equilibrium
The code of two functions written in Julia programming language for calculating the equilibrium composition of multicomponent heterogeneous thermodynamic systems and several examples are presented. 
To run the examples JuMP and Ipopt packages should be installed.

Substances and the data are arranged as follows. At the beginning there are condensed substances that form one-component phases, then gaseous substances, followed by substances that are part of condensed solutions.

Description of the arguments:  
m - number of elements;
k - number of substances; 
np - number of solutions (mixtures); 
nc - number of pure condensed subatances;
g - array of Gibbs energies of substances at temperature T and pressure P;
jj - the matrix of indices of substances in phases-solutions, jj[1,1] = nc+1, jj[2,1] = jj[1,1] + number of gaseous substances, jj[1,2] = jj[2,1]+1; jj[2,2] = jj[2,1] + number of substances in mixturе 2,... jj[1,np]=jj[2,np-1]+1; jj[2,np] = jj[2,np-1] + number of substances in mixturе np;
A - formula matrix, A[j,i] represents the number of atoms of the ith element in the jth species;
b - array of amounts of chemical elements that form the thermodynamic system;
ion = 1 if ions are present in the gas phase, else ion = 0.
