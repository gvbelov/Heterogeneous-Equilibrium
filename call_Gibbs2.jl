#Calculate equilibrium composition of thermodynamic system 
#composed by 1 mole of FeOC at pressure 0.1 MPa, temperature 1000K
#condensed mixture is formed by Fe oxides: FeO(c), Fe2O3(c), Fe3O4(c)
#=========================================================================#
# JuMP and Ipopt packages should be installed
using JuMP, Ipopt

function calc_Gibbs(m,k,np,nc,g,jj,A,b,ion)  # Gibbs energy min
# set the initial values
x_mol=fill(1.e-3,k)
y_mol=fill(1.e-3,np)
model = Model(Ipopt.Optimizer)  # define the Model for optimization
# set_optimizer_attributes(model, "tol" => 1e-8, "max_iter" => 1000)
MOI.set(model, MOI.RawOptimizerAttribute("print_level"), 0) # suppress output for Ipopt
#println("x_mol[1]=$x_mol[1]")
@variable(model, x[i=1:k] >= 0, start = x_mol[i])
#@variable(model, xc[1:nc] >= 0, start = 0.0)
@variable(model, y[i=1:np] >= 0, start = y_mol[i])
#@variable(model, yc >= 0, start = nc*1.e-3)
#set the objective function
@NLobjective(model, Min, sum(x[i]*g[i] for i in 1:nc)+sum(sum(x[i]*(log(x[i]) + g[i]) for i in jj[1,j]:jj[2,j]) - y[j]*log(y[j]) for j in 1:np))
#@NLobjective(model, Min, sum(x[i]*g[i] for i in 1:nc)+sum(sum(x[i]*(log(x[i]/y[j]) + g[i]) for i in jj[1,j]:jj[2,j]) for j in 1:np))
# set constraints
for j in 1:np
  @constraint(model, sum(x[i] for i in jj[1,j]:jj[2,j]) == y[j])
end
# material balance constraints
@constraint(model, con, A'*x .== b)

JuMP.optimize!(model)

#@show objective_value(model)
# extract the results of calculations
for i in 1:k x_mol[i]=value(x[i]) end  # equilibrium composition
for i in 1:np y_mol[i]=value(y[i]) end  # phase moles
sp=zeros(m)
for i in 1:m sp[i]=shadow_price(con[i]) end  # Lagrange multipliers
if ion != 0 sp[m]=-sp[m] end # change the sign of Lagrange multiplier for electroneutrality constraint
return x_mol, y_mol, sp
end

function get_data(fname)
txt=""
open(fname) do f1
  txt =readlines(f1)
end
return txt
end

txt=get_data("test11.dat") # get the data from the test11.dat file
# first line of datafile contains
# m - number of elements; 
# k - number of substances; 
# np - number of solutions (mixtures); 
# nc - number of pure condensed subatances;
# ion = 1 if ions are present in the gas phase, else ion = 0
m,k,np,nc,ion = parse.(Int32,split(txt[1]))

species=String[]; g=Float64[]; Nji=zeros(0,m);
jj=Int.(zeros(2,0)); b=Float64[]; 
for i in 1:k
# read the chemical formula, Gibbs energy value at temparature T and pressure p, Nji - number of atoms of element j in substance i
  st=split(txt[i+1])
  push!(species,st[1]);
  push!(g,parse(Float64,st[2]))
  tmp=parse.(Float64,st[3:m+2])
  Nji=vcat(Nji,tmp')
end

for i in 1:np
# read the array of indices of substances in phases-solutions 
  st=split(txt[i+1+k])
  tmp=parse.(Int32,st)
  jj=hcat(jj,tmp)
end

for i in 1:m
# read the amounts of chemical elements b
  st=txt[i+1+k+np]
  push!(b,parse(Float64,st))
end

# call calc_Gibbs function
x_mol, y_mol, sp = calc_Gibbs(m,k,np,nc,g,jj,Nji,b, ion)
# print the equilibrium composition, moles
for i in 1:k
  println(species[i],"\t", x_mol[i])
end
