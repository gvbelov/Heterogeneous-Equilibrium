using JuMP, Ipopt

function calc_Gibbs(m,k,np,nc,g,jj,A,b,ion)  # Gibbs energy min
x_mol=fill(1.e-3,k)
y_mol=fill(1.e-3,np)
model = Model(Ipopt.Optimizer)
# set_optimizer_attributes(model, "tol" => 1e-8, "max_iter" => 1000)
MOI.set(model, MOI.RawOptimizerAttribute("print_level"), 0) # suppress output for Ipopt
#println("x_mol[1]=$x_mol[1]")
@variable(model, x[i=1:k] >= 0, start = x_mol[i])
#@variable(model, xc[1:nc] >= 0, start = 0.0)
@variable(model, y[i=1:np] >= 0, start = y_mol[i])
#@variable(model, yc >= 0, start = nc*1.e-3)

@NLobjective(model, Min, sum(x[i]*g[i] for i in 1:nc)+sum(sum(x[i]*(log(x[i]) + g[i]) for i in jj[1,j]:jj[2,j]) - y[j]*log(y[j]) for j in 1:np))
#@NLobjective(model, Min, sum(x[i]*g[i] for i in 1:nc)+sum(sum(x[i]*(log(x[i]/y[j]) + g[i]) for i in jj[1,j]:jj[2,j]) for j in 1:np))

for j in 1:np
@constraint(model, sum(x[i] for i in jj[1,j]:jj[2,j]) == y[j])
end

@constraint(model, con, A'*x .== b)

JuMP.optimize!(model)

#@show objective_value(model)

for i in 1:k x_mol[i]=value(x[i]) end
for i in 1:np y_mol[i]=value(y[i]) end
sp=zeros(m)
for i in 1:m sp[i]=shadow_price(con[i]) end
if ion != 0 sp[m]=-sp[m] end 
return x_mol, y_mol, sp
end

function get_data(fname)
txt=""
open(fname) do f1
txt =readlines(f1)
end
return txt
end

txt=get_data("test1")
m,k,np,nc,ion = parse.(Int32,split(txt[1]))

species=String[]; g=Float64[]; Nji=zeros(0,m);
jj=Int.(zeros(2,0)); b=Float64[]; 
for i in 1:k
st=split(txt[i+1])
push!(species,st[1]);
push!(g,parse(Float64,st[2]))
tmp=parse.(Float64,st[3:m+2])
global Nji=vcat(Nji,tmp')
end

for i in 1:np
st=split(txt[i+1+k])
tmp=parse.(Int32,st)
global jj=hcat(jj,tmp)
end

for i in 1:m
st=txt[i+1+k+np]
push!(b,parse(Float64,st))
end

x_mol, y_mol, sp = calc_Gibbs(m,k,np,nc,g,jj,Nji,b, ion)
for i in 1:k
println(species[i],"\t", x_mol[i])
end
