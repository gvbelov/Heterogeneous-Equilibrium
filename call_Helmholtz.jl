using JuMP, Ipopt

function calc_Helmholtz(m,k,np,nc,g,jj,A,b,ion) # Helmholtz energy min
model = Model(Ipopt.Optimizer)
# set_optimizer_attributes(model, "tol" => 1e-8, "max_iter" => 1000)
MOI.set(model, MOI.RawOptimizerAttribute("print_level"), 0) 
@variable(model, x[1:k] >= 0, start = 1.e-3)
@variable(model, y[1:np] >= 0, start = 1.e-3)
@NLobjective(model, Min, sum(x[i]*g[i] for i in 1:nc)+
sum(x[i]*(log(x[i]) + g[i]) for i in jj[1,1]:jj[2,1])+
sum(sum(x[i]*(log(x[i]) + g[i]) for i in jj[1,j]:jj[2,j]) - y[j]*log(y[j]) for j in 2:np))
for j in 1:np
@constraint(model, sum(x[i] for i in jj[1,j]:jj[2,j]) == y[j])
end
@constraint(model, con, A'*x .== b)
JuMP.optimize!(model)
xx=zeros(k)
for i in 1:k xx[i]=value(x[i]) end
yy=zeros(np)
for i in 1:np yy[i]=value(y[i]) end
sp=zeros(m)
for i in 1:m sp[i]=shadow_price(con[i]) end
if ion != 0 sp[m]=-sp[m] end 
return xx, yy, sp
end

function get_data(fname)
txt=""
open(fname) do f1
txt =readlines(f1)
end
return txt
end

txt=get_data("test1vi")
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

x_mol, y_mol, sp = calc_Helmholtz(m,k,np,nc,g,jj,Nji,b, ion)
for i in 1:k
println(species[i],"\t", x_mol[i])
end
