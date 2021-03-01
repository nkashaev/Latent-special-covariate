using CSV, DataFrames
using LinearAlgebra
#using Optim
using Convex, ECOS
using Random
using Distributions, Statistics
using ForwardDiff
########################### Dir ################################################
rootdir="/Users/SSC4044-iMac27/Dropbox/PSU/Research/Random COeff/Simulations/Julia"
dirfunc=rootdir*"/functions"
dirresults=rootdir*"/results"
########################### Functions ##########################################
include(dirfunc*"/dgp.jl")
include(dirfunc*"/sieve_functions.jl")
include(dirfunc*"/beta_coeff.jl")
########################### Parameters #########################################
### Generation
samplesize=500
param=[-5 10 0.1 0.1]
order=4
M=1000
Beta_all=zeros(M,2)

for seed=1:M
ydata,xdata=DGP(samplesize,param,seed)
sieve_func=Sieve_func_gen(order,xdata)
eta=CondExpec2(sieve_func,ydata,xdata)

df11, df12, df21, df31, dfcros=derivatives(eta)
println(seed)
betahat=beta_coeff(eta,df11, df12, df21, df31, dfcros,xdata)
Beta_all[seed,:]=betahat'
CSV.write(dirresults*"/simul_n_$samplesize._param_$param.csv", DataFrame(Beta_all), writeheader=false)
end




CSV.write(dirresults*"/simul_n_$samplesize._param_$param.csv", DataFrame(Beta_all), writeheader=false)
Beta_all_cl=Beta_all[Beta_all[:,1].<3000,:]
histogram(Beta_all_cl[:,2].-param[2])
mean(Beta_all_cl[:,2])
BB=CSV.read(dirresults*"/simul_n_1000._param_[-5.0 10.0 0.1 0.1].csv")

BB_cl=BB[BB[:,1].<3000,:]
histogram(BB_cl[:,1].-param[1])
mean(BB_cl[:,1])


mean(BB[:,2])

using Plots
xx1=[minimum(xdata[:,1])+k*(maximum(xdata[:,1])-minimum(xdata[:,1]))/100 for k=0:100]
xx2=[minimum(xdata[:,2])+k*(maximum(xdata[:,2])-minimum(xdata[:,2]))/100 for k=0:100]
yy=[eta([xx1[i],0.4]) for i=1:length(xx1)]
plot(xx1,yy)

zz=[eta(xdata[i,:]) for i=1:length(ydata)]
zz=[eta([xx1[i],xx2[i]]) for i=1:length(xx1)]

xxx1=xdata[:,1]; xxx2=xdata[:,2];
xx1=rand(10000).*(maximum(xdata[:,1])-minimum(xdata[:,1])).+minimum(xdata[:,1])
xx2=rand(10000).*(maximum(xdata[:,2])-minimum(xdata[:,2])).+minimum(xdata[:,2])
zz=[eta([xx1[i],xx2[i]]) for i=1:length(xx1)]
surface( xx1, xx2, zz, size=[800,480] )


n=100; xr= 10*rand(n); yr= 10*rand(n); zr= Float64[ sin(xr[i])+cos(yr[i]) for i=1:n];
surface( xr, yr, zr, size=[800,480] )

LL(beta)=-sum([ydata[i]*log(cdf(Normal(),beta[1]+beta[2].*xdata[i,1]))+(1-ydata[i])*log(cdf(Normal(),-beta[1]-beta[2].*xdata[i,1])) for i=1:5000])
result = optimize(LL, zeros(2), BFGS())
Optim.minimizer(result)
