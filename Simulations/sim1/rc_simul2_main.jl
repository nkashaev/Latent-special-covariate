using CSV, DataFrames
using LinearAlgebra
#using Optim
using Plots
#using Convex, ECOS
using Random
using Distributions, Statistics
using ForwardDiff
########################### Dir ################################################
rootdir="/Users/SSC4044-iMac27/Dropbox/PSU/Research/Random COeff/Simulations/Julia"
rootdir="/Users/nailkashaev/Dropbox/PSU/Research/Random COeff/Simulations/Julia"

dirfunc=rootdir*"/functions"
dirresults=rootdir*"/results"
########################### Functions ##########################################
include(dirfunc*"/dgp.jl")
include(dirfunc*"/sieve_functions2.jl")
include(dirfunc*"/beta_coeff.jl")
########################### Parameters #########################################
### Generation
samplesize=5000
param=[-5.0 10.0 0.1 0.1]
order=8
M=1000
M1=M/1000;
Beta_all=zeros(M,2)
sieve_func=Cheb2_func_gen(order)
@time for seed=1:M
ydata,xdata=DGP(samplesize,param,seed)
#sieve_func=Sieve_func_gen(order,xdata)

#eta=EstCoef(sieve_func,ydata,xdata)
eta=EstCoef3(sieve_func,ydata,xdata)
#eta=CondExpec2(sieve_func,ydata,xdata)

df11, df12, df21, df31, dfcros=derivatives(eta)
println(seed)
#xdatatrim=xdata[(0.9.>xdata[:,1].>0.1).*(0.9.>xdata[:,2].>0.1),:]
betahat=beta_coeff(eta,df11, df12, df21, df31, dfcros,xdata)
Beta_all[seed,:]=betahat'
#CSV.write(dirresults*"/simul4noinv_n_$samplesize.order_$order._param_$param.csv", DataFrame(Beta_all), writeheader=false)
#CSV.write(dirresults*"/simul3noinv_n_$samplesize.order_$order.M1_$M1._param_$param.csv", DataFrame(Beta_all), writeheader=false)
end
sum(Beta_all[:,1].<1000)
mean(Beta_all[Beta_all[:,1].<1000,2])
std(Beta_all[Beta_all[:,1].<1000,2])
histogram(Beta_all[Beta_all[:,2].<1000,2])
maximum(Beta_all[Beta_all[:,2].<1000,2])

Beta_all_sh=CSV.read(dirresults*"/simul2_n_20000.order_7._param_[-5.0 10.0 0.1 0.1].csv", header=0)
sum(Beta_all_sh[:,1].<1000)
mean(Beta_all_sh[Beta_all_sh[:,1].<1000,2])
stdm(Beta_all_sh[Beta_all_sh[:,1].<1000,2],2)
histogram(Beta_all_sh[Beta_all_sh[:,2].<100,2])
