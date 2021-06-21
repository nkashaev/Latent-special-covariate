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
include(dirfunc*"/sieve_functions3.jl")
include(dirfunc*"/beta_coeff3.jl")
########################### Parameters #########################################
### Generation
samplesize=5000
param=[-5.0 10.0 0.1 0.1]
order=6
M=1000
M1=M/1000;
Beta_all=zeros(M,2)
x1=1.0;x2=1.0;
sieve_func=Cheb2_func_gen(order)
# #Chebishev polynomials and derivatives
# T=Array{Function}(undef,order+1)
# T[1]=x->1
# T[2]=x->x
# for k=1:(order-1)
# T[k+2]= x->2.0 .*x.*T[k+1](x)-T[k](x)
# end
# #First and Second Derivatives
# DT=Array{Function}(undef,order+1)
# DDT=Array{Function}(undef,order+1)
# for k=1:(order+1)
#     DT[k]=x->ForwardDiff.derivative(T[k],x)
#     DDT[k]=x->ForwardDiff.derivative(DT[k],x)
# end
# #Sieve derivatives
# sieve_func1=Array{Function}(undef,Int64((order+1)*(order+2)/2))
# sieve_func2=Array{Function}(undef,Int64((order+1)*(order+2)/2))
# sieve_func11=Array{Function}(undef,Int64((order+1)*(order+2)/2))
# sieve_func12=Array{Function}(undef,Int64((order+1)*(order+2)/2))
# sieve_func22=Array{Function}(undef,Int64((order+1)*(order+2)/2))
# sieve_func112=Array{Function}(undef,Int64((order+1)*(order+2)/2))
# l=1
# for j=1:(order+1)
#     #global l
#     for k=1:((order+1)-(j-1))
#         sieve_func[l]=x->T[k](x[1]) .* T[j](x[2])
#         sieve_func1[l]=x->DT[k](x[1]) .* T[j](x[2])
#         sieve_func2[l]=x->T[k](x[1]) .* DT[j](x[2])
#         sieve_func11[l]=x->DDT[k](x[1]) .* T[j](x[2])
#         sieve_func12[l]=x->DT[k](x[1]) .* DT[j](x[2])
#         sieve_func22[l]=x->T[k](x[1]) .* DDT[j](x[2])
#         sieve_func112[l]=x->DDT[k](x[1]) .* DT[j](x[2])
#         l=l+1
#     end
# end

for seed=1:M
ydata,xdata=DGP(samplesize,param,seed)

eta=EstCoef3(sieve_func,ydata,xdata)
#eta=CondExpec2(sieve_func,ydata,xdata)

f1,f2,f11,f12,f22,f121=derivatives3(eta)
f1(x1);f2(x1);f11(x1);f12(x1);f22(x1);f121(x1);
println(seed)
#xdatatrim=xdata[(0.9.>xdata[:,1].>0.1).*(0.9.>xdata[:,2].>0.1),:]
betahat=beta_coeff3(f1,f2,f11,f12,f22,f121,xdata)
Beta_all[seed,:]=betahat'
#CSV.write(dirresults*"/simul5ChebAndZ2_n_$samplesize.order_$order._param_$param.csv", DataFrame(Beta_all), writeheader=false)
#CSV.write(dirresults*"/simul3noinv_n_$samplesize.order_$order.M1_$M1._param_$param.csv", DataFrame(Beta_all), writeheader=false)
end
sum(Beta_all[:,2].<1000)
mean(Beta_all[Beta_all[:,2].<1000,2])
std(Beta_all[Beta_all[:,1].<1000,2])
histogram(Beta_all[Beta_all[:,2].<100,2])
maximum(Beta_all[Beta_all[:,2].<1000,2])

Beta_all_sh=CSV.read(dirresults*"/simul2_n_20000.order_7._param_[-5.0 10.0 0.1 0.1].csv", header=0)
sum(Beta_all_sh[:,1].<1000)
mean(Beta_all_sh[Beta_all_sh[:,1].<1000,2])
stdm(Beta_all_sh[Beta_all_sh[:,1].<1000,2],2)
histogram(Beta_all_sh[Beta_all_sh[:,2].<100,2])
