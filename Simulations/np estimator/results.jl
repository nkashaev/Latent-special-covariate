#Main file for simulations with full estimator
# tempdir1=@__DIR__
# rootdir=tempdir1[1:findfirst("HiddenMenus",tempdir1)[end]]
# cd(rootdir*"/Environment/Mac")
# using Pkg
# Pkg.activate(".")

using CSV, DataFrames
using Distributions, Random, LinearAlgebra, Plots
tempdir1=@__DIR__
rootdir=tempdir1[1:findfirst("Random-Coeff",tempdir1)[end]]
dir=rootdir*"/Simulations/np estimator/"
dirresults=dir*"Results/"
## ######################### Parameters #########################################
disz="normal"
disg="mixturenormal"

Sampsize=[1000,5000,10000]

TabMeannp=zeros(7,3)
TabMADnp=zeros(7,3)
for i in 1:3
    sampsize=Sampsize[i]
    for k in 1:6
        param=[-0.5,1.0,0.5,1.0*(k-1)]
        Bias=Matrix(CSV.read(dirresults*"bias_np_4_4_$(sampsize)_$(disz)_$(disg)_$(param[4]).csv", DataFrame))
        TabMeannp[k,i]=mean(Bias[:,1])
        TabMADnp[k,i]=mean(abs.(Bias[:,1]))
    end
    Bias=Matrix(CSV.read(dirresults*"bias_np_4_4_$(sampsize)_$(disz)_logistic_0.0.csv", DataFrame))
    TabMeannp[7,i]=mean(Bias[:,1])
    TabMADnp[7,i]=mean(abs.(Bias[:,1]))
end
CSV.write(rootdir*"/Simulations/tables/Table1.csv", DataFrame(TabMeannp',:auto))
CSV.write(rootdir*"/Simulations/tables/Table2.csv", DataFrame(TabMADnp',:auto))

dir=rootdir*"/Simulations/probit estimator/"
dirresults=dir*"Results/"

TabProbit=zeros(7,2)
sampsize=Sampsize[1]
for k in 1:6
    param=[-0.5,1.0,0.5,1.0*(k-1)]
    Bias=Matrix(CSV.read(dirresults*"biasKN_probit_$(sampsize)_$(disz)_$(disg)_$(param[4]).csv", DataFrame))
    TabProbit[k,1]=mean(Bias)
    TabProbit[k,2]=mean(abs.(Bias))
end
Bias=Matrix(CSV.read(dirresults*"biasKN_probit_$(sampsize)_$(disz)_logistic_0.0.csv", DataFrame))
TabProbit[7,1]=mean(Bias)
TabProbit[7,2]=mean(abs.(Bias))
CSV.write(rootdir*"/Simulations/tables/Table3.csv", DataFrame(TabProbit',:auto))

dir=rootdir*"/Simulations/logit estimator/"
dirresults=dir*"Results/"

TabLogit=zeros(7,2)
sampsize=Sampsize[1]
for k in 1:6
    param=[-0.5,1.0,0.5,1.0*(k-1)]
    Bias=Matrix(CSV.read(dirresults*"bias_logit_$(sampsize)_$(disz)_$(disg)_$(param[4]).csv", DataFrame))
    TabLogit[k,1]=mean(Bias)
    TabLogit[k,2]=mean(abs.(Bias))
end
Bias=Matrix(CSV.read(dirresults*"bias_logit_$(sampsize)_$(disz)_logistic_0.0.csv", DataFrame))
TabLogit[7,1]=mean(Bias)
TabLogit[7,2]=mean(abs.(Bias))
CSV.write(rootdir*"/Simulations/tables/Table4.csv", DataFrame(TabLogit',:auto))