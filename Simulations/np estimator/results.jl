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
## ######################### Functions #########################################
#include(dir*"common_functions.jl")
## ######################### Parameters #########################################
#disz="uniform"
disz="normal"
disg="logistic"
#disg="normal"
disg="mixturenormal"
t=1.0
param=[-0.5,1.0,0.5,t]

Sampsize=[1000,5000,10000]

TabMeannp=zeros(8,3)
TabMediannp=zeros(8,3)
for i in 1:3
    sampsize=Sampsize[i]
    for k in 1:6
        param=[-0.5,1.0,0.5,1.0*(k-1)]
        Bias=Matrix(CSV.read(dirresults*"bias_np_$(sampsize)_$(disz)_$(disg)_$(param[4]).csv", DataFrame))
        TabMeannp[k,i]=mean(Bias[:,1])
        TabMediannp[k,i]=median(Bias[:,1])
    end
    Bias=Matrix(CSV.read(dirresults*"bias_np_$(sampsize)_$(disz)_normal_0.0.csv", DataFrame))
    TabMeannp[7,i]=mean(Bias[:,1])
    TabMediannp[7,i]=median(Bias[:,1])
    Bias=Matrix(CSV.read(dirresults*"bias_np_$(sampsize)_$(disz)_logistic_0.0.csv", DataFrame))
    TabMeannp[8,i]=mean(Bias[:,1])
    TabMediannp[8,i]=median(Bias[:,1])
end

dir=rootdir*"/Simulations/probit estimator/"
dirresults=dir*"Results/"

TabMeanP=zeros(8,3)
TabMedianP=zeros(8,3)
for i in 1:3
    sampsize=Sampsize[i]
    for k in 1:6
        param=[-0.5,1.0,0.5,1.0*(k-1)]
        Bias=Matrix(CSV.read(dirresults*"bias_probit_$(sampsize)_$(disz)_$(disg)_$(param[4]).csv", DataFrame))
        TabMeanP[k,i]=mean(Bias)
        TabMedianP[k,i]=median(Bias)
    end
    Bias=Matrix(CSV.read(dirresults*"bias_probit_$(sampsize)_$(disz)_normal_0.0.csv", DataFrame))
    TabMeanP[7,i]=mean(Bias)
    TabMedianP[7,i]=median(Bias)
    Bias=Matrix(CSV.read(dirresults*"bias_probit_$(sampsize)_$(disz)_logistic_0.0.csv", DataFrame))
    TabMeanP[8,i]=mean(Bias)
    TabMedianP[8,i]=median(Bias)
end


