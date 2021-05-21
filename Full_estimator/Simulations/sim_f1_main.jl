#Main file for simulations with full estimator
tempdir1=@__DIR__
rootdir=tempdir1[1:findfirst("HiddenMenus",tempdir1)[end]]
cd(rootdir*"/Environment/Mac")
using Pkg
Pkg.activate(".")

using Distributed
addprocs(7)
using CSV, DataFrames
@everywhere begin
  using QuadGK#, SpecialFunctions
  using LinearAlgebra
  using Random
  #using Distributions, Statistics, StatsBase
  using Combinatorics
  using JuMP, KNITRO
## ######################### Dir ################################################
tempdir1=@__DIR__
rootdir=tempdir1[1:findfirst("Random-Coeff",tempdir1)[end]]
dir=rootdir*"/Full_estimator/Simulations/"
dirresults=dir*"Results/"
## ######################### Functions #########################################
include(dir*"common_functions.jl")
end
## ######################### Parameters #########################################
seed=10
param=[-0.5,1,0.5,1]
beta=param[1:3]
sampsize=2000
degree=2
inttol=1e-8
gamma=0.01*ones(degree+1)
cdf=ncdf
pdf=npdf
function 

function callbackEvalF(kc, cb, evalRequest, evalResult, userParams)
  x = evalRequest.x
  evalResult.obj[1] = profLike(x,[])
  return 0
end



@time logL(param,gamma,Y,X,degree,cdf,inttol);## Passing parameters for optimization to different procs
@everywhere begin
  tau=0.01; Pd=$Pd; Pyd=$Pyd; sampsize=$sampsize
end
## Optimization
println("Starting optimization")
numMC=1000
@time begin Outcome=pmap(oneMC,1:numMC); end
## Preparing the output of MC simulation
Ms1=DataFrame([Outcome[s][4] for s in 1:numMC],:auto)
Ms2=DataFrame([Outcome[s][8] for s in 1:numMC],:auto)
## Saving the results
CSV.write(dirresults*"Ms1_$(model)_$(sampsize).csv", Ms1)
CSV.write(dirresults*"Ms2_$(model)_$(sampsize).csv", Ms2)
