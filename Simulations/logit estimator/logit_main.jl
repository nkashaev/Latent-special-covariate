using Pkg
Pkg.activate(".")

using Distributed
addprocs(7)
using CSV, DataFrames
@everywhere begin
    using Distributions, KNITRO, JuMP, SpecialFunctions, Random, LinearAlgebra, Optim, QuadGK
    ## ######################### Dir ################################################
    tempdir1=@__DIR__
    rootdir=tempdir1[1:findfirst("Random-Coeff",tempdir1)[end]]
    dir=rootdir*"/Simulations/logit estimator/"
    dirresults=dir*"Results/"
    ## ######################### Functions #########################################
    include(dir*"common_functions.jl")
end

@everywhere begin
    disz="normal"
    sampsize=1000 
    #Set disg to "logistic" or mixturenormal
    
    #disg="logistic"
    disg="mixturenormal"
    
    #Set t to 0.0, 1.0, 2.0, 3.0, 4.0, or 5.0
    t=0.0
    param=[-0.5,1.0,0.5,t]
end

@time begin Result=pmap(oneSim,1:1000); end
Bias=[Result[i][1][2] - param[2] for i=1:1000]
CSV.write(dirresults*"bias_logit_$(sampsize)_$(disz)_$(disg)_$(param[4]).csv", DataFrame(Bias',:auto))


