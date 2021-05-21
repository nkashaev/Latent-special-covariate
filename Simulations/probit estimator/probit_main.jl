#Main file for simulations with full estimator
# tempdir1=@__DIR__
# rootdir=tempdir1[1:findfirst("HiddenMenus",tempdir1)[end]]
# cd(rootdir*"/Environment/Mac")
# using Pkg
# Pkg.activate(".")

using Distributed
addprocs(7)
using CSV, DataFrames
@everywhere begin
    using Distributions, KNITRO, JuMP, SpecialFunctions, Random, LinearAlgebra, Optim
    ## ######################### Dir ################################################
    tempdir1=@__DIR__
    rootdir=tempdir1[1:findfirst("Random-Coeff",tempdir1)[end]]
    dir=rootdir*"/Simulations/probit estimator/"
    dirresults=dir*"Results/"
    ## ######################### Functions #########################################
    include(dir*"common_functions.jl")
    ## ######################### Parameters #########################################
    #disz="uniform"
    disz="normal"
    disg="logistic"
    param=[-0.5,1,0.5,5.0]
    sampsize=10000
end

@time oneSim(8)
@time oneSimJ(2)
@time begin for i in 1:10
    oneSimJ(i)
end
end    
@time begin Result=pmap(oneSim,1:1000); end
@time begin ResultJ=pmap(oneSimJ,1:1000); end
a=[Result[i][1][2] - 1.0 for i=1:1000]
aJ=[ResultJ[i][1][2] - 1.0 for i=1:1000]
histogram(a)
histogram(aJ)
histogram(a[abs.(a).<10])
mean(a[abs.(a).<2])
histogram(aJ[abs.(aJ).<10])


oneSimJ(seed)
oneSim(seed)
