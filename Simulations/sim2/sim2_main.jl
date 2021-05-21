## Simulations 2
using Statistics
using DataFrames, CSV
addprocs(7)

@everywhere begin
  using Random
  using Distributions
  using ForwardDiff
  using Combinatorics
  using LinearAlgebra
  using JuMP
  using Gurobi
  using KNITRO
end
## Defining the file directories
tempdir1=@__DIR__
rootdir=tempdir1[1:findfirst("Random-Coeff",tempdir1)[end]]
simdir=rootdir*"/Simulations/sim2"
simresults=simdir*"/results"

## Functions
@everywhere include($(sim2dir)*"/functions_common_sim2.jl")

param=[1,1,1,1]
hdegree=3
N=10000
seed=1
y,z= dgp(N,param,seed)