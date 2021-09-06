using CSV, DataFrames
using LinearAlgebra
using Random
using Distributions, Statistics
using QuadGK, ForwardDiff
#using ForwardDiff
using Optim
## Directories
tempdir1=@__DIR__
rootdir=tempdir1[1:findfirst("Random-Coeff",tempdir1)[end]]
dir=rootdir*"/Application/parametric_estimation/"
dirresults=dir*"results/"
dirdata=rootdir*"/Application/data/"

## Functions
include(dir*"/parametric_functions.jl")

## Data
data=CSV.read(dirdata*"margarine_subset.csv", DataFrame)
id=data[:,1];
choice=data[:,2];
z21=data[:,3]; #Blue Bonnet
z22=data[:,4]; #Fleischmann’s
z23=data[:,5]; #House Brand
z24=data[:,6]; #Generic
z25=data[:,7]; #Shed Spread
z1=data[:,8];

# Choice data
choice=data[:,2]; #2=Blue Bonnet, 3=Fleischman's,
                  #4=House brand, 5=Generic brand, 7=Shed Spread
n=length(choice); # Sample size

## Estimating beta
# For estimation of beta I use y=5 (Generic) as the outside good
ydata=zeros(length(choice))
ydata[choice.==2].=1
ydata[choice.==3].=2
ydata[choice.==4].=3
ydata[choice.==5].=5
ydata[choice.==7].=4
#Setting Generic as an outside good
xdata=hcat(z1,z21-z24,z22-z24,z23-z24,z25-z24) # z24 -- price of Generic


func2(vars) = -LogitL(ydata, xdata, vars);
opt = optimize(func2, -0.1*ones(7))
sol=Optim.minimizer(opt)
V=ForwardDiff.hessian(func2, sol)./length(ydata)
V=V[1:3,1:3]
beta0_param=sol[1]*exp(-sol[3])
beta1_param=sol[2]*exp(-sol[3])
se_beta1_param=sqrt([0.0,exp(-sol[3]),-sol[2]*exp(-sol[3])]'*V*[0.0, exp(-sol[3]),-sol[2]*exp(-sol[3])]/length(ydata))
se_beta0_param=sqrt([exp(-sol[3]), 0.0, -sol[1]*exp(-sol[3])]'*V*[exp(-sol[3]), 0.0, -sol[1]*exp(-sol[3])]/length(ydata))
se=sqrt.(abs.(diag(ForwardDiff.hessian(func2, sol))./length(ydata)^2))

