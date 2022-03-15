using Pkg
Pkg.activate(".")
using CSV, DataFrames
using LinearAlgebra
using Random
using Distributions, Statistics
## Directories
tempdir1=@__DIR__
rootdir=tempdir1[1:findfirst("Random-Coeff",tempdir1)[end]]
dir=rootdir*"/Application/semiparametric_estimation/"
dirresults=rootdir*"/Application/results/"
dirdata=rootdir*"/Application/data/"

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
choice=data[:,2]; #2=Blue Bonnet, 3=Fleischmann's, 
                  #4=House brand, 5=Generic brand, 7=Shed Spread
n=length(choice); # Sample size
# Income
IncomeSum=DataFrame(Average=round(mean(z1),digits=2),Median=median(z1),Min=minimum(z1),Max=maximum(z1))
#Table 2
Sumstat=describe(data[:,3:end-1])[:,1:5]
Table2=DataFrame(Brand=["Generic", "Blue Bonnet", "House Brand", "Shed Spread", "Fleischmann's"],Shares=[mean(data[:,2].==5), mean(data[:,2].==2), mean(data[:,2].==4), mean(data[:,2].==7), mean(data[:,2].==3)])
Table2[!, :MeanPrice]=Sumstat[[4,1,3,5,2],2]
Table2[!, :MedianPrice]=Sumstat[[4,1,3,5,2],4]
Table2[!, :MinPrice]=Sumstat[[4,1,3,5,2],3]
Table2[!, :MaxPrice]=Sumstat[[4,1,3,5,2],5]
Table2[:,2:end]=round.(Table2[:,2:end],digits=2)

## Saving the results
Table2=CSV.write(dirresults*"table2.csv", Table2)
