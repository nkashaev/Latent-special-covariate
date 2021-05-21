# Functions used in sim1

# This function generates data
# sampsize -- number of observations
# param    -- Y=I(param_3+[param_2+param_1 z_1+e_1]z_2+param_4 g >0)
# seed     -- random seed
# y,z      -- Output: binary y and z=[z1 z2]
function dgp(sampsize,param,seed)
    rng1 = MersenneTwister(seed)
    rng2 = MersenneTwister(seed+2)
    Sigma_chol=cholesky(0.90*[1.0 0.0; 0.0 1.0]+0.10*ones(2,2)) # variance matrix of z1 and z2
    ze=randn(rng1,Float64,(sampsize,3))
    g=-log.(1.0./rand(rng2,sampsize) .- 1.0)
    z= atan.(ze[:,2:3]*Sigma_chol.L').*2/pi .+ 1.0
    u=param[3] .+ (param[2] .+ param[1].*z[:,1] .+ ze[:,1]).*z[:,2] .+ (param[4].*g)
    y= (u.>0.0)*1.0
    return y, z
end

