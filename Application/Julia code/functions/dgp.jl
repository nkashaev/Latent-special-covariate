function DGP(sampsize,param,seed)
rng = MersenneTwister(seed);
Sigma_chol=cholesky(0.90*[1.0 0.0; 0.0 1.0]+0.10*ones(2,2))
g12=randn(rng,Float64,(sampsize,4))
z=atan.(g12[:,3:4]*Sigma_chol.L')*1/pi.+0.5
u=param[4].+(param[1].+param[2].*z[:,1].+g12[:,1]).*z[:,2].+(param[3].*g12[:,2])
#u=param[1].+param[2].*z[:,1].+g[:,1]
y= (u.>0.0)*1.0
#return y, hcat(ones(sampsize,1),z)
return y, z
end
