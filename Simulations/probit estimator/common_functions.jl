function DGP(sampsize,param,seed,disg)
    Random.seed!(seed)
    e=randn(Float64,sampsize)
    Sigma_chol=cholesky(0.90*[1.0 0.0; 0.0 1.0]+0.10*ones(2,2))
    gz=randn(Float64,(sampsize,2))
    z=5*(atan.(gz*Sigma_chol.L')*1/pi.+0.5)

    if disg=="mixturenormal"
        g=Random.rand(MixtureModel(Normal, [(-param[4], 1.0), (0.0, 1.0), (param[4], 1.0)]),sampsize) .+ param[3]
    elseif disg=="logistic"
        g=Random.rand(Logistic(0,param[4]+1.0),sampsize) .+ param[3]
    end

    u=(param[1] .+ param[2].*z[:,1] .+ e).*z[:,2].- g
    y= (u.>0.0)*1.0
    
    return y, z
end


ncdf(x)=(1.0+erf(x/sqrt(2.0)))/2.0  # Normal CDF

function ProbitL(y,z,theta)
    beta0=theta[1]
    beta1=theta[2]
    mu=theta[3]
    sigma2=exp(theta[4])
    f=0.0
    n=length(y)
    for i in 1:n
        t=((beta0+beta1*z[i,1])*z[i,2] - mu)/sqrt(sigma2+z[i,2]^2)
        y[i]==1.0 ? ll=log(ncdf(t)) : ll=log(1.0 - ncdf(t))
        f=f + ll
    end
    return f
end


# function oneSim(seed)
#     y,z=DGP(sampsize,param,seed,disz,disg)
#     func = TwiceDifferentiable(vars -> -ProbitL(y, z, vars), ones(4); autodiff=:forward);
#     opt = optimize(func, ones(4))
#     return Optim.minimizer(opt), Optim.minimum(opt)
# end

function oneSimJ(seed)
    y,z=DGP(sampsize,param,seed,disg)
    model = Model(KNITRO.Optimizer)
    set_optimizer_attribute(model,"outlev",0)              # Turning off the ouput
    register(model, :ncdf, 1, ncdf; autodiff = true)
    set_optimizer_attribute(model,"ms_enable",1)      # Multistart option
    @variable(model, theta[1:4])
    set_start_value.(theta[1:4], vcat(param[1:3],0.0))
    @NLobjective(model, Max, sum( y[i]*log(ncdf(((theta[1]+theta[2]*z[i,1])*z[i,2] - theta[3])/sqrt(exp(theta[4])+z[i,2]^2)))+(1.0 - y[i])*log(1.0 - ncdf(((theta[1]+theta[2]*z[i,1])*z[i,2] - theta[3])/sqrt(exp(theta[4])+z[i,2]^2))) for i = 1:length(y)))
    optimize!(model)

    return JuMP.value.(theta), JuMP.objective_value(model)
end

