# Common functions

function DGP(sampsize,param,seed,disz,disg)
    Random.seed!(seed)
    e=randn(Float64,sampsize)
    if disz=="normal"
        Sigma_chol=cholesky(0.90*[1.0 0.0; 0.0 1.0]+0.10*ones(2,2))
        gz=randn(Float64,(sampsize,2))
        z=5*(atan.(gz*Sigma_chol.L')*1/pi.+0.5)
    elseif disz=="uniform" 
        z=5*rand(Float64,(sampsize,2))
    end

    if disg=="normal"
        g=randn(sampsize) .+ param[3]
    elseif disg=="mixturenormal"
        g=Random.rand(MixtureModel(Normal, [(-param[4], 1.0), (0.0, 1.0), (param[4], 1.0)]),sampsize) .+ param[3]
    elseif disg=="logistic"
        g=Random.rand(Logistic(0,param[4]+1.0),sampsize) .+ param[3]
    end

    u=(param[1] .+ param[2].*z[:,1] .+ e).*z[:,2].- g
    y= (u.>0.0)*1.0
    
    return y, z
end

lcdf(x)=1.0/(1.0+exp(-x))           # Logistic CDF
npdf(x)=exp(-x^2/2.0)/sqrt(2.0*pi)  # Normal PDF

function LogitL(y,z,theta)
    beta0=theta[1]
    beta1=theta[2]
    mu=theta[3]
    sigma2=exp(theta[4])

    f=0.0
    n=length(y)
    for i in 1:n
        t=quadgk(e->lcdf(((beta0+beta1*z[i,1]+e)*z[i,2]-mu)/sigma2)*npdf(e),-20.0,20.0)[1]
        y[i]==1.0 ? ll=log(t) : ll=log(1.0 - t)
        f=f + ll
    end
    return f
end


function oneSim(seed)
    y,z=DGP(sampsize,param,seed,disz,disg)
    func2(vars) =-LogitL(y, z, vars);
    opt = optimize(func2, 0.1*ones(4))
    return Optim.minimizer(opt), Optim.minimum(opt)
end
