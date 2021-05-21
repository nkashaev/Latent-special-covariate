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

## Computing ∫(a+be)^pϕ(e)de
function momentN(a,b,p)
    mN=0.0
    for k in 0:2:p
        mN=mN+a^(p-k) * b^(k) * 2.0^(-k/2) * factorial(k) ./ factorial(Int(k/2))
    end
    return mN
end


profObj=Model(KNITRO.Optimizer)
set_optimizer_attribute(profObj,"outlev",0)
@variable(profObj, alphaparam[1:hdegree+1])

function profLike(beta::Vector, grad::Vector)
    if length(grad) > 0
    end
    MN=zeros(N,hdegree+1)
    for k in 1:N, i in 1:hdegree+1
        MN[k,i]=momentN(beta[1]*z[k,2]+beta[2]*z[k,1]*z[k,2],z[k,2],i-1)
    end
    @NLobjective(profObj,Max, sum(y[k]*log(1/(1+exp(-sum(alphaparam[i]*MN[k,i] for i in 1:hdegree+1))))+(1.0-y[k])*(1.0-log(1/(1+exp(-sum(alphaparam[i]*MN[k,i] for i in 1:hdegree+1))))) for k=1:N))
    JuMP.optimize!(profObj)

    return objective_value(profObj)
end


function profLikeall(beta,profObj,y,z)
    N=
    MN=zeros(N,hdegree+1)
    for k in 1:N, i in 1:hdegree+1
        MN[k,i]=momentN(beta[1]*z[k,2]+beta[2]*z[k,1]*z[k,2],z[k,2],i-1)
    end
    @NLobjective(profObj,Max, sum(y[k]*log(1/(1+exp(-sum(alphaparam[i]*MN[k,i] for i in 1:hdegree+1))))+(1.0-y[k])*(1.0-log(1/(1+exp(-sum(alphaparam[i]*MN[k,i] for i in 1:hdegree+1))))) for k=1:N))
    JuMP.optimize!(profObj)

    return objective_value(profObj)
end




function callbackEvalF(kc, cb, evalRequest, evalResult, userParams)
    x = evalRequest.x
    evalResult.obj[1] = profLike(x,[])
    return 0
end
# Create a new Knitro solver instance.
kc = KNITRO.KN_new()
KNITRO.KN_add_vars(kc, 2) # number of variables
#Upper and lower bounds on Q
# KNITRO.KN_set_var_lobnds(kc, zeros(dx))
# KNITRO.KN_set_var_upbnds(kc, ones(dx))
#Objective
cb = KNITRO.KN_add_objective_callback(kc, callbackEvalF)
#No output on the screen
KNITRO.KN_set_param(kc, KNITRO.KN_PARAM_OUTLEV, KNITRO.KN_OUTLEV_NONE)
#KNITRO.KN_set_param(kc, KNITRO.KN_PARAM_OUTLEV, KNITRO.KN_OUTLEV_SUMMARY)
# Algorithm
KNITRO.KN_set_param(kc, KNITRO.KN_PARAM_ALGORITHM, KNITRO.KN_ALG_ACT_SQP)

function minimize_obj(j)
global param_pref=Pref_all[j,:,:]
#KNITRO.KN_set_var_primal_init_values(kc, PC_ini[:,j])
KNITRO.KN_set_var_primal_init_values(kc, 0.5*ones(1))
KNITRO.KN_solve(kc)
return KNITRO.KN_get_solution(kc)
end

@time nStatus, objSol, x, lambda_ =minimize_obj(2)


result1 = optimize(profLike, zeros(2), KNITRO(); autodiff = :forward)
Optim.minimum(result1)
betahat1=Optim.minimizer(result1) # Estimates of beta
profLike(beta)



function onerepl(seed)
    y,z =dgp(N,param,seed)
    profObj=Model(KNITRO.Optimizer)
    set_optimizer_attribute(profObj,"outlev",0)
    @variable(profObj, alphaparam[1:hdegree+1])
