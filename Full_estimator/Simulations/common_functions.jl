# Common functions

function DGP(sampsize,param,seed)
    rng = MersenneTwister(seed);
    Sigma_chol=cholesky(0.90*[1.0 0.0; 0.0 1.0]+0.10*ones(2,2))
    g12=randn(rng,Float64,(sampsize,4))
    z=atan.(g12[:,3:4]*Sigma_chol.L')*1/pi.+0.5
    u=param[3].+(param[1].+param[2].*z[:,1].+g12[:,1]).*z[:,2].-(param[4].*g12[:,2])
    #u=param[1].+param[2].*z[:,1].+g[:,1]
    y= (u.>0.0)*1.0
    #return y, hcat(ones(sampsize,1),z)
    return y, z
end

ncdf(x)=(1.0+erf(x/sqrt(2.0)))/2.0  # Normal CDF
lcdf(x)=1.0/(1.0+exp(-x))           # Logistic CDF
ucdf(x)=x                           # Uniform CDF
npdf(x)=exp(-x^2/2.0)/sqrt(2.0*pi)  # Normal PDF
lpdf(x)=exp(-x)/(1.0+exp(-x))^2     # Logistic PDF
updf(x)=1.0                         # Uniform PDF

# function inth(beta,gamma,X,degree,cdf,inttol)
#     d=size(X,1)
#     A=zeros(d)
#     for i in 1:size(X,1)
#         b=X[i,2]
#         a=beta[3] + (beta[1]+beta[2]*X[i,1])*b
#         t=quadgk(e->cdf(sum(gamma[k]*(a+e*b)^(k-1) for k in 1:degree+1))*npdf(e), -Inf, +Inf, rtol=inttol)[1]
#         A[i]=maximum([minimum([t,1.0]),0.0])
#     end
#     return A
# end

function intF(a,b,γ,cdf,inttol)
    t=quadgk(e->cdf(sum(γ[k]*(a+e*b)^(k-1) for k in 1:length(γ)))*npdf(e), -Inf, +Inf, rtol=inttol)[1]
    return maximum([minimum([t,1.0]),0.0])
end


# function intF(β,γ,x,degree,cdf,inttol)
#     b=x[2]
#     a=beta[3] + (β[1]+β[2]*x[1])*x[2]
#     t=quadgk(e->cdf(sum(γ[k]*(a+e*b)^(k-1) for k in 1:degree+1))*npdf(e), -Inf, +Inf, rtol=inttol)[1]
#     return maximum([minimum([t,1.0]),0.0])
# end

function intG(a,b,γ,pdf,inttol)
    dg=length(γ)
    t=zeros(dg+1)
    ftemp(e)=pdf(sum(γ[k]*(a+e*b)^(k-1) for k in 1:dg))*npdf(e)
    t[1]=quadgk(e->ftemp(e)*sum(γ[k]*(k-1)*(a+e*b)^(k-2) for k in 2:dg), -Inf, +Inf, rtol=inttol)[1]
    for i in 1:length(γ)
        t[1+i]=quadgk(e->ftemp(e)*(a+e*b)^(i-1), -Inf, +Inf, rtol=inttol)[1]
    end
    return t
end

# @time inth1(beta,gamma,X[1,:],degree,cdf,inttol);
# A=zeros(size(X,1))
# @time begin for i in 1:size(X,1) 
#     A[i]=inth1(beta,gamma,X[i,:],degree,cdf,inttol);
# end end

# function momN(order)
#     t=0.0
#     if order==0 
#         t=1.0 
#     elseif iseven(order) 
#         t=Float64(doublefactorial(order-1))
#     else
#         t=0.0
#     end
#     return t
# end


# function intL(beta,gamma,X,degree)
#     d=size(X,1)
#     A=zeros(d)
#     for i in 1:size(X,1)
#         b=X[i,2]
#         a=beta[3] + (beta[1]+beta[2]*X[i,1])*b
#         t=sum(gamma[k]*sum(binomial(k-1, l)*a^(k-1-l)*b^l*momN(l) for l in 0:k-1) for k in 1:degree+1)
#         A[i]=maximum([minimum([t,1.0]),0.0])
#     end
#     return A
# end

# function logL(beta,gamma,Y,X,degree,cdf,inttol)
#     if length(gamma)!= degree+1
#         return "dimensionality of gamma is wrong"
#     else
#         d=size(X,1)
#         f=0.0
#         #P=intL(beta,gamma,X,degree)
#         P=inth(beta,gamma,X,degree,cdf,inttol)
#         for i in 1:d
#             Y[i]==1.0 ? f=f+log(P[i]) : f=f+log(1.0-P[i])
#             # if 0<P[i]<=1
#             # Y[i]==1.0 ? f=f+log(P[i]) : f=f+log(1.0-P[i])
#             # else
#             #  println([param,gamma])
#             # end
#         end
#         #println(f)
#         return f
#     end
# end


# function pLL(paramlong::Vector, grad::Vector)
#     if length(grad) > 0
#     end
#     beta1=paramlong[1:3]
#     gamma1=paramlong[4:end]
#     return logL(beta1,gamma1,Y,X,degree,cdf,inttol)
# end


function callbackEvalF(kc, cb, evalRequest, evalResult, userParams)
    x = evalRequest.x
    β=x[1:3]
    γ=x[4:end]
    f=0.0
    for i in 1:size(X,1)
        b=X[i,2]
        a=β[3] + (β[1]+β[2]*X[i,1])*b
        Pi=intF(a,b,γ,cdf,inttol)
        t =intG(a,b,γ,pdf,inttol)
        if Y[i]==1.0
            f=f + log(Pi)
        else
            f=f + log(1.0-Pi)
        end
    end
    evalResult.obj[1] = -f
    return 0
end

function callbackEvalG!(kc, cb, evalRequest, evalResult, userParams)
    x = evalRequest.x
    β=x[1:3]
    γ=x[4:end]
    g=zeros(length(x))
    for i in 1:size(X,1)
        b=X[i,2]
        a=β[3] + (β[1]+β[2]*X[i,1])*b
        Pi=intF(a,b,γ,cdf,inttol)
        t =intG(a,b,γ,pdf,inttol)
        if Y[i]==1.0
            g=g + vcat(t[1]*b,t[1]*b*X[i,1],t[1],t[2:end])/Pi           
        else
            g=g - vcat(t[1]*b,t[1]*b*X[i,1],t[1],t[2:end])/(1.0-Pi)
        end
    end
    for i in 1:length(g)
        evalResult.objGrad[i] = -g[i]
    end
    return 0
end



function oneMC(seed)
    Y,X=DGP(sampsize,param,seed)
    #mean(Y)
    kc = KNITRO.KN_new()
    KNITRO.KN_add_vars(kc, degree+4) # number of variables
    #Upper and lower bounds on Q
    # KNITRO.KN_set_var_lobnds(kc, zeros(dx))
    # KNITRO.KN_set_var_upbnds(kc, ones(dx))
    #Objective
    cb = KNITRO.KN_add_objective_callback(kc, callbackEvalF)
    #Gradient
    KNITRO.KN_set_cb_grad(kc, cb, callbackEvalG!)
    #No output on the screen
    #KNITRO.KN_set_param(kc, KNITRO.KN_PARAM_OUTLEV, KNITRO.KN_OUTLEV_NONE)
    KNITRO.KN_set_param(kc, KNITRO.KN_PARAM_OUTLEV, KNITRO.KN_OUTLEV_SUMMARY)
    # Algorithm
    KNITRO.KN_set_param(kc, KNITRO.KN_PARAM_ALGORITHM, KNITRO.KN_ALG_ACT_SQP)
    KNITRO.KN_set_var_primal_init_values(kc, vcat(beta,0.01*ones(degree+1)))
    KNITRO.KN_solve(kc)
    return paramlong=KNITRO.KN_get_solution(kc)[3]

end

