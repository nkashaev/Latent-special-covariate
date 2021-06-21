function derivatives(f)

    gradf(x12)=ForwardDiff.gradient(f,x12)

    df11(x1)=gradf(vcat(x1,x2))[1]
    df12(x2)=gradf(vcat(x1,x2))[2]

    df21(x1)=ForwardDiff.derivative(df11,x1)
    df31(x1)=ForwardDiff.derivative(df21,x1)
    dfcros(x1)=ForwardDiff.derivative(x1->gradf(vcat(x1,x2))[2],x1)

    return df11, df12, df21, df31, dfcros
end

# function derivatives2(sieve_func, deltahat)
#     L(t)=1.0/(1.0+exp(-t))
#     L1(t)=L(t)*(1.0-L(t))
#     L2(t)=L(t)*(1.0-L(t))*(1.0-2.0*L(t))
#     L3(t)=L(t)*(1.0-L(t))*(1.0-6.0*L(t)+6.0*L(t)^2)
#     p(x)=sum([sieve_func[j](x).*deltahat[j] for j=1:k])
#     eta(x)=L(p(x))
#     deta1(x)=L1(p(x))
#     gradf(x12)=ForwardDiff.gradient(f,x12)
#
#     df11(x1)=gradf(vcat(x1,x2))[1]
#     df12(x2)=gradf(vcat(x1,x2))[2]
#
#     df21(x1)=ForwardDiff.derivative(df11,x1)
#     df31(x1)=ForwardDiff.derivative(df21,x1)
#     dfcros(x1)=ForwardDiff.derivative(x1->gradf(vcat(x1,x2))[2],x1)
#
#     return df11, df12, df21, df31, dfcros
# end


#Array of functions
# function Sieve_func_gen(order,xdata)
#     b=Array{Function,2}(undef,order+1,order+1)
#     for j=1:(order+1), k=1:(order+1)
#         g(x)=(x[1]^(j-1)) * (x[2]^(k-1))
#         b[j,k]=g
#     end
#     sieve_func_temp=b[:]
#     samplsize=size(xdata,1)
#     V=zeros(length(sieve_func_temp),length(sieve_func_temp))
#     for k=1:length(sieve_func_temp), l=1:length(sieve_func_temp)
#     V[k,l]=mean(sieve_func_temp[k](xdata[i,:])*sieve_func_temp[l](xdata[i,:]) for i=1:samplesize)
#     end
#     T=cholesky(V)
#     A=T.L^(-1)
#     sieve_func=Array{Function}(undef,length(sieve_func_temp))
#     for j=1:length(sieve_func_temp)
#         temp(x)=dot(A[j,:],[sieve_func_temp[k](x) for k=1:length(sieve_func_temp)])
#         sieve_func[j]=temp
#     end
#
#     return sieve_func
# end
function Sieve_func_gen(order,xdata)
    sieve_func_temp=Array{Function}(undef,Int64((order+1)*(order+2)/2))
    l=1
    for j=1:(order+1)
        #global l
        for k=1:((order+1)-(j-1))
            tempfun(x)=(x[1]^(k-1)) * (x[2]^((j-1)))
            sieve_func_temp[l]=tempfun
            l=l+1
        end
    end
    return sieve_func_temp

    # samplsize=size(xdata,1)
    #
    # SF=[sieve_func_temp[j](xdata[i,:]) for i=1:samplsize, j=1:length(sieve_func_temp)]
    # V=Hermitian(SF'*SF/samplesize)
    # A=inv(cholesky(V).L)
    # sieve_func=Array{Function}(undef,length(sieve_func_temp))
    # for j=1:length(sieve_func_temp)
    #     temp(x)=dot(A[j,:],[sieve_func_temp[k](x) for k=1:length(sieve_func_temp)])
    #     sieve_func[j]=temp
    # end
    # return sieve_func
end

#delta
function EstCoef(sieve_func,ydata,xdata)
n=size(ydata,1)
k=length(sieve_func)
X=zeros(n,k)
for i=1:n, j=1:k
    X[i,j]=sieve_func[j](xdata[i,:])
end
delta=X\ydata
return condexpec(x)=sum([sieve_func[j](x).*delta[j] for j=1:k])
end

# function CondExpec1(sieve_func,ydata,xdata)
# #deltahat_ini=EstCoef(sieve_func,ydata.-0.5,xdata)
# deltahat_ini=zeros(length(sieve_func))
# F(t)=cdf(Logistic(),t)
# etatemp(x,delta)=sum([sieve_func[j](x).*delta[j] for j=1:length(sieve_func)])
# NLLS(delta)=-sum(ydata[i]*log(F(etatemp(xdata[i,:],delta)))+(1-ydata[i])*(1-log(F(etatemp(xdata[i,:],delta)))) for i=1:length(ydata))
# #NLLS(delta)=sum([(ydata[i]-F(etatemp(xdata[i,:],delta)))^2 for i=1:length(ydata)])
# result = optimize(NLLS, deltahat_ini)
# deltahat=Optim.minimizer(result)
# #return x -> etatemp(x,deltahat)
# return condexpec(x)=F(etatemp(x,deltahat))
# end

function CondExpec2(sieve_func,ydata,xdata)
n=length(ydata)
k=length(sieve_func)
delta=Variable(k)
SF=[sieve_func[j](xdata[i,:]) for i=1:n, j=1:k]
modvex=minimize(logisticloss(-ydata.*(SF*delta)))
solve!(modvex,ECOS.Optimizer(verbose=0))
#solve!(modvex,ECOSSolver())
return condexpec(x)=1.0/(1.0+exp(-sum([sieve_func[j](x).*delta.value[j] for j=1:k])))
end
