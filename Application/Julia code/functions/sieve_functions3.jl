function derivatives3(f)

    gradf(x12)=ForwardDiff.gradient(f,x12)
    hessf(x12)=ForwardDiff.hessian(f,x12)

    f1(x1)=gradf(vcat(x1,x2))[1]
    f2(x2)=gradf(vcat(x1,x2))[2]

    f11(x1)=ForwardDiff.derivative(f1,x1)
    f12(x1)=ForwardDiff.derivative(x1->gradf(vcat(x1,x2))[2],x1)
    f22(x2)=ForwardDiff.derivative(f2,x2)
    f121(x1)=ForwardDiff.derivative(f12,x1)

    # f11(x1)=hessf(vcat(x1,x2))[1,1]
    # f12(x1)=hessf(vcat(x1,x2))[1,2]
    # f22(x1)=hessf(vcat(x1,x2))[2,2]
    #
    # f121(x1)=ForwardDiff.derivative(f12,x1)


    return f1,f2,f11,f12,f22,f121
end

function Cheb2_func_gen(order)
    sieve_func=Array{Function}(undef,Int64((order+1)*(order+2)/2))
    T=Array{Function}(undef,order+1)
    T[1]=x->1
    T[2]=x->x
    for k=1:(order-1)
    T[k+2]= x->2.0 .*x.*T[k+1](x)-T[k](x)
    end
    l=1
    for j=1:(order+1)
        #global l
        for k=1:((order+1)-(j-1))
            sieve_func[l]=x->T[k](x[1]) .* T[j](x[2])
            l=l+1
        end
    end
    return sieve_func
end




# function Sieve_func_gen(order,xdata)
#     sieve_func_temp=Array{Function}(undef,Int64((order+1)*(order+2)/2))
#     l=1
#     for j=1:(order+1)
#         #global l
#         for k=1:((order+1)-(j-1))
#             tempfun(x)=(x[1]^(k-1)) * (x[2]^((j-1)))
#             sieve_func_temp[l]=tempfun
#             l=l+1
#         end
#     end
#     return sieve_func_temp
# end

#delta
function EstCoef3(sieve_func,ydata,xdata)
n=size(ydata,1)
k=length(sieve_func)
X=zeros(n,k)
for i=1:n, j=1:k
    X[i,j]=sieve_func[j](2.0.*xdata[i,:] .-1.0)
end
delta=X\ydata
return condexpec(x)=sum([sieve_func[j](2.0.*x.-1.0).*delta[j] for j=1:k])
end
