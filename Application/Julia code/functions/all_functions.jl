# This function returns a tensor product of Chebyshev polynomials
# of order1 and order 2. scale_ap and shift_ap are used to adjust
# for the actuall support of z
function Cheb2_func_gen_app(order1,order2,scale_ap,shift_ap)
    sieve_func=Array{Function}(undef,Int64((order2+1)^4*(order1+1)))
    T=Array{Function}(undef,maximum([order1,order2])+1)
    T[1]=x->1 # First polynomial
    T[2]=x->x # Second polynomial
    for k=1:(length(T)-2)
    T[k+2]= x->2.0 .*x.*T[k+1](x)-T[k](x) # Recursive formula for higher order polynomials
    end
    l=1
    # Taking tenzor products
    for j1=1:(order2+1)
        for j2=1:(order2+1)
            for j3=1:(order2+1)
                for j4=1:(order2+1)
                    for k=1:(order1+1)
                        sieve_functemp(x)=T[k](x[1]) .* T[j1](x[2]).* T[j2](x[3]).* T[j3](x[4]).* T[j4](x[5])
                        sieve_func[l]=x->sieve_functemp(scale_ap*x .-shift_ap)
                        l=l+1
                    end
                end
            end
        end
    end
    return sieve_func
end

# This function regresses ydata on polynomials.
function EstCoef_app(sieve_func,ydata,xdata)
n=size(ydata,1)
ds=length(sieve_func)
X=zeros(n,ds)
for i=1:n, j=1:ds
    X[i,j]=sieve_func[j](xdata[i,:])
end
delta=X\ydata
condexpec(x)=sum([sieve_func[j](x).*delta[j] for j=1:ds])
return condexpec, delta, X
end

# This function computes partial deriavatives of a given function f
function derivatives_app(f)

    gradf(x_all)=ForwardDiff.gradient(f,x_all)

    f1(x1)=gradf(vcat(x1,x2_all))[1]
    f2(x2_all)=dot(x2_all,gradf(vcat(x1,x2_all))[2:end])
    f2temp(x1)=dot(x2_all,gradf(vcat(x1,x2_all))[2:end])
    f11(x1)=ForwardDiff.derivative(f1,x1)
    f12(x1)=ForwardDiff.derivative(f2temp,x1)
    f111(x1)=ForwardDiff.derivative(f11,x1)

    return f1,f2,f11,f12,f111 # z1, z2, z1^2, z1z2, z1^3
end

# This function computes the sign of beta1
function betasign1(f1,xdata)
    signbeta=0.0
    # I need to call to global variables since f1 implicitly takes x2_all as
    # an argument. See the definition of f1 (line 48 of this file)
    for i=1:n
        global x1=xdata[i,1]
        global x2_all=xdata[i,2:end]
    signbeta=signbeta+f1(x1)/n
    end
return sign(signbeta)
end

# This function estimates beta
function beta_coeff_app(f1,f2,f11,f12,f111,xdata,signbeta1)
    num=0.0
    denom=0.0
    a=0.0
    b=0.0
    c=0.0
    # Again, I need to call to global variables since all derivative functions
    # implicitly take x2_all as an argument.
@inbounds for i=1:size(xdata,1)
    global x1=xdata[i,1]
    global x2_all=xdata[i,2:end]
    num=num+f111(x1)*f1(x1)-f11(x1)^2
    denom=denom+f12(x1)*f1(x1)-f2(x2_all)*f11(x1)-f1(x1)^2
    a=a+f2(x2_all)-x1*f1(x1)
    b=b+f11(x1)
    c=c+f1(x1)
end
beta_sq1=num/denom
if beta_sq1>0
    beta1=signbeta1*sqrt(beta_sq1)
    beta0=(a*beta1-b/beta1)/c
    return [beta1,beta0], denom/size(xdata,1),beta1*c/size(xdata,1)
else
return [10000,10000], 0 , 0 #In finite samples the estimator may not be positive
end
end

# This function compute the derivatives of polynomials.
function derivativesofsieves(sieve_func)
sfunc1=Array{Function}(undef,length(sieve_func))
sfunc11=Array{Function}(undef,length(sieve_func))
sfunc111=Array{Function}(undef,length(sieve_func))
sfunc2=Array{Function}(undef,length(sieve_func))

for i=1:length(sieve_func)
    gradftemp(x_all)=ForwardDiff.gradient(sieve_func[i],x_all)
    sfunc1[i]=x1->gradftemp(vcat(x1,x2_all))[1]
    sfunc2[i]=x1->dot(x2_all,gradftemp(vcat(x1,x2_all))[2:end])
    sfunc11[i]=x1->ForwardDiff.derivative(sfunc1[i],x1)
    sfunc111[i]=x1->ForwardDiff.derivative(sfunc11[i],x1)
end
return sfunc1, sfunc11, sfunc111, sfunc2
end

# THis function computes the asymptotic variance of betahat
function Asy_variance(sieve_func,eta, sfunc1, sfunc11, sfunc111, sfunc2,beta1,gamma, ydata,xdata)
n=size(xdata,1)
ds=length(sieve_func)
Sigmahat=zeros(ds,ds)
for j=1:ds, k=1:ds
    for i=1:n
        Sigmahat[j,k]=Sigmahat[j,k]+sieve_func[j](xdata[i,:])*sieve_func[k](xdata[i,:]).*(ydata[i]-eta(xdata[i,:])).^2/n
    end
end
da1=zeros(ds)
da2=zeros(ds)
Psi1=zeros(ds)
Psi11=zeros(ds)
Psi111=zeros(ds)
Psi2=zeros(ds)
for i=1:n
    global x1=xdata[i,1]
    global x2_all=xdata[i,2:end]
    for j=1:ds
        Psi111[j]=sfunc111[j](x1)
        Psi11[j]=sfunc11[j](x1)
        Psi1[j]=sfunc1[j](x1)
        Psi2[j]=sfunc2[j](x1)
    end
    da1=da1+2.0*(Psi111*Psi1'-Psi11*Psi11')*gamma./n
    da2=da2+(beta1^2*(Psi2-x1.*Psi1)-Psi11)./n
end
return hcat(da1,da2), Sigmahat
end
