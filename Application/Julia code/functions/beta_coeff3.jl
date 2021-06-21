function beta_coeff3(f1,f2,f11,f12,f22,f121,xdata)
    num=0.0
    denom=0.0
    a=0.0
    b=0.0
    c=0.0
    #f1,f2,f11,f12,f22,f121=derivatives(f)
@inbounds for i=1:size(xdata,1)
    #global num,denom,a,b
    global x1=xdata[i,1]
    global x2=xdata[i,2]
    num=num+f121(x1)*f1(x1)-f11(x1)*f12(x1)
    denom=denom+x2*(f22(x1)*f1(x1)-f2(x1)*f12(x1))+f1(x1)*f2(x1)

    a=a+x2*f2(x1)-x1*f1(x1)
    b=b+f11(x1)
    c=c+f1(x1)
end
beta_sq1=num/denom
if beta_sq1>0
    beta1=sqrt(beta_sq1)
    beta0=(a*beta1-b/beta1)/c
    return [beta0,beta1]
else
return [10000,10000]
end
end
