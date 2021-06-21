function beta_coeff(f,df11, df12, df21, df31, dfcros,xdata)
    num=0.0
    denom=0.0
    a=0.0
    b=0.0
    c=0.0
    #DF11=zeros(size(xdata,1))
    #beta_sq2=0.0
@inbounds for i=1:size(xdata,1)
    #global num,denom,a,b
    global x1=xdata[i,1]
    global x2=xdata[i,2]
    num=num+df31(x1)*df11(x1)-df21(x1)^2
    denom=denom+x2*(dfcros(x1)*df11(x1)-df12(x2)*df21(x1))-df11(x1)^2
    #beta_sq2=beta_sq2+((df31(x1)*df11(x1)-df21(x1)^2)/(x2*(dfcros(x1)*df11(x1)-df12(x2)*df21(x1))-df11(x1)^2))/size(xdata,1)
    #num=num+df31(x1)/df11(x1)-(df21(x1)/df11(x1))^2
    #denom=denom+(x2*dfcros(x1))/df11(x1)-1-x2*df12(x2)*df21(x1)/(df11(x1))^2
    a=a+x2*df12(x2)-x1*df11(x1)
    b=b+df21(x1)
    c=c+df11(x1)
    #a=a+x2*df12(x2)/df11(x1)-x1
    #b=b+df21(x1)/df11(x1)
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
