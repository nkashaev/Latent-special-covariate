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

#polynomials
function polyn(order1,order2)
    pol1=Array{Function}(undef,order1+1)
    Dpol1=Array{Function}(undef,order1+1)
    Dpol11=Array{Function}(undef,order1+1)
    Dpol111=Array{Function}(undef,order1+1)
    for i in 1:length(pol1)
        pol1[i]=x->x[1]^(i-1)
        i<2 ? Dpol1[i]=x->0.0 : Dpol1[i]=x->(i-1)*x[1]^(i-2)
        i<3 ? Dpol11[i]=x->0.0 : Dpol11[i]=x->(i-1)*(i-2)*x[1]^(i-3)
        i<4 ? Dpol111[i]=x->0.0 : Dpol111[i]=x->(i-1)*(i-2)*(i-3)*x[1]^(i-4)
    end

    pol2=Array{Function}(undef,order2+1)
    Dpol2=Array{Function}(undef,order2+1)
    Dpol22=Array{Function}(undef,order2+1)
    for i in 1:length(pol2)
        pol2[i]=x->x[2]^(i-1)
        i<2 ? Dpol2[i]=x->0.0 : Dpol2[i]=x->(i-1)*x[2]^(i-2)
        i<3 ? Dpol22[i]=x->0.0 : Dpol22[i]=x->(i-1)*(i-2)*x[2]^(i-3)
    end

    Pol=Array{Function}(undef,(order1+1)*(order2+1))
    DPol1=Array{Function}(undef,(order1+1)*(order2+1))
    DPol2=Array{Function}(undef,(order1+1)*(order2+1))
    DPol12=Array{Function}(undef,(order1+1)*(order2+1))
    DPol11=Array{Function}(undef,(order1+1)*(order2+1))
    DPol22=Array{Function}(undef,(order1+1)*(order2+1))
    DPol111=Array{Function}(undef,(order1+1)*(order2+1))
    DPol112=Array{Function}(undef,(order1+1)*(order2+1))
    k=1
    for i in 1:length(pol1), j in 1:length(pol2)
        Pol[k]=x->pol1[i](x)*pol2[j](x)
        DPol1[k]=x->Dpol1[i](x)*pol2[j](x)
        DPol2[k]=x->pol1[i](x)*Dpol2[j](x)
        DPol12[k]=x->Dpol1[i](x)*Dpol2[j](x)
        DPol11[k]=x->Dpol11[i](x)*pol2[j](x)
        DPol22[k]=x->pol1[i](x)*Dpol22[j](x)
        DPol111[k]=x->Dpol111[i](x)*pol2[j](x)
        DPol112[k]=x->Dpol11[i](x)*Dpol2[j](x)
        k=k+1
    end
    return Pol, DPol1, DPol2, DPol12, DPol11, DPol22, DPol111, DPol112
end

function oneSim(seed)
    y,x=DGP(sampsize,param,seed,disz,disg)
    X=zeros(length(y),length(Pol))
    for i in 1:length(y), j in 1:length(Pol)
        X[i,j]=Pol[j](x[i,:])
    end
    gamma=X\y
    phat1(x)=sum(gamma[i]*DPol1[i](x) for i in 1:length(Pol))
    phat2(x)=sum(gamma[i]*DPol2[i](x) for i in 1:length(Pol))
    phat12(x)=sum(gamma[i]*DPol12[i](x) for i in 1:length(Pol))
    phat11(x)=sum(gamma[i]*DPol11[i](x) for i in 1:length(Pol))
    phat22(x)=sum(gamma[i]*DPol22[i](x) for i in 1:length(Pol))
    phat111(x)=sum(gamma[i]*DPol111[i](x) for i in 1:length(Pol))
    phat112(x)=sum(gamma[i]*DPol112[i](x) for i in 1:length(Pol))
    
    num1=sum(phat111(x[i,:])*phat1(x[i,:])-phat11(x[i,:])^2 for i in 1:size(x,1))
    denum1=sum(x[i,2]*phat12(x[i,:])*phat1(x[i,:])-x[i,2]*phat2(x[i,:])*phat11(x[i,:])-phat1(x[i,:])^2 for i in 1:size(x,1))

    num2=sum(phat112(x[i,:])*phat1(x[i,:])-phat12(x[i,:])*phat11(x[i,:]) for i in 1:size(x,1))
    denum2=sum(phat1(x[i,:])*phat2(x[i,:])-x[i,2]*phat22(x[i,:])*phat11(x[i,:])-x[i,2]*phat2(x[i,:])*phat12(x[i,:]) for i in 1:size(x,1))

    b1=sum( (phat111(x[i,:])*phat1(x[i,:])-phat11(x[i,:])^2)/(x[i,2]*phat12(x[i,:])*phat1(x[i,:])-x[i,2]*phat2(x[i,:])*phat11(x[i,:])-phat1(x[i,:])^2) for i in 1:size(x,1))/size(x,1)
    b2=sum( (phat112(x[i,:])*phat1(x[i,:])-phat12(x[i,:])*phat11(x[i,:]))/(phat1(x[i,:])*phat2(x[i,:])-x[i,2]*phat22(x[i,:])*phat11(x[i,:])-x[i,2]*phat2(x[i,:])*phat12(x[i,:])) for i in 1:size(x,1))/size(x,1)
    return sqrt(abs(num1/denum1)), sqrt(abs(num2/denum2)), sqrt(abs(b1)), sqrt(abs(b2))
end




