# Pure postprocess: no Mats iteration, new physical Φ, μ fit or file access.
module MatsubaraMoments
using LinearAlgebra
const STATIC=(:acoustic,:impurity,:IFR,:alloy)
relative(a,b)=iszero(b) ? (iszero(a) ? 0.0 : Inf) : a/b
herm(A)=(A+A')/2
field(A)=reshape(ComplexF64.(A),1,size(A)...)
K(c,s,X)=Main.MatsubaraAll5.applyK(c,s,X)
weightednorm(c,X)=sqrt(sum(c.wk[k]*sum(abs2,@view X[1,k,:,:]) for k in axes(X,2)))
blocknorm(X)=maximum(opnorm(Matrix(X[1,k,:,:])) for k in axes(X,2))
function tails(N,beta)
    N>0 && beta>0 || error("tail domain")
    t=setprecision(256) do
        b=BigFloat(beta);piB=big(pi)
        s2=sum(inv(BigFloat(2n+1)^2) for n in 0:N-1)
        s4=sum(inv(BigFloat(2n+1)^4) for n in 0:N-1)
        (Float64(2/b*(b/piB)^2*(piB^2/8-s2)),Float64(2/b*(b/piB)^4*(piB^4/96-s4)))
    end
    scale=2/beta*(beta/pi)^4;lo=scale/(6*(2N+1)^3);hi=lo+scale/(2N+1)^4
    lo<=t[2]<=hi || error("t4 inconsistent with positive integral bounds")
    return (;t2=t[1],t4=t[2],t4_integral_lower=lo,t4_integral_upper=hi)
end
function moments(c)
    I4=zeros(ComplexF64,1,c.nk,c.nb,c.nb);h=similar(I4)
    for k in 1:c.nk
        I4[1,k,:,:]=Matrix{ComplexF64}(I,c.nb,c.nb)
        h[1,k,:,:]=c.full[k][c.centre,c.centre]-c.mu*I4[1,k,:,:]
    end
    KI=K(c,:LO,I4);S0=sum(K(c,s,I4) for s in STATIC)+(2*c.NB+1).*KI
    B0=sum(K(c,s,h) for s in STATIC)+K(c,:LO,(2*c.NB+1).*h+c.Omega.*I4)
    for k in 1:c.nk;B0[1,k,:,:].+=h[1,k,:,:]*S0[1,k,:,:]+S0[1,k,:,:]*h[1,k,:,:];end
    return (;h,KI,S0,B0)
end
function pauli(rho)
    return (;minimum=minimum(eigmin(Hermitian(herm(Matrix(rho[k,:,:])))) for k in axes(rho,1)),
        maximum=maximum(eigmax(Hermitian(herm(Matrix(rho[k,:,:])))) for k in axes(rho,1)),
        hermitian_defect=norm(rho-permutedims(conj(rho),(1,3,2))))
end
number(c,rho)=c.spin*sum(c.wk[k]*real(tr(rho[k,:,:])) for k in 1:c.nk)
function correct(c,rhoN,N;check=()->nothing)
    m=moments(c);t=tails(N,c.beta);a=2*c.Omega*t.t4;rhs=field(rhoN)+t.t4.*m.B0
    eta=a*blocknorm(m.KI);eta<1 || error("conditional Neumann bound unavailable; stop limited prototype")
    one=field(rhoN)+t.t4.*(m.B0-2*c.Omega.*K(c,:LO,field(rhoN)))
    d=length(rhs);L=Matrix{ComplexF64}(I,d,d);basis=zeros(ComplexF64,size(rhs))
    for j in 1:d
        check();fill!(basis,0);basis[j]=1;L[:,j].+=a.*vec(K(c,:LO,basis))
    end
    solved=reshape(L\vec(rhs),size(rhs));res=solved-rhs+a.*K(c,:LO,solved)
    rho=reshape(solved,size(rhoN));oneshot=reshape(one,size(rhoN))
    err=eta<1 ? blocknorm(res)/(1-eta) : Inf
    return (;rho,oneshot,feedback=rho-oneshot,oneshot_number=number(c,oneshot),oneshot_pauli=pauli(oneshot),S0=m.S0,B0=m.B0,correction=rho-rhoN,tails=t,
        linear_residual_absolute=norm(res),linear_residual_scale=norm(rhs),linear_residual_relative=relative(norm(res),norm(rhs)),
        linear_weighted_absolute=weightednorm(c,res),linear_weighted_scale=weightednorm(c,rhs),
        linear_weighted_relative=relative(weightednorm(c,res),weightednorm(c,rhs)),
        eta,conditional_linear_error_bound=err,linear_condition=cond(L),number=number(c,rho),pauli=pauli(rho),
        note="conditional linear moment-equation bound only; no finite-N physical remainder certificate")
end
"Classic M1 reported separately; correction is applied to exact-bare-subtracted rhoN."
function classic_M1(c,G,N)
    rho=zeros(ComplexF64,c.nk,c.nb,c.nb);m=moments(c);t=tails(N,c.beta)
    for k in 1:c.nk
        rho[k,:,:]=Matrix{ComplexF64}(I,c.nb,c.nb)/2-t.t2.*m.h[1,k,:,:]
        for n in 1:N;rho[k,:,:].+=(2/c.beta).*herm(Matrix(G[n,k,:,:]));end
    end
    return rho
end
"Independent auxiliary-pole oracle: exact moments/thermal residues, no K helper."
function scalar_oracle(beta;epsilon=0.1,b=-0.15,a=0.02)
    split=sqrt((epsilon-b)^2+4a);lp=(epsilon+b+split)/2;lm=(epsilon+b-split)/2
    wp=(lp-b)/split;wm=1-wp;f(x)=inv(1+exp(beta*x))
    exact=wp*f(lp)+wm*f(lm);B=a*(2epsilon+b);frompoles=wp*lp^3+wm*lm^3-epsilon^3
    rows=[]
    for N in (128,256)
        bare=f(epsilon)
        for n in 0:N-1
            z=im*(2n+1)*pi/beta;bare+=(2/beta)*real(wp/(z-lp)+wm/(z-lm)-inv(z-epsilon))
        end
        corrected=bare+tails(N,beta).t4*B
        push!(rows,(;N,exact,bare,corrected,bare_error=abs(bare-exact),M3_error=abs(corrected-exact)))
    end
    return (;B_formula=B,B_independent_poles=frompoles,moment_absolute_error=abs(B-frompoles),rows)
end
end
