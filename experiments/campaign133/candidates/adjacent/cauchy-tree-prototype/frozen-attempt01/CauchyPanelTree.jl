module CauchyPanelTree
# WRITE-ONLY prototype: no Core dependency, no Green/Dyson, no density quadrature.
using LinearAlgebra
export build_tree, real_boundary, complex_cauchy, linear_gamma, raw_moments

struct Node
    lo::Int; hi::Int; center::Float64; radius::Float64
    children::Tuple{Int,Int}
    moments::Matrix{ComplexF64}
    reverse_moments::Matrix{ComplexF64}
    mass::Float64                       # convex Frobenius majorant, NOT trace
end
struct PanelTree
    x::Vector{Float64}
    gamma::Matrix{ComplexF64}            # (physical + zero ghosts, nb²)
    nb::Int
    nodes::Vector{Node}
    root::Int
    pmax::Int
    input_witness::NamedTuple
    geometry_witness::NamedTuple
end
finiteall(A) = all(z -> isfinite(real(z)) && isfinite(imag(z)), A)
function gauss_rule(n)
    beta = [j / sqrt(4j^2 - 1) for j in 1:n-1]
    f = eigen(SymTridiagonal(zeros(n), beta))
    return f.values, 2 .* abs2.(f.vectors[1, :])
end
function leaf_moments(x, g, lo, hi, c, r, p, tq, wq; reverse=false)
    mu = zeros(ComplexF64, p, size(g, 2))
    js = reverse ? (hi:-1:lo) : (lo:hi)
    ks = reverse ? (length(tq):-1:1) : (1:length(tq))
    for j in js, k in ks
        theta = (tq[k] + 1) / 2
        v = (1-theta) .* view(g,j,:) .+ theta .* view(g,j+1,:)
        t = ((x[j]-c) + theta*(x[j+1]-x[j])) / r
        a = (x[j+1]-x[j]) * wq[k] / (4pi)
        for q in 1:p
            view(mu,q,:) .+= a .* v
            a *= t
        end
    end
    return mu
end
function translated(children, nodes, c, r, p, C; reverse=false)
    mu = zeros(ComplexF64, p, C)
    for index in (reverse ? Base.reverse(children) : children)
        child = nodes[index]
        a, b = (child.center-c)/r, child.radius/r
        source = reverse ? child.reverse_moments : child.moments
        for q in 0:p-1
            js = reverse ? (q:-1:0) : (0:q)
            for j in js
                f = binomial(q,j) * a^(q-j) * b^j
                view(mu,q+1,:) .+= f .* view(source,j+1,:)
            end
        end
    end
    return mu
end
function build_tree(E, Gamma; ghosts, pmax=24, leaf_panels=8,
                    hermitian_atol=0.0, psd_atol=0.0)
    pmax isa Integer && 1 <= pmax <= 24 || throw(ArgumentError("pmax must be integer in 1:24"))
    leaf_panels isa Integer && leaf_panels >= 1 || throw(ArgumentError("positive integer leaf size required"))
    isfinite(hermitian_atol) && hermitian_atol >= 0 || error("Hermiticity tolerance")
    isfinite(psd_atol) && psd_atol >= 0 || error("PSD tolerance")
    eltype(E)===Float64 && eltype(Gamma)===ComplexF64 || error("prototype requires Float64 energy / ComplexF64 Gamma")
    all(v -> v isa Float64,ghosts) || error("Float64 ghosts required; no hidden coordinate conversion")
    ndims(Gamma)==3 && size(Gamma,2)==size(Gamma,3) || error("Gamma axes (E,a,b)")
    length(E)>=2 && size(Gamma,1)==length(E) || error("source dimensions")
    x = Float64[ghosts[1]; E; ghosts[2]]
    all(isfinite,x) && all(>(0),diff(x)) || error("strict finite source/ghost order")
    finiteall(Gamma) || error("nonfinite Gamma")
    nb = size(Gamma,2); g = zeros(ComplexF64,length(x),nb^2)
    g[2:end-1,:] .= reshape(Gamma,length(E),:)
    hdef, mineig = 0.0, Inf
    for i in 2:length(x)-1
        A = reshape(view(g,i,:),nb,nb)
        hdef = max(hdef,norm(A-A'))
        # Hermitian wrapper is used ONLY to measure a witness, never to alter g.
        mineig = min(mineig,eigmin(Hermitian((A+A')/2)))
    end
    hdef <= hermitian_atol && mineig >= -psd_atol || error("raw Hermitian/PSD witness fails")
    tq,wq = gauss_rule(cld(pmax+1,2)) # degree >= pmax: exact polynomial rule in exact arithmetic
    nodes = Node[]; depthmax = Ref(0)
    function make(lo,hi,depth)
        depthmax[] = max(depthmax[],depth)
        c = x[lo]+(x[hi+1]-x[lo])/2
        r = max(c-x[lo],x[hi+1]-c)
        isfinite(c) && isfinite(r) && r>0 || error("unusable interval geometry")
        if hi-lo+1 <= leaf_panels
            children = (0,0)
            mu = leaf_moments(x,g,lo,hi,c,r,pmax,tq,wq)
            rev = leaf_moments(x,g,lo,hi,c,r,pmax,tq,wq;reverse=true)
            mass = sum((x[j+1]-x[j])*(norm(view(g,j,:))+norm(view(g,j+1,:)))/(4pi) for j in lo:hi)
        else
            mid = (lo+hi)÷2
            children = (make(lo,mid,depth+1),make(mid+1,hi,depth+1))
            mu = translated(children,nodes,c,r,pmax,nb^2)
            rev = translated(children,nodes,c,r,pmax,nb^2;reverse=true)
            mass = sum(nodes[k].mass for k in children)
        end
        finiteall(mu) && finiteall(rev) && isfinite(mass) && mass>=0 || error("nonfinite node")
        push!(nodes,Node(lo,hi,c,r,children,mu,rev,mass)); return length(nodes)
    end
    root = make(1,length(x)-1,0)
    gaps = diff(x); ratio = maximum(max(gaps[i]/gaps[i+1],gaps[i+1]/gaps[i]) for i in 1:length(gaps)-1)
    isfinite(ratio) || error("unusable gap ratio")
    input = (;hermitian_defect=hdef,minimum_eigenvalue=mineig,hermitian_atol,psd_atol,
             representation_changed=false,formal_psd_certificate=false)
    geometry = (;max_depth=depthmax[],max_adjacent_gap_ratio=ratio,leaf_panels,
                balanced_by_panel_count=true,geometric_balance_claim=false)
    return PanelTree(x,g,nb,nodes,root,pmax,input,geometry)
end
function linear_gamma(T,E::Real)
    isfinite(E) || error("finite target required")
    x,g = T.x,T.gamma
    (E<x[1] || E>x[end]) && return zeros(ComplexF64,T.nb^2)
    E==x[end] && return copy(g[end,:])
    j = clamp(searchsortedlast(x,E),1,length(x)-1)
    theta = (E-x[j])/(x[j+1]-x[j])
    return (1-theta).*g[j,:] .+ theta.*g[j+1,:]
end
function panel_integral(T,j,z; real_pv=false)
    l,u = T.x[j],T.x[j+1]; width=u-l
    (real_pv && (z==l || z==u)) && error("touching panel must use grouped PV")
    # Remote direct leaves use high precision SCALAR coefficients to avoid cancellation.
    a,b = if abs(z-(l+width/2)) > 4width
        setprecision(BigFloat,192) do
            zz = real_pv ? BigFloat(z) : Complex{BigFloat}(z)
            ll,uu = BigFloat(l),BigFloat(u)
            L = real_pv ? log(abs(zz-ll))-log(abs(zz-uu)) : log(zz-ll)-log(zz-uu)
            B = (zz-ll)/(uu-ll)*L-1
            ComplexF64(L-B),ComplexF64(B)
        end
    else
        L = real_pv ? log(abs(z-l))-log(abs(z-u)) : log(z-l)-log(z-u)
        B = (z-l)/width*L-1
        ComplexF64(L-B),ComplexF64(B)
    end
    return (a .* T.gamma[j,:] .+ b .* T.gamma[j+1,:]) ./ (2pi)
end
function near_pv(T,E,lo,hi)
    lo>hi && return zeros(ComplexF64,T.nb^2),0.0
    gE = linear_gamma(T,E); pieces = Vector{ComplexF64}[]
    if any(v -> !iszero(v),gE)
        a,b = T.x[lo],T.x[hi+1]
        a<E<b || error("nonzero density requires two-sided PV region")
        push!(pieces,gE .* (log(E-a)-log(b-E)))
    end
    for j in lo:hi
        l,u=T.x[j],T.x[j+1]; B=(T.gamma[j+1,:]-T.gamma[j,:])/(u-l)
        v = -B*(u-l)
        if E!=l && E!=u
            v += (T.gamma[j,:]+B*(E-l)-gE) .* (log(abs(E-l))-log(abs(E-u)))
        end
        push!(pieces,v)
    end
    f = reduce(+,pieces); b = reduce(+,Base.reverse(pieces))
    return f/(2pi),norm(f-b)/(2pi)
end
function multipole(B,z,p)
    t=z-B.center; ratio=B.radius/t
    v=zeros(ComplexF64,size(B.moments,2))
    for q in p:-1:1; v = v.*ratio .+ B.moments[q,:]; end
    v ./= t
    alt=zeros(ComplexF64,length(v)); sensitivity=0.0
    for q in 1:p
        a=ratio^(q-1)/t
        alt .+= a.*B.moments[q,:]
        sensitivity += abs(a)*norm(B.moments[q,:]-B.reverse_moments[q,:])
    end
    return v,sensitivity,norm(v-alt)
end
function evaluate(T,z; real_pv, tau_far, reverse_order=false)
    isfinite(tau_far) && tau_far>0 || error("positive finite far tolerance required")
    (isfinite(real(z)) && isfinite(imag(z))) || error("finite target required")
    real_pv ? (z isa Real || error("PV requires real E")) : (imag(z)>0 || error("complex path requires Im z>0"))
    x=T.x; C=T.nb^2; plo,phi=1,0
    if real_pv && x[1]<=z<=x[end]
        j=clamp(searchsortedlast(x,z),1,length(x)-1)
        plo,phi=max(1,j-1),min(length(x)-1,j+1)
    end
    pieces=Vector{ComplexF64}[]; near_round=0.0
    if plo<=phi
        value,near_round=near_pv(T,z,plo,phi); push!(pieces,value)
    end
    stack=[T.root]; nf=0; nd=0; zero=0; far=0.0; moment=0.0; evalround=0.0
    visited=0; accepted=0; maxorder=0
    rootmass=T.nodes[T.root].mass; maxratio=0.0
    while !isempty(stack)
        index=pop!(stack); B=T.nodes[index]; visited+=1
        if plo<=B.lo && B.hi<=phi; continue; end
        overlap=!(B.hi<plo || phi<B.lo)
        if B.mass==0 && !overlap; zero+=B.hi-B.lo+1; continue; end
        d=abs(z-B.center); admissible=!overlap && d>B.radius && B.radius/d<=0.5
        p=1; trunc=Inf
        if admissible
            q=B.radius/d; allowance=tau_far*(B.mass/rootmass)
            while p<=T.pmax
                trunc=B.mass*q^p/(d-B.radius)
                trunc<=allowance && break
                p+=1
            end
        end
        if admissible && p<=T.pmax
            v,mw,ew=multipole(B,z,p); push!(pieces,v)
            far+=trunc; moment+=mw; evalround+=ew; nf+=B.hi-B.lo+1
            accepted+=1; maxorder=max(maxorder,p)
            maxratio=max(maxratio,B.radius/d)
        elseif B.children==(0,0)
            for j in (reverse_order ? (B.hi:-1:B.lo) : (B.lo:B.hi))
                plo<=j<=phi && continue
                push!(pieces,panel_integral(T,j,z;real_pv)); nd+=1
            end
        else
            append!(stack,reverse_order ? Base.reverse(B.children) : B.children)
        end
    end
    covered=max(0,phi-plo+1)+nf+nd+zero
    covered==length(x)-1 || error("near/far/direct partition does not cover source")
    value=isempty(pieces) ? zeros(ComplexF64,C) : reduce(+,pieces)
    reversevalue=isempty(pieces) ? copy(value) : reduce(+,Base.reverse(pieces))
    ledger=(;far_truncation_estimate=far,far_requested=tau_far,
        moment_reverse_sensitivity=moment,polynomial_order_sensitivity=evalround,
        near_sum_reverse_sensitivity=near_round,final_sum_reverse_sensitivity=norm(value-reversevalue),
        near_panels=max(0,phi-plo+1),far_panels=nf,direct_panels=nd,zero_panels=zero,
        visited_nodes=visited,accepted_far_nodes=accepted,maximum_far_order=maxorder,
        covered_panels=covered,max_far_radius_over_distance=maxratio,
        geometry_input_roundoff_bound="not_certified",near_log_roundoff_bound="not_certified",
        formal_whole_certificate=false,far_arithmetic_bound="not_outward_rounded")
    finiteall(value) || error("nonfinite evaluator result")
    return reshape(value,T.nb,T.nb),ledger
end
function real_boundary(T,E::Real; tau_far=1e-10,reverse_order=false)
    L,ledger=evaluate(T,E;real_pv=true,tau_far,reverse_order)
    return L-0.5im*reshape(linear_gamma(T,E),T.nb,T.nb),ledger
end
function complex_cauchy(T,z::Complex; tau_far=1e-10,reverse_order=false)
    return evaluate(T,z;real_pv=false,tau_far,reverse_order) # NO extra -iGamma/2
end
function raw_moments(T)
    S0=zeros(ComplexF64,T.nb^2); S1=copy(S0)
    for j in 1:length(T.x)-1
        l,u=T.x[j],T.x[j+1]; g,h=T.gamma[j,:],T.gamma[j+1,:]
        S0 .+= (u-l)/(4pi).*(g+h)
        S1 .+= (u-l)/(12pi).*((2l+u).*g+(l+2u).*h)
    end
    return reshape(S0,T.nb,T.nb),reshape(S1,T.nb,T.nb)
end
end
