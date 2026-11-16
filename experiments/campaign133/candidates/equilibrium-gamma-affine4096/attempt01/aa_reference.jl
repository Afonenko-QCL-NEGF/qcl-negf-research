module ParentAAOracle
using LinearAlgebra
import ..EquilibriumGammaPilot
const Field=EquilibriumGammaPilot.Field
const inner=EquilibriumGammaPilot.inner
const wnorm=EquilibriumGammaPilot.wnorm
function parent_aa_trial(states::Vector{Field}, residuals::Vector{Field}, we, wk, cfg)
    x,f = last(states),last(residuals)
    beta,alpha = Float64(cfg["aa_beta"]),Float64(cfg["linear_alpha"])
    linear() = x .+ alpha.*f
    length(states) < 2 && return linear(),(kind="linear",rank=0,condition=1.0,l1=1.0,coefficients=Float64[],firstpair=0,R_diagonal=Float64[],initial_scales=Float64[])
    for firstpair in 1:length(states)-1
        m = length(states)-firstpair
        columns = Field[]
        R = zeros(Float64,m,m)
        initial_scales = Float64[]
        valid = true
        for j in 1:m
            k = firstpair+j-1
            v = residuals[k+1] .- residuals[k]
            initial = wnorm(v,we,wk)
            push!(initial_scales,initial)
            for pass in 1:2, i in 1:j-1
                projection = inner(columns[i],v,we,wk)
                R[i,j] += projection
                v .-= projection.*columns[i]
            end
            R[j,j] = wnorm(v,we,wk)
            if initial == 0 || R[j,j] <= 64eps(Float64)*initial
                valid = false; break
            end
            v ./= R[j,j]
            push!(columns,v)
        end
        valid || continue
        condition = cond(R)
        isfinite(condition) && condition <= cfg["qr_condition_max"] || continue
        rhs = [inner(q,f,we,wk) for q in columns]
        coeff = UpperTriangular(R)\rhs
        # Expansion in past (x_i+beta*f_i); constrained affine coefficients.
        affine = zeros(Float64,m+1); affine[end] = 1.0
        for j in 1:m
            affine[j] += coeff[j]; affine[j+1] -= coeff[j]
        end
        l1 = sum(abs,affine)
        isfinite(l1) && l1 <= cfg["affine_l1_max"] || continue
        trial = x .+ beta.*f
        for j in 1:m
            k = firstpair+j-1
            # Fused broadcast; no retained DeltaGamma/DeltaF arrays here.
            d = coeff[j]
            xnew,xold,fnew,fold = states[k+1],states[k],residuals[k+1],residuals[k]
            @. trial -= d*((xnew-xold)+beta*(fnew-fold))
        end
        return trial,(kind="AAII",rank=m,condition=condition,l1=l1,coefficients=coeff,
            firstpair=firstpair,R_diagonal=diag(R),initial_scales=initial_scales)
    end
    return linear(),(kind="linear_rank_reset",rank=0,condition=Inf,l1=1.0,coefficients=Float64[],firstpair=0,R_diagonal=Float64[],initial_scales=Float64[])
end

end
