using Distributions

struct momentum

end

function ZigZag(N::Normal, ξ0 , T)
    t = 0; θ = 1; ξ = ξ0
    Ξ = [0, ξ0]
    while t<T
        τ = sqrt.(-2*N.σ*log(rand()) + max.(0 ,ξ - N.μ)^2) - θ(ξ-N.μ)
        τ, i0 = findmin(τ)
        t = t + τ
        push!(Ξ,[t, ξ + θ.*τ])
        θ[i0] = -θ[i0]
    end
    return(Ξ)
end
