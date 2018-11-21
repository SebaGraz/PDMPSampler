using Distributions
struct momentum

end

function ZigZag(N::Normal, ξ0, T)
    d = length(N.μ)
    t = zeros(d) ; θ = 1; ξ = ξ0
    Ξ = [0, ξ0]
    while t<T
        τ = sqrt.(-2*N.Σ*log(rand(d) + max.(0 ,ξ - N.μ).^2) - θ(ξ-N.μ)
        t = t + τ
        push!(Ξ ,[t , ξ + θ*τ])
        θ = -θ
    end
    return(Ξ)
end
