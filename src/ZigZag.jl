module ZigZag
    export MultiNormal, Normal, ZigZagGaussian
    struct MultiNormal
        μ::Array{Float64, 1}
        Σ::Array{Float64, 2}
    end
    struct Normal
        μ::Float64
        σ::Float64
    end
    function ZigZagGaussian(N::Normal, ξ0::Float64, T::Float64)
        t = 0.0 ; θ = 1; ξ = ξ0
        Ξ = [ξ0]
        Γ = [t]
        while (t<T)
            τ = sqrt.(-2*N.σ*log(rand()) + max(0, θ*(ξ - N.μ))^2) - θ*(ξ - N.μ)
            t = t + τ
            ξ = ξ + θ*τ
            push!(Ξ, ξ)
            push!(Γ, t)
            θ = -θ
        end
        return(Γ,Ξ)
    end
    function ZigZagGaussian(N::MultiNormal, ξ0::Array{Float64, 1} , T::Float64)
        d = length(N.μ)
        t = 0.0 ; θ = ones(d); ξ = ξ0
        Ξ = [ξ]
        Γ = [t]
        while t<T
            τ = sqrt.(-2*N.Σ*log.(rand(d)) + max.(0 ,θ.*(ξ - N.μ)).^2) - θ.*(ξ - N.μ)
            τ, i0 = findmin(τ)
            t = t + τ
            ξ = ξ + θ.*τ
            push!(Ξ, ξ)
            push!(Γ, t)
            θ[i0] = -θ[i0]
        end
        return(Γ, Ξ)
    end
end  # module ZigZag
