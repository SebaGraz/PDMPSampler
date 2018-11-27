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
        τ = zeros(d)
        V = inv(N.Σ)
        while t<T
            b = θ.*(V*θ)
            a = θ.*(V*ξ)
            println(a)
            for i in 1:2
                if b[i] > 0
                    if a[i] < 0
                        τ[i] = sqrt(-log(rand())*2.0/b[i]) - a[i]/b[i]
                    else #a[i]>0
                        τ[i] = sqrt((a[i]/b[i])^2 - log(rand())*2.0/b[i]) - a[i]/b[i]
                    end
                elseif  b[i] == 0
                    if a[i] > 0
                        τ[i] = -log(rand())/a[i]
                    else #a[i] <= 0
                        τ[i] = Inf
                    end
                else #b[i] < 0
                    if a[i] <= 0
                        τ[i] = Inf
                    else
                        u = rand()
                        if -log(u) <= -a[i]^2/b[i] + a[i]^2/(2*b[i])
                            τ[i] = - sqrt((a[i]/b[i])^2 - log(u)*2.0/b[i]) - a[i]/b[i]
                        else
                            τ[i] = Inf
                        end
                    end
                end
            end
            τ1, i0 = findmin(τ)
            t = t + τ1
            ξ = ξ + θ*τ1
            push!(Ξ, ξ)
            push!(Γ, t)
            θ[i0] = -θ[i0]
        end
        return(Γ, Ξ)
    end
end  # module ZigZag
