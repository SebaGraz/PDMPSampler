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
    struct Skeleton
        t::Float64
        ξ::Array{Float64, 1}
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
        d = length(N.μ) ; t = 0.0 ; θ = ones(d);
        ξ = ξ0
        Ξ = [Skeleton(t,ξ)]
        τ = zeros(d)
        V = inv(N.Σ)
        while t<T
            b = θ.*(V*θ)
            a = θ.*(V*ξ)
            for i in 1:2
                τ = GetTime.(a,b,rand(2))
            end
            τ1, i0 = findmin(τ)
            t = t + τ1
            ξ = ξ + θ*τ1
            push!(Ξ, Skeleton(t, ξ))
            θ[i0] = -θ[i0]
        end
        return(Ξ)
    end
    function GetTime(a,b,u)
        if b > 0
            if a < 0
                τ = sqrt(-log(u)*2.0/b) - a/b
            else #a[i]>0
                τ = sqrt((a/b)^2 - log(u)*2.0/b) - a/b
            end
        elseif  b == 0
            if a > 0
                τ = -log(u)/a
            else #a[i] <= 0
                τ = Inf
            end
        else #b[i] < 0
            if a <= 0
                τ = Inf
            else
                if -log(u) <= -a^2/b + a^2/(2*b)
                    τ = - sqrt((a/b)^2 - log(u)*2.0/b) - a/b
                else
                    τ = Inf
                end
            end
        end
    end
    function Plots.plot(Ξ::Array{Skeleton,1},dims::Vector)
        d = length(Ξ)
        x = zeros(d)
        y = zeros(d)
        [x[i] = Ξ[i].ξ[dims[1]] for i in 1:d]
        [y[i] = Ξ[i].ξ[dims[2]] for i in 1:d]
        plot(x,y)
    end
end  # module ZigZag
