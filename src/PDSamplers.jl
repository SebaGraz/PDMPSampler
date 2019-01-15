struct Skeleton
    t::Float64
    ξ::Array{Float64, 1}
end
#ZigZag one dimensional unitary velocity
function ZigZagSampler(N::Normal, ξ0::Float64, T::Float64)
    t = 0.0 ; θ = 1; ξ = ξ0
    Ξ = [ξ]
    Δ = [t]
    while (t<T)
        τ = sqrt.(-2*N.σ*log(rand()) + max(0, θ*(ξ - N.μ))^2) - θ*(ξ - N.μ)
        t = t + τ
        ξ = ξ + θ*τ
        push!(Ξ, ξ)
        push!(Δ, t)
        θ = -θ
    end
    return(Δ, Ξ)
end

#ZigZag one dimensional with velocity
function ZigZagSampler(N::Normal, ξ0::Float64, v::Float64, T::Float64)
    t = 0.0 ; θ = v; ξ = ξ0
    Ξ = [ξ]
    Δ = [t]
    while (t<T)
        a= θ*ξ/N.σ
        b= θ^2/N.σ
        if a\b<0
            τ = sqrt(-2*log(rand())/b) - a/b
        else
            τ = sqrt(2/b*(log(rand()) + a^2/(2*b))) - a/b
        end
        t = t + τ
        ξ = ξ + θ*τ
        push!(Ξ, ξ)
        push!(Δ, t)
        θ = -θ
    end
    return(Δ, Ξ)
end

#ZigZag multidimensional
function ZigZagSampler(N::MultiNormal, ξ0::Array{Float64, 1} , T::Float64)
    d = length(N.μ) ; t = 0.0 ; θ = ones(d);
    ξ = ξ0
    Ξ = [Skeleton(t,ξ)]
    V = inv(N.Σ)
    while t<T
        b = θ.*(V*θ)
        a = θ.*(V*ξ)
        τ = GetTime.(a,b,rand(d))
        τ1, i0 = findmin(τ)
        t = t + τ1
        ξ = ξ + θ*τ1
        push!(Ξ, Skeleton(t, ξ))
        flip!(θ,i0)
    end
    return(Ξ)
end
#ZigZag n dimensional independent components
function ZigZagSampler(N::IDNormal, ξ0::Array{Float64, 1} , T::Float64)
    d = length(N.μ) ; t = 0.0 ; θ = ones(d);
    ξ = ξ0
    Ξ = [Skeleton(t,ξ)]
    V = inv(N.Σ)
    τ = zeros(d) #lazy initialization
    i0 = 1      #lazy inizialization
    while t<T
        b = θ[i0]* LinearAlgebra.dot(V[i0,:],θ)
        a = θ[i0]* LinearAlgebra.dot(V[i0,:],ξ)
        τ[i0] = GetTime(a, b, rand()) #ok
        τ1, i0 = findmin(τ) #ok
        τ[1:end .!= i0] .-= τ1 #ok
        t = t + τ1  #ok
        ξ = ξ + θ*τ1 #not necessary
        push!(Ξ, Skeleton(t, ξ)) #not necessary: create queue class as bouncy particle sampler
        flip!(θ,i0)
    end
    return(Ξ)
end

function BouncySampler(N::MultiNormal, ξ0::Array{Float64, 1} , T::Float64, λref)
    d = length(N.μ) ; t = 0.0 ; θ = ones(d);
    ξ = ξ0
    Ξ = [Skeleton(t,ξ)]
    V = inv(N.Σ) #precision matrix
    while t<T
        b = θ'*V*θ
        a = ξ'*V*θ
        τ = GetTime(a,b,rand())
        τref = GetTime(λref,0,rand())
        τ1 = min(τ, τref)
        t = t + τ1
        ξ = ξ + θ*τ1
        push!(Ξ, Skeleton(t, ξ))
        if (τ1 == τ)
            θ = reflection(V*ξ, θ)
        else
            θ = randn(d)
        end
    end
    return(Ξ)
end

function CoordinateSampler(N::MultiNormal, ξ0::Array{Float64, 1} , T::Float64)
    t = 0.0 ; θ = VelVect(1,1,N.d); #lazy initialization
    ξ = ξ0
    Ξ = [Skeleton(t,ξ)]
    while t<T
        b = θ'*N.V*θ #you are not using the fact that the velocity are zero eweryehre
        a = ξ'*N.V*θ
        τ = GetTime(a, b, rand()) #the same as Bouncy
        t = t + τ
        ξ = ξ + θ*τ
        push!(Ξ, Skeleton(t, ξ))
        θ = kernel(N, ξ, θ)
    end
    return(Ξ)
end

#obtaining waiting time for Inhomogeneous Poisson Process with rate of the form (a + b*t)^+
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
        elseif -log(u) <= -a^2/b + a^2/(2*b)
            τ = - sqrt((a/b)^2 - log(u)*2.0/b) - a/b
        else
            τ = Inf
        end
    end
end


#Plot Skeleton and choose the dimension
function Plots.plot(Ξ::Array{Skeleton,1},dims::Vector)
    x = [Ξ[i].ξ[dims[1]] for i in 1 : length(Ξ)]
    y = [Ξ[i].ξ[dims[2]] for i in 1 : length(Ξ)]
    plot(x,y)
end
#find coordinates for a fixed t
function FindCoordinates(Ξ::Array{Skeleton,1}, t::Float64)
    tt = [Ξ[i].t for i in 1:length(Ξ)]
    #you can throw the error if t is greater than the final point of the sequence.
    i = findfirst(x -> x>t ,tt)
    return(Skeleton(t, interpolate(Ξ[i-1].t, Ξ[i-1].ξ, Ξ[i].t, Ξ[i].ξ, t)))
end
