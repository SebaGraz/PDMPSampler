using Ciesielski
using LinearAlgebra
using Plots
using SparseArrays
include("Models/Gaussian.jl")
include("Faber1.jl")
include("PDSamplers.jl")

#L = levels of the Faber-Schauder basis
#a1,b1 parameters of the sde dX_t = (a + b X_t)dt + dW_t
#T timer of the ZigZagSampler before exiting
#ξ0 Initial points
function ZigZagExperiment( ξ0::Array{Float64, 1}, a1, b1, L, T::Float64)
    B = sparse(FillFingerMatrix(L, a1, b1))
    A = vector1(L, a1, b1)
    d = length(A) ; t = 0.0 ; θ = ones(d);
    ξ = ξ0
    Ξ = [Skeleton(t,ξ)]
    τ = zeros(d)
    c = []
    while t<T
        a = θ.*(-A - 2(B*ξ))
        b =-2θ.*(B*θ)
        τ = GetTime.(a,b,rand(d))
        τ1, i0 = findmin(τ)
        k, kk = Faber(i0)
        t = t + τ1
        ξ = ξ + θ*τ1
        push!(c, xy(k,kk))
        push!(Ξ, Skeleton(t, ξ))
        θ[i0] = -θ[i0]
    end
    return(Ξ, c)
end

#experiment
T = 21.0
L = 10
d = 2^(L+1) - 1 #dimensions given by the level L
a = 10.0
b = -1.0

#if you want to see the matrix
#V = FillFingerMatrix(L, a, b)

#ZigZag process having as posterior OU(a,b,1)
ZZ, c = ZigZagExperiment(zeros(d), a, b, L, T)

#Plot Path at time t
t = 21.0
plot([cies(ti, FindCoordinates(ZZ, t).ξ ,  L , 1.0) for ti in 0:0.001:1])


#Path updates while process running
plot()
for i in 1.0:5.0:21.0
    plot!([cies(ti, FindCoordinates(ZZ, i).ξ ,  L , 1.0) for ti in 0:0.001:1])
end
plot!()


#How the process mix
plot()
    for i in 1.0:0.5:21.0
        plot!([cies(ti, FindCoordinates(ZZ, i).ξ ,  L , 1.0) for ti in 0:0.001:1])
    end
plot!(leg=false)



# Analysis of the dimensions selected from the algorithm
x = zeros(length(c))
[x[i] = c[i].x[1] for i in 1:length(c)]
histogram(x)


# Average number of change of coordinates each basis
APB = []
for i in 1:L
    push!(APB, sum(x .== i )/2^(i -1))
end
scatter(APB)
