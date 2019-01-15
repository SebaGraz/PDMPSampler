using Plots
using LinearAlgebra
include("Models/Gaussian.jl")
include("PDSamplers.jl")
include("Kernels.jl")

# Comparison of the three algorithms
μ = zeros(2)
Σ = [1.0 -0.5 ; -0.5 1.0 ]
N = MultiNormal(μ,Σ)

a = BouncySampler(N, zeros(2) , 1000.0, 1.0)
plot(a,[1,2])

b = CoordinateSampler(N, zeros(2) , 1000.0)
plot(b,[1,2])

c = ZigZagSampler(N, zeros(2) , 10000.0)
plot(c,[1,2])

# Faster algorithm for Independent ZigZag
M = IDNormal(zeros(2),[1.0 0.0 ; 0.0 1.0])
d = ZigZagSampler(M, zeros(2) , 10000.0)
plot(d ,[1,2])

# one variable ZigZag varying the velocity
Z = Normal(0.0, 1.0)
x = ZigZagSampler(Z, -10.0, 10.0, 100.0)
y = ZigZagSampler(Z, -10.0, 1000.0)
plot(x[1],x[2])
plot(y[1],y[2])
