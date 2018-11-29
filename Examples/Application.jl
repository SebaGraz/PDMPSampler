push!(LOAD_PATH, "C:\\Users\\sebas\\.julia\\dev\\ZigZag\\src")
using ZigZag
using Plots
using Statistics, Random, LinearAlgebra
ξ0 = zeros(2)
M = Normal(5.0, 1.0)
b = ZigZagGaussian(M, -10.0, 1000.0)
#Plot one dimensional Gaussian
plot(b[1],b[2])
inv(b[1][end]-b[1][1])*sum(diff(b[1]).* [(b[2][i]+b[2][i-1])/2 for i in 2:length(b[2])])

N= MultiNormal([1.0, 1.0, 1.0],[10.0 -9.0 9;
                            -9.0 10.0 -9.0;
                            9 -9.0 10.0])
b = ZigZagGaussian(N, zeros(3) , 1000.0)
plot(b,[1,3])
