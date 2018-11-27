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

#Plot 2 coordinates of d-dimensional Gaussian
N= MultiNormal([0.0, 0.0],[10.0 -9.0; -9.0 10.0])
a = ZigZagGaussian(N, ξ0 , 10000.0)
d = length(a[2])
x = zeros(d)
y = zeros(d)
[x[i] = a[2][i][1] for i in 1:d]
[y[i] = a[2][i][2] for i in 1:d]
plot(x,y)
