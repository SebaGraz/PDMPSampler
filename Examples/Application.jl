#push!(LOAD_PATH, "C:\\Users\\sebas\\.julia\\dev\\ZigZag\\src")
using ZigZag
using Plots
ξ0 = zeros(2)
N= MultiNormal([1.0, 1.0],[1.0 0.5; 0.5 1.0])
M = Normal(1.0, 1.0)
a = ZigZagGaussian(N, ξ0 , 100.0)
b = ZigZagGaussian(M, -10.0, 100.0)

plot(b[1],b[2])
