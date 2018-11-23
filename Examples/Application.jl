push!(LOAD_PATH, "C:\\Users\\sebas\\.julia\\dev\\ZigZag\\src")
using ZigZag
using Plots
ξ0 = zeros(2)
N= MultiNormal([5.0, 1.0],[10.0 0.5; 0.5 1.0])
M = Normal(5.0, 1.0)
a = ZigZagGaussian(N, ξ0 , 1000.0)
b = ZigZagGaussian(M, -10.0, 1000.0)
#Plot one dimensional Gaussian
plot(b[1],b[2])

#Plot 2 coordinates of d-dimensional Gaussian
d = length(a[2])
x = zeros(d)
y = zeros(d)
[x[i] = a[2][i][1] for i in 1:d]
[y[i] = a[2][i][2] for i in 1:d]
plot(x,y)
