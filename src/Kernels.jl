#Change of velocity at the poissont events
#ZigZag flips one coordinate
function flip!(v::Array{Float64,1}, i0::Int64)
    v[i0] = -v[i0]
    return(v)
end

#Bouncy Particle changes all the coordinates
function reflection(g,v)
    v .-= 2.0*dot(g,v)*g./dot(g,g)
end

#CoordinateSampler change its coordinate
#but the change its stochastic

function kernel(N::MultiNormal,x,θ)
    p = zeros(N.d)  #inizialization probability vector
    s = fill(1 ,N.d)   #inizialization signs
    precmu = N.V*x  #precomputation precision matrix * position
    num = -precmu  #precomputation
    den = sum(abs.(precmu))
    for i = 1:N.d
        if num[i]<= 0
            s[i] = -1
            p[i] = -num[i]/den
        else
            p[i] = num[i]/den
        end
    end
    p = cumsum(p)
    println(p) #DEBUG should be monotone and ending around 1
    i0 = findfirst(x -> x>rand(), p)
    return (VelVect(s[i0],i0,N.d))
end

#Help functions
function VelVect(s,i,d)
    v = zeros(d)
    v[i] = s*1
    return v
end
