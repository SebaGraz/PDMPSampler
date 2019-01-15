struct MultiNormal{T}
    μ::Vector{T}    #mean
    Σ::Matrix{T}    #Coavar Matrix
    V::Matrix{T}    #Precision Matrix Σ^-1
    d::Int64        #Dimension
    function MultiNormal(μ, Σ)
        new{Float64}(μ, Σ, inv(Σ), length(μ))
    end
end

struct Normal
    μ::Float64
    σ::Float64
    v::Float64
    function Normal(μ, σ)
        new(μ, σ, inv(σ))
    end
end

struct IDNormal{T}
    μ::Vector{T}    #mean
    Σ::Matrix{T}    #Coavar Matrix
    V::Matrix{T}    #Precision Matrix Σ^-1
    d::Int64        #Dimension
    function IDNormal(μ, Σ)
        new{Float64}(μ, Σ, inv(Σ), length(μ))
    end
end
