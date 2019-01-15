#the whole script is used for constructing the Sparse Matrix and the vector
#used for simulating the Bridge
function Faber(i, j)
    return 2^i + j
end
function Faber(n)
    return (floor(log2(n)), n - 2^(floor(log2(n))))
end

#integral of quadratic function1
#int^n1_n2 a*t^2 + b*t dt
function IntQuad1(a, b, n1, n2)
    a*(n1^3 - n2^3)/3 + b*(n1^2 - n2^2)/2
end

#ordering from the larger basis to the smaller
function OrderBasis(i, j)
    if i <= j
        return (i, j)
    else
        return (j, i)
    end
end

function OrderBasis(i1 ,j1 ,i2 ,j2)
    if i1 <= i2
        return (i1, j1, i2, j2)
    else
        return (i2, j2, i1, j1)
    end
end

#check when the funcitons share their support
function ShareDomain(i,j)
    if i == j return 1
    else
        (i, j) = OrderBasis(i, j)
        i1, j1 = Faber(i)
        i2, j2 = Faber(j)
        if j1*2^(i2 - i1) <= j2 < (j1 + 1)*2^(i2 - i1)
            return 1
        else
            return 0
        end
    end
end
function ShareDomain(i1,j1,i2,j2)
    if i1 == i2 && j1 == j2  return 1
    else
        (i1, j1, i2, j2) = OrderBasis(i1, j1, i2, j2)

        if j1*2^(i2 - i1) <= j2 < (j1 + 1)*2^(i2 - i1)
            return 1
        else
            return 0
        end
    end
end

#flip i,j with respect to psi(0,0). It makes the integration easy
function flip(i,j)
    bound = ((1/2)*2^i - 1/2)
    if bound < j
        j = 2*bound - j
    end
    return Integer(j)
end

#int_0^1 psi_i (t) psi_j(t) dt with i<j
function DproductSchauder(i,j)
    (i, j) = OrderBasis(i, j)
    i1, j1 = Faber(i) #ordering  i < j --> i1 < i2
    i2, j2 = Faber(j)
    i3 = i2 - i1
    j3 = j2 - 2^(i2-i1)*j1
    j3 = flip(i3,j3)
    A = j3/2^i3
    B = (j3 + 1/2)/2^i3
    C = (j3 + 1)/ 2^i3
    return 2^((i2 - 5*i1)/2)*(IntQuad1(1,-j3/2^i3, B, A) - IntQuad1(1,-(j3+1)/2^i3, C, B))
end

#int_0^1 psi_i (t) d(psi_j(t)) with i<j
function secondterm(i,j)
    (i, j) = OrderBasis(i, j)
    i1, j1 = Faber(i) #ordering  i1<i2 --> i < j
    i2, j2 = Faber(j)
    i3 = i2 - i1
    j3 = j2 - 2^(i2-i1)*j1
    j4 = flip(i3,j3)
    A = j4/2^i3
    B = (j4 + 1/2)/2^i3
    C = (j4 + 1)/ 2^i3
    if j4 == j3
        return 2^((i2 - 5*i1)/2)*(IntQuad1(0, 1, B, A) - IntQuad1(0, 1, C, B))
    else
        return -2^((i2 - 5*i1)/2)*(IntQuad1(0, 1, B, A) - IntQuad1(0, 1, C, B))
    end
end

#construction of the matrix as pdf
function FillFingerMatrix(L::Int64)
    dim = 2^(L+1) -1
    V = zeros(dim,dim)
    for i in 1:dim
        for j in 1:i
            if i == j
                i1, j1 = Faber(i)
                V[i,i] = 2^(-2*i1 - 2)/3
            else
                if ShareDomain(i,j)==1
                    res = DproductSchauder(i, j) - secondterm(i, j)
                    V[i,j], V[j,i] = (res, res)
                else
                    V[i,j], V[j,i] = (0 , 0)
                end
            end
        end
    end
    return V
end

#coordinate Faber Schauder basis
struct xy
    x::Int64
    y::Int64
end

#process dX_t = (a + bX_t)dt + dW_t

#finger matrix case specific for Ornstein Uhlembeck model
#See pdf for more information

# int_0^1 psi_i (t) psi_j (t) dt for i,t different
function DproductSchauder(i,j)
    (i, j) = OrderBasis(i, j)
    i1, j1 = Faber(i) #ordering  i < j --> i1 < i2
    i2, j2 = Faber(j)
    i3 = i2 - i1
    j3 = j2 - 2^(i2-i1)*j1
    j3 = flip(i3,j3)

    A = j3/2^i3
    B = (j3 + 1/2)/2^i3
    C = (j3 + 1)/ 2^i3
    return 2^((i2 - 5*i1)/2)*(IntQuad1(1,-j3/2^i3, B, A) - IntQuad1(1,-(j3+1)/2^i3, C, B))
end

#int_0^1 psi_i (t) d(psi_j(t)) with i<j
# when i = j  secondterm(i,i) != 0 --> you cant use this function
function secondterm(i,j)
    (i, j) = OrderBasis(i, j)
    i1, j1 = Faber(i) #ordering  i1<i2 --> i < j
    i2, j2 = Faber(j)
    i3 = i2 - i1
    j3 = j2 - 2^(i2-i1)*j1
    j4 = flip(i3,j3)
    A = j4/2^i3
    B = (j4 + 1/2)/2^i3
    C = (j4 + 1)/ 2^i3
    if j4 == j3
        return 2^((i2 - 5*i1)/2)*(IntQuad1(0, 1, B, A) - IntQuad1(0, 1, C, B))
    else
        return -2^((i2 - 5*i1)/2)*(IntQuad1(0, 1, B, A) - IntQuad1(0, 1, C, B))
    end
end

#construction of the matrix
function FillFingerMatrix(L::Int64,a,b)
    dim = 2^(L+1) -1
    V = zeros(dim,dim)
    for i in 1:dim
        for j in 1:i
            if i == j
                i1, j1 = Faber(i)
                V[i,i] = -b^2/2*2^(-2*i1 - 2)/3  - 1/2
            else
                if ShareDomain(i,j)==1
                    res = b*secondterm(i, j) - b^2/2*DproductSchauder(i, j)
                    V[i,j], V[j,i] = (res, res)
                else
                    V[i,j], V[j,i] = (0 , 0) #puo1 essere evitato
                end
            end
        end
    end
    return V
end

#consructing vecor
#int_0^1 psi_i(t) dt (not j dependent)
function integrateSchauder(i)
    return 2^((-3*i - 4 )/2)
end
function vector1(L, a, b)
    A = zeros(2^(L+1) - 1)
    for i in 1:length(A)
        A[i] = -a*b*integrateSchauder(floor(log2(i)))
    end
    return A
end
