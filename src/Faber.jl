function Faber(i,j)
    return 2^i + j
end
function Faber(n)
    return (floor(log2(n)), n - 2^(floor(log2(n))) )
end
%check
Faber(1,1)
Faber(1)

function FillFingerMatrix(L::Int64)
    dim = 2^(L+1) -1
    V = zeros(dim,dim)
    for i in 1:dim
        for j in 1:i
            if i == j
                i1, j1 = Faber(i)
                V[i,i] = 1
            else
                i1, j1 = Faber(j)
                i2, j2 = Faber(i)
                if j1*2^(i2 - i1) <= j2 < (j1 + 1)*2^(i2 - i1)     #i1 = i2 è già considerato quà? penso di si
                    #se non è 0, calcola
                    V[i,j], V[j,i] = (1,1) #qualcosa
                else
                    V[i,j], V[j,i] = (0 , 0 )
                end
            end
        end
    end
    return V
end


#Check
L = 2
V = FillFingerMatrix(L)
#check the symmetry and length
V' == V
2^(L+1)-1
