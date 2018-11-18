module Hamiltonians
export pairing
using ..ManyBody
using ..ManyBody.Operators: @A

pairing(δ, g) = function(X, Y)
    X = supelem(X)
    Y = supelem(Y)
    SPB = spbasis(Y)

    f = sum(SPB) do p
        δ*(p.level-1)*@A(p', p)(X, Y)
    end
    
    V = 0.0
    # Avoiding Iterators.product() here
    for p in SPB, q in SPB, r in SPB, s in SPB
        if p.level == q.level && r.level == s.level #=
        =# && spinup(p) && spindown(q) && spinup(r) && spindown(s)
        
            V += -g/2*@A(p', q', s, r)(X, Y)
        end
    end

    f + V
end

end # Hamiltonians
