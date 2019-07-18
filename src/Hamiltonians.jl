module Hamiltonians
export pairing
using ..ManyBody
using ..ManyBody.Operators: @A

function pairing(δ, g)
    ret(X::DualBasis, Y) = ret(X', Y)
    function ret(X, Y)
        X, Y = supelem.((X, Y))
        SPB = spbasis(Y)

        f = sum(+SPB) do p
            δ*(p.level-1)*@A(p', p)(X, Y)
        end
        
        V = 0.0
        # Avoiding Iterators.product() here
        for p in +SPB, q in +SPB, r in +SPB, s in +SPB
            if p.level == q.level && r.level == s.level #=
            =# && spinup(p) && spindown(q) && spinup(r) && spindown(s)
            
                V += -g/2*@A(p', q', s, r)(X, Y)
            end
        end

        f + V
    end

    ret
end

end # Hamiltonians
