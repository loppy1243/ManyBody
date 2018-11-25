using ManyBody
using ManyBody.Operators: @A
using ManyBody.Hamiltonians: pairing

_pairing44H_direct(δ, g, i, j) = if i == j
    [2, 4, 6, 6, 8, 10][i]*δ - g
elseif i + j - 1 == 6
    0.0
else
    -g/2
end
pairing44H_direct(δ, g) = [_pairing44H_direct(δ, g, i, j) for i = 1:6, j = 1:6]

pairingtest(; g_samples, atol) = @testset "Pairing Hamiltonian" begin
    MBBASIS = Bases.Paired{4, 4}
    for g in (g_samples == 1 ? 0.5 : range(-1.0, stop=1.0, length=g_samples))
        @debug "Running pairing matrix test" g
        mat = tabulate(pairing(1, g), Array{Float64, 2}, MBBASIS)
        diff =  mat - pairing44H_direct(1, g)
        @test all(abs.(diff) .< atol)
    end
end
