using ManyBody
using ManyBody.Operators: @A
using ManyBody.Hamiltonians: pairing

const LEVEL_SPACING = 1

_pairingH_direct(δ, g, i, j) = if i == j
    [2, 4, 6, 6, 8, 10][i]*δ - g
elseif i + j - 1 == 6
    0.0
else
    -g/2
end
pairingH_direct(δ, g) = [_pairingH_direct(δ, g, i, j) for i = 1:6, j = 1:6]

pairingtest(; g_samples, atol, mbbasis) = @testset "Pairing Hamiltonian" begin
    for g in (g_samples == 1 ? 0.5 : range(-1.0, stop=1.0, length=g_samples))
        @debug "Running pairing matrix test" g
        diff = tabulate(pairing(1, g), Float64, mbbasis) - pairingH_direct(1, g)
        @test all(abs.(diff) .< atol)
    end
end
