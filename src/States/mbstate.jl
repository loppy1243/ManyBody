struct MBState{R<:RefState, SP<:SPState}
    parts::BitVector
    holes::BitVector

    function MBState{R, SP}(phs) where {R, SP}
        ps = falses(nholes(R))
        hs = falses(nparts(R))

        for (p, h) in phs
            pn = pnum(R, p)
            hn = hnum(R, h)
            @assert pn != hn
            @assert !ps[pn] && !hs[hn]

            ps[pn] = true
            hs[hn] = true
        end

        new(ps, hs)
    end
end
MBState{R, SP}(phs::Vararg{Tuple{SP, SP}}) where {R, SP} = MBState{R, SP}(phs)
MBState{R, SP}(phs::Vararg{Tuple{Int, Int}} where {R, SP} = MBState{R, SP}(phs)
MBState(phs::Vararg{Tuple{SP, SP}}) where SP = MBState{Vaccuum{SP}, SP}(phs)

Base.:(==)(::MBState, ::MBState) = false
Base.:(==)(s1::MBState{R, S}, s2::MBState{R, S}) =
    s1.parts == s2.parts && s1.holes == s2.holes

function snum(s::MBState)
    bs = [s.holes; s.parts]

    sum(bs[i]*2^i for i = 0:length(bs)-1)
end

nstates(::Type{MBState{RS, SP}}) where {RS, SP} = 2^nstates(SP) - 1
