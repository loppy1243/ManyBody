struct MBState{R<:RefState, SP<:SPState}
    parts::IntSet
    holes::IntSet

    function MBState{R, SP}(phs) where {R, SP}
        ps = Set(pnum(p) for (p, _) in phs)
        hs = Set(pnum(h) for (_, h) in phs)

        new(ps, hs)
    end
end
MBState{R, SP}(phs::Vararg{Tuple{SP, SP}}) where {R, SP} = MBState{R, SP}(phs)
MBState{R, SP}(phs::Vararg{Tuple{Int, Int}} where {R, SP} = MBState{R, SP}(phs)
MBState(phs::Vararg{Tuple{SP, SP}}) where SP = MBState{Vaccuum{SP}, SP}(phs)
