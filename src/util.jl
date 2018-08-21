struct MethodDisallowedError <: Exception
    me::MethodError
end
MethodDisallowedError(f, args) = MethodDisallowedError(MethodError(f, args))
function Base.showerror(io::IO, ex::MethodDisallowedError)
    println(io, "Method explicitly disallowed")
    showerror(io, ex.me)
    nothing
end
macro disallow(expr::Expr)
    ## Make sure this is a function expr
    @assert expr.head === :call || expr.head === :where #=
                             =# && expr.args[1].head === :call

    _get_sym!(x::Symbol) = x
    _get_sym!(x::Expr) = if x.head == :(::) && length(x.args) == 1
        ret = gensym()
        x.args = vcat(ret, x.args)
        ret
    else
        _get_sym!(x.args[1])
    end

    func_sym, arg_syms = if expr.head === :where
        (_get_sym!(expr.args[1].args[1]), map(_get_sym!, expr.args[1].args[2:end]))
    else
        (_get_sym!(expr.args[1]), map(_get_sym!, expr.args[2:end]))
    end

    func_sym_esc = esc(func_sym)
    arg_sym_escs = map(esc, arg_syms)

    :($(esc(expr)) = throw(MethodDisallowedError($func_sym_esc, ($(arg_sym_escs...),))))
end

### Better name: @symmetric
macro commutes(expr::Expr)
    _commutes(:(:), :identity, expr)
end
macro commutes(ns, expr::Expr)
    _commutes(ns, :identity, expr)
end
macro commutes(ns, f, expr::Expr)
    _commutes(ns, f, expr)
end
## Add more input validation
function _commutes(ns, f, expr::Expr)
    # Just assume we have a function exp

    _get_nargs(expr) = if expr.args[1].head === :where
        length(expr.args[1].args[1].args) - 1
    else
        length(expr.args[1].args) - 1
    end

    ixs = if ns === :(:)
        2:_get_nargs(expr)+1
    elseif ns isa Expr && ns.head === :tuple
        if all(x -> x isa Int, ns.args)
            map(x -> x+1, ns.args)
        else
            error("Expected literal NTuple{<:Any, Int}, found $ns")
        end
    elseif ns isa Expr && ns.head === :call && ns.args[1] === :(:)
        if length(ns.args[2:end]) == 2 && all(x -> x isa Int, ns.args[2:end])
            (ns.args[2]+1):(ns.args[3]+1)
        else
            error("Expected UnitRange{Int}, found $ns")
        end
    else
        error("Expected literal range or tuple, found $ns")
    end

    commed_expr = deepcopy(expr)

    if commed_expr.args[1].head === :where
        reverse!(@view commed_expr.args[1].args[1].args[ixs])
    else
        reverse!(@view commed_expr.args[1].args[ixs])
    end

    if f !== :identity
        commed_expr.args[2] = :($f($(commed_expr.args[2])))
    end

    quote
        $(esc(expr))
        $(esc(commed_expr))
    end
end
