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

macro commutes(expr::Expr)
    _commutes(:identity, expr)
end
macro commute(f, expr::Expr)
    _commutes(f, expr)
end
function _commutes(f, expr::Expr)
    # Just assume we have a function expr

    commed_expr = deepcopy(expr)

    if commed_expr.args[1].head === :where
        reverse!(@view commed_expr.args[1].args[1].args[2:end])
    else
        reverse!(@view commed_expr.args[1].args[2:end])
    end

    if f !== :identity
        commed_expr.args[2] = :($f(commed_expr.args[2]))
    end

    quote
        $(esc(expr))
        $(esc(commed_expr))
    end
end
