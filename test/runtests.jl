using Test, Logging
using ManyBody

function get_ty_name(t)
    x = Meta.parse(string(t))
    if x isa Expr
        if x.head === :curly
            Expr(:curly,
                 x.args[1] isa Symbol ? x.args[1] : x.args[1].args[2].value,
                 x.args[2:end]...)
        else
            x.args[2].value
        end
    else
        x
    end |> string
end

with_logger(ConsoleLogger(stdout, Logging.Info)) do
    for file in readdir(@__DIR__)
        file == basename(@__FILE__) && continue
        include(file)
    end
end
