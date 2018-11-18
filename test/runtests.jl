using Test, Logging
using ManyBody

### Why is this here...?
#function get_ty_name(t)
#    x = Meta.parse(string(t))
#    if x isa Expr
#        if x.head === :curly
#            Expr(:curly,
#                 x.args[1] isa Symbol ? x.args[1] : x.args[1].args[2].value,
#                 x.args[2:end]...)
#        else
#            x.args[2].value
#        end
#    else
#        x
#    end |> string
#end

for file in readdir(@__DIR__)
    file == basename(@__FILE__) && continue
    include(file)
end

with_logger(ConsoleLogger(stdout, Logging.Info)) do
    rlapplytest()
    normordtest()
    pairingtest(g_samples=5, atol=1e-5, mbbasis=Bases.Paired{4, 4})
end
