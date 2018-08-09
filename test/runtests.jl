using Test
using ManyBody

for file in readdir(@__DIR__)
    file == basename(@__FILE__) && continue
    include(file)
end
