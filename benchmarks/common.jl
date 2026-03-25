using DelimitedFiles
using Statistics

const OUTPUT_DIR = joinpath(@__DIR__, "output")
mkpath(OUTPUT_DIR)

function timed_median(f; warmup::Int = 2, repeats::Int = 9)
    for _ in 1:warmup
        f()
    end
    times = [@elapsed f() for _ in 1:repeats]
    return median(times)
end

function write_table(path::AbstractString, header::Vector{String}, rows::Vector{<:Tuple})
    open(path, "w") do io
        println(io, join(header, ","))
        for row in rows
            println(io, join(row, ","))
        end
    end
end

