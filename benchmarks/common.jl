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

function maybe_plot(path::AbstractString, x, y; xlabel::String, ylabel::String, title::String,
                    yscale::Symbol = :identity, marker = :circle)
    try
        @eval import Plots
    catch
        @warn "Plots.jl not available; skipping plot output at $path"
        return false
    end

    plt = Plots.plot(x, y;
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        marker = marker,
        linewidth = 2,
        legend = false,
        yscale = yscale,
    )
    Plots.savefig(plt, path)
    return true
end
