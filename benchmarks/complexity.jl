using MeijerG
using Plots

include("common.jl")

"""
Run time-complexity style measurements:
1) scaling with order p=q for a fixed (m,n) and Float64 arithmetic
2) scaling with BigFloat precision for a fixed parameter set
"""
function run_complexity_benchmarks()
    println("Running complexity benchmarks...")

    orders = collect(2:12)
    z = 0.65
    rows_order = Tuple{Int, Float64}[]

    for p in orders
        a = Tuple(0.11 + 0.23 * i for i in 1:p)
        b = Tuple(0.07 + 0.31 * i for i in 1:p)
        t = timed_median(() -> meijerg(a, b, 1, 1, z); repeats = 11)
        push!(rows_order, (p, t))
        println("  order p=q=$p: median time = $(round(t * 1e3, digits=3)) ms")
    end

    order_csv = joinpath(OUTPUT_DIR, "complexity_order.csv")
    write_table(order_csv, ["order_pq", "median_seconds"], rows_order)

    order_plot = Plots.plot(
        first.(rows_order),
        last.(rows_order);
        xlabel = "Order (p=q)",
        ylabel = "Median runtime (s)",
        title = "MeijerG runtime vs order",
        marker = :circle,
        linewidth = 2,
        legend = false,
    )
    Plots.savefig(order_plot, joinpath(OUTPUT_DIR, "complexity_order.png"))

    precisions = [64, 96, 128, 192, 256, 384, 512]
    rows_precision = Tuple{Int, Float64}[]
    a_big = (big"0.25", big"0.75", big"1.35")
    b_big = (big"0.1", big"0.6", big"1.6")
    z_big = big"0.8"

    for prec in precisions
        t = setprecision(prec) do
            timed_median(() -> meijerg(a_big, b_big, 1, 1, z_big); repeats = 9)
        end
        push!(rows_precision, (prec, t))
        println("  precision=$prec bits: median time = $(round(t * 1e3, digits=3)) ms")
    end

    precision_csv = joinpath(OUTPUT_DIR, "complexity_precision.csv")
    write_table(precision_csv, ["precision_bits", "median_seconds"], rows_precision)

    precision_plot = Plots.plot(
        first.(rows_precision),
        last.(rows_precision);
        xlabel = "BigFloat precision (bits)",
        ylabel = "Median runtime (s)",
        title = "MeijerG runtime vs precision",
        marker = :circle,
        linewidth = 2,
        legend = false,
    )
    Plots.savefig(precision_plot, joinpath(OUTPUT_DIR, "complexity_precision.png"))

    println("Complexity benchmark outputs written to: $OUTPUT_DIR")
    return nothing
end

run_complexity_benchmarks()

