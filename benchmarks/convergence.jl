using MeijerG

include("common.jl")

"""
Run convergence-focused measurements:
1) error vs BigFloat precision using a 1024-bit reference value
2) stability scan near the |z|=1 branch-selection boundary (p=q)
"""
function run_convergence_benchmarks()
    println("Running convergence benchmarks...")

    a = (big"0.25", big"0.75", big"1.25")
    b = (big"0.1", big"0.6", big"1.6")
    z = big"0.73"

    reference = setprecision(1024) do
        meijerg(a, b, 1, 1, z)
    end

    precisions = [64, 96, 128, 192, 256, 384, 512, 768]
    rows_precision = Tuple{Int, BigFloat}[]

    for prec in precisions
        err = setprecision(prec) do
            val = meijerg(a, b, 1, 1, z)
            abs(val - reference) / max(abs(reference), eps(BigFloat))
        end
        push!(rows_precision, (prec, err))
        println("  precision=$prec bits: relative error = $err")
    end

    precision_csv = joinpath(OUTPUT_DIR, "convergence_precision.csv")
    write_table(precision_csv, ["precision_bits", "relative_error"], rows_precision)

    maybe_plot(
        joinpath(OUTPUT_DIR, "convergence_precision.png"),
        first.(rows_precision),
        Float64.(last.(rows_precision)),
        xlabel = "BigFloat precision (bits)",
        ylabel = "Relative error vs 1024-bit reference",
        title = "MeijerG convergence with precision",
        yscale = :log10,
    )

    a64 = (0.25, 0.75)
    b64 = (0.1, 0.6)
    z_values = [0.90, 0.99, 0.999, 1.001, 1.01, 1.10]
    rows_boundary = Tuple{Float64, Float64}[]

    for z0 in z_values
        val64 = meijerg(a64, b64, 1, 1, z0)
        val_ref = setprecision(1024) do
            meijerg((big"0.25", big"0.75"), (big"0.1", big"0.6"), 1, 1, BigFloat(z0))
        end
        err = abs(BigFloat(val64) - val_ref) / max(abs(val_ref), eps(BigFloat))
        push!(rows_boundary, (z0, Float64(err)))
    end

    boundary_csv = joinpath(OUTPUT_DIR, "convergence_boundary.csv")
    write_table(boundary_csv, ["z_value", "relative_error_vs_1024bit"], rows_boundary)

    maybe_plot(
        joinpath(OUTPUT_DIR, "convergence_boundary.png"),
        first.(rows_boundary),
        last.(rows_boundary),
        xlabel = "z (real axis, p=q branch switch near |z|=1)",
        ylabel = "Relative error vs 1024-bit reference",
        title = "MeijerG stability near |z|=1",
        yscale = :log10,
    )

    println("Convergence benchmark outputs written to: $OUTPUT_DIR")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_convergence_benchmarks()
end
