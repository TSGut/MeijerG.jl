using MeijerG
using Plots

include("common.jl")

"""
Run convergence-focused measurements:
1) error vs BigFloat precision using a 1024-bit reference value
2) stability scan near the |z|=1 branch-selection boundary (p=q)
"""
function run_convergence_benchmarks()
    println("Running convergence benchmarks...")

    setups = [
        (name = "setup1", a = (0.25, 0.75, 1.25), b = (0.1, 0.6, 1.6), m = 1, n = 1, z_precision = 0.73),
        (name = "setup2", a = (0.30, 1.10), b = (0.2, 1.4), m = 1, n = 1, z_precision = 0.67),
        (name = "setup3", a = (0.45, 0.95, 1.55, 2.05), b = (0.12, 0.72, 1.22, 1.82), m = 1, n = 1, z_precision = 0.81),
        (name = "setup4", a = (0.18, 0.88, 1.48), b = (0.05, 0.55, 1.05), m = 1, n = 1, z_precision = 0.59),
        (name = "setup5", a = (0.40, 1.20, 2.00), b = (0.15, 0.85, 1.65), m = 1, n = 1, z_precision = 0.92),
    ]

    precisions = [64, 96, 128, 192, 256, 384, 512, 768]
    rows_precision = Tuple{String, Int, Float64, Float64}[]

    for setup in setups
        a_big = map(BigFloat, setup.a)
        b_big = map(BigFloat, setup.b)
        z_big = BigFloat(setup.z_precision)

        reference = setprecision(1024) do
            meijerg(a_big, b_big, setup.m, setup.n, z_big)
        end

        println("  $(setup.name):")
        for prec in precisions
            abs_err, rel_err = setprecision(prec) do
                val = meijerg(a_big, b_big, setup.m, setup.n, z_big)
                ae = abs(val - reference)
                re = ae / max(abs(reference), eps(BigFloat))
                ae, re
            end
            abs64 = Float64(abs_err)
            rel64 = Float64(rel_err)
            push!(rows_precision, (setup.name, prec, rel64, abs64))
            println("    precision=$prec bits: relative error = $rel64, absolute error = $abs64")
        end
    end

    precision_csv = joinpath(OUTPUT_DIR, "convergence_precision.csv")
    write_table(precision_csv, ["setup", "precision_bits", "relative_error", "absolute_error"], rows_precision)

    precision_rel_plot = Plots.plot(
        xlabel = "BigFloat precision (bits)",
        ylabel = "Relative error vs 1024-bit reference",
        title = "MeijerG convergence with precision",
        legend = :topright,
        yscale = :log10,
    )
    precision_abs_plot = Plots.plot(
        xlabel = "BigFloat precision (bits)",
        ylabel = "Absolute error vs 1024-bit reference",
        title = "MeijerG convergence with precision (absolute)",
        legend = :topright,
        yscale = :log10,
    )
    for setup in setups
        setup_rows = filter(row -> row[1] == setup.name, rows_precision)
        bits = [row[2] for row in setup_rows]
        Plots.plot!(precision_rel_plot, bits, [row[3] for row in setup_rows];
            marker = :circle, linewidth = 2, label = setup.name)
        Plots.plot!(precision_abs_plot, bits, [row[4] for row in setup_rows];
            marker = :circle, linewidth = 2, label = setup.name)
    end
    Plots.savefig(precision_rel_plot, joinpath(OUTPUT_DIR, "convergence_precision.png"))
    Plots.savefig(precision_abs_plot, joinpath(OUTPUT_DIR, "convergence_precision_absolute.png"))

    z_values = [
        0.90,
        0.94,
        0.96,
        0.98,
        0.99,
        0.995,
        0.9975,
        0.999,
        0.9995,
        0.9999,
        1.0001,
        1.0005,
        1.001,
        1.0025,
        1.005,
        1.01,
        1.02,
        1.04,
        1.06,
        1.08,
        1.10,
    ]
    rows_boundary = Tuple{String, Float64, Float64, Float64}[]

    for setup in setups
        a64 = setup.a
        b64 = setup.b
        println("  $(setup.name): stability around |z|=1")
        for z0 in z_values
            val64 = meijerg(a64, b64, setup.m, setup.n, z0)
            val_ref = setprecision(1024) do
                meijerg(map(BigFloat, a64), map(BigFloat, b64), setup.m, setup.n, BigFloat(z0))
            end
            abs_err = abs(BigFloat(val64) - val_ref)
            rel_err = abs_err / max(abs(val_ref), eps(BigFloat))
            abs64 = Float64(abs_err)
            rel64 = Float64(rel_err)
            push!(rows_boundary, (setup.name, z0, rel64, abs64))
            println("    z=$z0: relative error = $rel64, absolute error = $abs64")
        end
    end

    boundary_csv = joinpath(OUTPUT_DIR, "convergence_boundary.csv")
    write_table(boundary_csv, ["setup", "z_value", "relative_error_vs_1024bit", "absolute_error_vs_1024bit"], rows_boundary)

    boundary_rel_plot = Plots.plot(
        xlabel = "z (real axis, p=q branch switch near |z|=1)",
        ylabel = "Relative error vs 1024-bit reference",
        title = "MeijerG stability near |z|=1",
        legend = :topright,
        yscale = :log10,
    )
    boundary_abs_plot = Plots.plot(
        xlabel = "z (real axis, p=q branch switch near |z|=1)",
        ylabel = "Absolute error vs 1024-bit reference",
        title = "MeijerG stability near |z|=1 (absolute)",
        legend = :topright,
        yscale = :log10,
    )
    for setup in setups
        setup_rows = filter(row -> row[1] == setup.name, rows_boundary)
        zs = [row[2] for row in setup_rows]
        Plots.plot!(boundary_rel_plot, zs, [row[3] for row in setup_rows];
            marker = :circle, linewidth = 2, label = setup.name)
        Plots.plot!(boundary_abs_plot, zs, [row[4] for row in setup_rows];
            marker = :circle, linewidth = 2, label = setup.name)
    end
    Plots.savefig(boundary_rel_plot, joinpath(OUTPUT_DIR, "convergence_boundary.png"))
    Plots.savefig(boundary_abs_plot, joinpath(OUTPUT_DIR, "convergence_boundary_absolute.png"))

    println("Convergence benchmark outputs written to: $OUTPUT_DIR")
    return nothing
end

run_convergence_benchmarks()
