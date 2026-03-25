using MeijerG
using Statistics, Random
using Plots

include("common.jl")

"""
Run threading benchmarks by exercising `meijerg` at multiple orders and
testing all thread counts from 1 up to the system maximum. The script uses a
controller/worker pattern: the controller spawns Julia worker processes with
different `-t` thread counts and the workers print `RESULT:<threads>:<order>:<time>`
lines which the controller parses.

The benchmark accepts the following optional environment variables:
- `BENCHMARK_TYPE` = "Float64" (default) or "BigFloat" to run BigFloat workloads.
- `BENCHMARK_PREC` = precision (bits) used when `BENCHMARK_TYPE=BigFloat` (default: 256).

Output: `benchmarks/output/threading.csv` and
`benchmarks/output/threading.png`.
"""
function run_threading_benchmarks()
    # Worker mode: run a benchmark for a single thread count and print results
    if haskey(ENV, "BENCHMARK_THREAD_COUNT")
        n_threads = parse(Int, ENV["BENCHMARK_THREAD_COUNT"])

        # Test type selection
        bench_type = get(ENV, "BENCHMARK_TYPE", "Float64")
        is_big = lowercase(bench_type) == "bigfloat"
        bench_prec = tryparse(Int, get(ENV, "BENCHMARK_PREC", "256")) === nothing ? 256 : parse(Int, get(ENV, "BENCHMARK_PREC", "256"))

        # Use larger orders and more repeats to reduce noise in timing
        orders = [6, 8, 10, 12, 14]

        results = Dict{Int,Float64}()
        for p in orders
            if is_big
                a = Tuple(big(0.11 + 0.23 * i) for i in 1:p)
                b = Tuple(big(0.07 + 0.31 * i) for i in 1:p)
                z = big(0.65)
                # run timing under a fixed BigFloat precision
                t = setprecision(bench_prec) do
                    timed_median(() -> meijerg(a, b, 1, 1, z); warmup = 3, repeats = 11)
                end
            else
                a = Tuple(0.11 + 0.23 * i for i in 1:p)
                b = Tuple(0.07 + 0.31 * i for i in 1:p)
                z = 0.65
                t = timed_median(() -> meijerg(a, b, 1, 1, z); warmup = 3, repeats = 11)
            end
            results[p] = t
        end

        # Emit parseable lines for the controller to collect (include bench type)
        for p in orders
            println("RESULT:$n_threads:$bench_type:$p:$(results[p])")
        end
        return nothing
    end

    # Controller mode: spawn a Julia process per thread count and collect results
    max_threads_system = Sys.CPU_THREADS
    max_threads_bench = tryparse(Int, get(ENV, "BENCHMARK_MAX_THREADS", ""))
    max_threads = max_threads_bench !== nothing ? max_threads_bench : max_threads_system

    println("\nTHREADING BENCHMARK — testing 1..$max_threads threads (system: $max_threads_system)\n")

    # all_results[threads][bench_type][order] = time
    all_results = Dict{Int,Dict{String,Dict{Int,Float64}}}()
    bench_types = ["Float64", "BigFloat"]
    project_root = normpath(joinpath(@__DIR__, ".."))
    script_path = joinpath(project_root, "benchmarks", "threading.jl")
    julia_exe = Base.julia_cmd().exec[1]

    for n_threads in 1:max_threads
        print("Benchmarking with $n_threads thread(s)... ")
        flush(stdout)

        # Run both numeric types for each thread count
        all_results[n_threads] = Dict{String,Dict{Int,Float64}}()
        for bench_type in bench_types
            # Propagate BENCHMARK_PREC from controller env if present, but force BENCHMARK_TYPE
            pre_assign = "ENV[\"BENCHMARK_TYPE\"]=\"" * bench_type * "\"; "
            if get(ENV, "BENCHMARK_PREC", "") != ""
                pre_assign *= "ENV[\"BENCHMARK_PREC\"]=\"" * get(ENV, "BENCHMARK_PREC", "") * "\"; "
            end

            cmd = Cmd([julia_exe, "--project=$project_root", "-t", "$n_threads", "-e", string(pre_assign, "ENV[\"BENCHMARK_THREAD_COUNT\"]=\"$n_threads\"; include(\"$script_path\")")])
            output = read(cmd, String)

            all_results[n_threads][bench_type] = Dict{Int,Float64}()
            for line in split(output, '\n')
                if startswith(line, "RESULT:")
                    parts = split(line, ":")
                    if length(parts) == 5
                        n_t = parse(Int, parts[2])
                        type_str = parts[3]
                        p = parse(Int, parts[4])
                        time_val = parse(Float64, parts[5])
                        # initialize container if necessary
                        if !haskey(all_results, n_t)
                            all_results[n_t] = Dict{String,Dict{Int,Float64}}()
                        end
                        if !haskey(all_results[n_t], type_str)
                            all_results[n_t][type_str] = Dict{Int,Float64}()
                        end
                        all_results[n_t][type_str][p] = time_val
                    end
                end
            end
        end

        println("✓")
    end

    # Flatten to rows and write CSV (include numeric type)
    rows = Vector{Tuple{Int,String,Int,Float64}}()
    for t in sort(collect(keys(all_results)))
        for type_str in sort(collect(keys(all_results[t])))
            for p in sort(collect(keys(all_results[t][type_str])))
                push!(rows, (t, type_str, p, all_results[t][type_str][p]))
            end
        end
    end
    csv_path = joinpath(OUTPUT_DIR, "threading.csv")
    write_table(csv_path, ["threads", "type", "order_pq", "median_seconds"], rows)
    println("Threading CSV written to: $csv_path")

    # Produce one plot per numeric type
    for type_str in bench_types
        # find orders from the first thread that has this type
        sample_thread = first(sort(collect(keys(all_results))))
        orders = sort(collect(keys(all_results[sample_thread][type_str])))
        plt = Plots.plot(xlabel = "Order (p=q)", ylabel = "Median runtime (s)", title = "Threading benchmark ($type_str)", legend = :topleft)
        for t in sort(collect(keys(all_results)))
            # if this thread didn't record the type, skip
            if !haskey(all_results[t], type_str)
                continue
            end
            times = [all_results[t][type_str][o] for o in orders]
            Plots.plot!(plt, orders, times; marker = :circle, label = "$(t) threads")
        end
        savepath = joinpath(OUTPUT_DIR, "threading_$(type_str).png")
        Plots.savefig(plt, savepath)
        println("Threading plot saved to: $savepath")
    end

    return nothing
end

run_threading_benchmarks()