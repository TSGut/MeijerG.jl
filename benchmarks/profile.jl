using MeijerG
using Profile

include("common.jl")


"""
Run a short CPU profile of a representative heavy Meijer-G evaluation and
write a human-readable profile dump to `benchmarks/output/profile.txt`.

This avoids extra dependencies and produces a text profile suitable for
quick inspection. For interactive flamegraphs, use `ProfileView.jl` locally.
"""
function run_profile_benchmarks()
    println("Running profile benchmark...")

    # Choose a representative heavy configuration (BigFloat, moderate order)
    a = (big"0.25", big"0.75", big"1.25")
    b = (big"0.1", big"0.6", big"1.6")
    z = big"0.73"

    PROFILE_PATH = joinpath(OUTPUT_DIR, "profile.txt")

    # Clear any previous profile state
    Profile.clear()

    # Warm up and run a profiled evaluation at higher precision
    setprecision(256) do
        # quick warmup
        meijerg(a, b, 1, 1, z)

        # run under the profiler; use multiple evals to get meaningful samples
        Profile.@profile for _ in 1:6
            meijerg(a, b, 1, 1, z)
        end
    end

    # Write a textual profile summary
    open(PROFILE_PATH, "w") do io
        Profile.print(io)
    end

    println("Profile written to: $PROFILE_PATH")
    return nothing
end

# Execute the profile benchmark
run_profile_benchmarks()

# If you have ProfileView.jl installed, you can also view the profile interactively:
try 
    @eval import ProfileView
catch
    @warn "ProfileView.jl not available; skipping interactive profile view"
else
    println("Opening interactive profile view with ProfileView.jl...")
    ProfileView.view()
end