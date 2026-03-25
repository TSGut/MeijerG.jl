using Random: MersenneTwister, rand

const RANDOM_SAMPLE_COUNT = 5

function sample_real_points(fixed_points; seed::Int, lo::Float64, hi::Float64, count::Int=RANDOM_SAMPLE_COUNT)
    rng = MersenneTwister(seed)
    samples = Float64[fixed_points...]
    append!(samples, lo .+ (hi - lo) .* rand(rng, count))
    return samples
end

function sample_complex_points(
    fixed_points;
    seed::Int,
    real_lo::Float64,
    real_hi::Float64,
    imag_lo::Float64,
    imag_hi::Float64,
    count::Int=RANDOM_SAMPLE_COUNT,
)
    rng = MersenneTwister(seed)
    samples = ComplexF64[fixed_points...]
    append!(samples,
            [Complex(real_lo + (real_hi - real_lo) * rand(rng),
                    imag_lo + (imag_hi - imag_lo) * rand(rng)) for _ in 1:count])
    return samples
end
