# External Backends from Julia

Beyond the native Julia implementation in `MeijerG.jl`, there are a few other ways to evaluate Meijer G-functions from within Julia.
Meijer G-functions can be subtle, so having independent backends for sanity checks is often useful.

Of the available options, we primarily recommend the following two:
1. Calling Python [`mpmath`](https://mpmath.org/) via [`PyCall.jl`](https://github.com/JuliaPy/PyCall.jl) (free and open source).
2. Calling Wolfram Mathematica via [`MathLink.jl`](https://github.com/JuliaInterop/MathLink.jl) (requires a Wolfram Mathematica license).

We show how to get the mpmath option running since it's simpler, free and open source. See the [`MathLink.jl` docs](https://github.com/JuliaInterop/MathLink.jl) for hints on how to set up the Mathematica version.

## mpmath via PyCall

A helper file to set up mpmath's Meijer-G function to be callable in a similar style to the rest of this package (including both calling APIs) is provided as part of `MeijerG.jl` in [`extras/meijerg_mpmath.jl`](https://github.com/TSGut/MeijerG.jl/blob/main/extras/meijerg_mpmath.jl).
