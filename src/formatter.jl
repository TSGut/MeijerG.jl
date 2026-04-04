"""
    @formatter expr
    @formatter expr verbose=false

Apply the `@formatter` macro to a `meijerg` function call to automatically format its output.

# Arguments
- `expr`: A `meijerg` function call with its arguments.
- `verbose`: (Optional, default `true`) Pass `false` to suppress printed output.

# Returns
A tuple `(plaintext_str, latex_str, mathematica_str)` containing the Plaintext, LaTeX, and Mathematica representations of the Meijer G-function.

# Examples
```julia
@formatter meijerg((1, 2), (3, 4), 1, 1, 5)
```
"""
macro formatter(expr, verbose=true)
    # Ensure the expression is a function call
    if !(expr isa Expr && expr.head == :call)
        throw(ArgumentError("@formatter must be applied to a function call"))
    end

    # Extract the function name and arguments
    func_name = expr.args[1]
    args = expr.args[2:end]

    # Check if the function is `meijerg`
    if func_name != :meijerg
        throw(ArgumentError("@formatter can only be applied to `meijerg` function calls"))
    end

    esc(:( meijerg_to_latex_and_math($(args...), verbose=$(verbose)) ))
end

"""
    meijerg_to_latex_and_math(args...; verbose=true)

Convert Meijer G-function syntax to LaTeX and Mathematica string formats.

# Arguments
- `args`: Either complete or split parameter forms of the Meijer G-function.
- `verbose`: (Optional) If true (default), prints the representations. Set to false to suppress output.

# Returns
A tuple `(plaintext_str, latex_str, mathematica_str)` containing the Plaintext, LaTeX, and Mathematica representations.
"""
function meijerg_to_latex_and_math(args...; verbose=true)
    # Normalize tuple inputs to arrays
    args = map(x -> isa(x, Tuple) ? collect(x) : x, args)

    if length(args) != 5
        throw(ArgumentError("Invalid number of arguments; expected 5"))
    end

    # Distinguish complete form (a, b, m, n, z) from split form
    if isa(args[3], Number)
        # Complete parameter form: (a, b, m, n, z)
        a, b, m, n, z = args
        a = isa(a, AbstractVector) ? collect(a) : a
        b = isa(b, AbstractVector) ? collect(b) : b
        m = Int(m)
        n = Int(n)

        a_str = isempty(a) ? "" : join(a, ", ")
        b_str = isempty(b) ? "" : join(b, ", ")

        plaintext = "G_{" * string(length(b)) * "," * string(length(a)) * "}^{" * string(m) * "," * string(n) * "}(" * string(z) * " | " * a_str * " ; " * b_str * ")"

        # LaTeX: use an array/matrix for the stacked parameters so the output is copy-pasteable
        a_tex = isempty(a) ? "" : join(map(string, a), ", ")
        b_tex = isempty(b) ? "" : join(map(string, b), ", ")
        latex = "G_{" * string(length(b)) * "," * string(length(a)) * "}^{" * string(m) * "," * string(n) * "}\\left(" * string(z) * "\\;\\middle|\\begin{matrix}" * a_tex * "\\\\" * b_tex * "\\end{matrix}\\right)"
        mathematica = "MeijerG[{{" * join(a[1:n], ", ") * "}, {" * join(a[n+1:end], ", ") * "}}, {{" * join(b[1:m], ", ") * "}, {" * join(b[m+1:end], ", ") * "}}, " * string(z) * "]"
    else
        # Split parameter form: (a_left, a_right, b_left, b_right, z)
        a_left, a_right, b_left, b_right, z = args
        a_left = isa(a_left, AbstractVector) ? collect(a_left) : a_left
        a_right = isa(a_right, AbstractVector) ? collect(a_right) : a_right
        b_left = isa(b_left, AbstractVector) ? collect(b_left) : b_left
        b_right = isa(b_right, AbstractVector) ? collect(b_right) : b_right

        a = vcat(a_left, a_right)
        b = vcat(b_left, b_right)
        m = length(b_left)
        n = length(a_left)

        a_str = isempty(a) ? "" : join(a, ", ")
        b_str = isempty(b) ? "" : join(b, ", ")

        plaintext = "G_{" * string(length(b)) * "," * string(length(a)) * "}^{" * string(m) * "," * string(n) * "}(" * string(z) * " | " * a_str * " ; " * b_str * ")"

        # LaTeX output using matrix for stacked parameters (split-form)
        a_tex = isempty(a) ? "" : join(map(string, a), ", ")
        b_tex = isempty(b) ? "" : join(map(string, b), ", ")
        latex = "G_{" * string(length(b)) * "," * string(length(a)) * "}^{" * string(m) * "," * string(n) * "}\\left(" * string(z) * "\\;\\middle|\\begin{matrix}" * a_tex * "\\\\" * b_tex * "\\end{matrix}\\right)"
        mathematica = "MeijerG[{{" * join(a_left, ", ") * "}, {" * join(a_right, ", ") * "}}, {{" * join(b_left, ", ") * "}, {" * join(b_right, ", ") * "}}, " * string(z) * "]"
    end

    if verbose
        println("\n==============================================")
        println("Plaintext: ", plaintext)
        println("\nLaTeX: ", latex)
        println("\nMathematica: ", mathematica)
        println("==============================================\n")
    end

    return (plaintext, latex, mathematica)
end
