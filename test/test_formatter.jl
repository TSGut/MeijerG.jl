using Test
using MeijerG
import MeijerG: meijerg_to_latex_and_math

@testset "Formatter Tests" begin
    @testset "Complete Parameter Form" begin
        plain, latex, mathematica = meijerg_to_latex_and_math((1, 2), (3, 4), 1, 1, 5, verbose=false)
        @test plain == "G_{2,2}^{1,1}(5 | 1, 2 ; 3, 4)"
        @test occursin("\\begin{matrix}", latex)
        @test occursin("\\end{matrix}", latex)
        @test occursin("\\middle|", latex)
        @test occursin("G_{2,2}^{1,1}", latex)
        @test occursin("5", latex)
        @test occursin("1, 2", latex)
        @test occursin("3, 4", latex)
        @test mathematica == "MeijerG[{{1}, {2}}, {{3}, {4}}, 5]"
    end

    @testset "Split Parameter Form" begin
        plain, latex, mathematica = meijerg_to_latex_and_math((1,), (2,), (3,), (4,), 5, verbose=false)
        @test plain == "G_{2,2}^{1,1}(5 | 1, 2 ; 3, 4)"
        @test occursin("\\begin{matrix}", latex)
        @test occursin("\\end{matrix}", latex)
        @test occursin("\\middle|", latex)
        @test occursin("G_{2,2}^{1,1}", latex)
        @test occursin("5", latex)
        @test occursin("1, 2", latex)
        @test occursin("3, 4", latex)
        @test mathematica == "MeijerG[{{1}, {2}}, {{3}, {4}}, 5]"
    end

    @testset "Macro Usage" begin
        result = @formatter meijerg((1, 2), (3, 4), 1, 1, 5);
        @test result[1] == "G_{2,2}^{1,1}(5 | 1, 2 ; 3, 4)"  # plaintext
        @test occursin("\\begin{matrix}", result[2])  # LaTeX
        @test occursin("\\end{matrix}", result[2])
        @test occursin("\\middle|", result[2])
        @test occursin("G_{2,2}^{1,1}", result[2])
        @test occursin("5", result[2])
        @test occursin("1, 2", result[2])
        @test occursin("3, 4", result[2])
        @test result[3] == "MeijerG[{{1}, {2}}, {{3}, {4}}, 5]"  # Mathematica

        result = @formatter meijerg((1, 2), (3, 4), 1, 1, 5) verbose=false;
        @test result[1] == "G_{2,2}^{1,1}(5 | 1, 2 ; 3, 4)"
        @test occursin("\\begin{matrix}", result[2])
        @test occursin("\\end{matrix}", result[2])
        @test occursin("\\middle|", result[2])
        @test occursin("G_{2,2}^{1,1}", result[2])
        @test occursin("5", result[2])
        @test occursin("1, 2", result[2])
        @test occursin("3, 4", result[2])
        @test result[3] == "MeijerG[{{1}, {2}}, {{3}, {4}}, 5]"

        # Macro with split-parameter form
        result = @formatter meijerg((1,), (2,), (3,), (4,), 5);
        @test result[1] == "G_{2,2}^{1,1}(5 | 1, 2 ; 3, 4)"
        @test occursin("\\begin{matrix}", result[2])
        @test occursin("\\end{matrix}", result[2])
        @test occursin("\\middle|", result[2])
        @test occursin("G_{2,2}^{1,1}", result[2])
        @test occursin("5", result[2])
        @test occursin("1, 2", result[2])
        @test occursin("3, 4", result[2])
        @test result[3] == "MeijerG[{{1}, {2}}, {{3}, {4}}, 5]"

        # Macro with Float inputs
        result = @formatter meijerg((1.0, 2.0), (3.0, 4.0), 1, 1, 5.0);
        @test result[1] == "G_{2,2}^{1,1}(5.0 | 1.0, 2.0 ; 3.0, 4.0)"
        @test occursin("\\begin{matrix}", result[2])
        @test occursin("\\end{matrix}", result[2])
        @test occursin("\\middle|", result[2])
        @test occursin("G_{2,2}^{1,1}", result[2])
        @test occursin("5.0", result[2])
        @test occursin("1.0, 2.0", result[2])
        @test occursin("3.0, 4.0", result[2])
        @test result[3] == "MeijerG[{{1.0}, {2.0}}, {{3.0}, {4.0}}, 5.0]"

        # Macro with Complex Float inputs
        result = @formatter meijerg((1.0+2.0im, 2.0+3.0im), (3.0, 4.0), 1, 1, 5);
        @test result[1] == "G_{2,2}^{1,1}(5 | 1.0 + 2.0im, 2.0 + 3.0im ; 3.0, 4.0)"
        @test occursin("\\begin{matrix}", result[2])
        @test occursin("\\end{matrix}", result[2])
        @test occursin("\\middle|", result[2])
        @test occursin("G_{2,2}^{1,1}", result[2])
        @test occursin("5", result[2])
        @test occursin("1.0 + 2.0im", result[2])
        @test occursin("2.0 + 3.0im", result[2])
        @test occursin("3.0, 4.0", result[2])
        @test result[3] == "MeijerG[{{1.0 + 2.0im}, {2.0 + 3.0im}}, {{3.0}, {4.0}}, 5]"

        # Macro with BigFloat inputs
        result = @formatter meijerg((BigFloat(1.5), BigFloat(2.25)), (BigFloat(3.0), BigFloat(4.0)), 1, 1, BigFloat(5.0));
        @test result[1] == "G_{2,2}^{1,1}(5.0 | 1.5, 2.25 ; 3.0, 4.0)"
        @test occursin("\\begin{matrix}", result[2])
        @test occursin("\\end{matrix}", result[2])
        @test occursin("\\middle|", result[2])
        @test occursin("G_{2,2}^{1,1}", result[2])
        @test occursin("5.0", result[2])
        @test occursin("1.5, 2.25", result[2])
        @test occursin("3.0, 4.0", result[2])
        @test result[3] == "MeijerG[{{1.5}, {2.25}}, {{3.0}, {4.0}}, 5.0]"

        # Macro with Complex BigFloat inputs
        a1 = Complex(BigFloat(1.5), BigFloat(2.5))
        a2 = Complex(BigFloat(2.75), BigFloat(1.25))
        result = @formatter meijerg((a1, a2), (BigFloat(3.0), BigFloat(4.0)), 1, 1, BigFloat(5.0));
        @test result[1] == "G_{2,2}^{1,1}(5.0 | 1.5 + 2.5im, 2.75 + 1.25im ; 3.0, 4.0)"
        @test occursin("\\begin{matrix}", result[2])
        @test occursin("\\end{matrix}", result[2])
        @test occursin("\\middle|", result[2])
        @test occursin("G_{2,2}^{1,1}", result[2])
        @test occursin("5.0", result[2])
        @test occursin("1.5 + 2.5im", result[2])
        @test occursin("2.75 + 1.25im", result[2])
        @test occursin("3.0, 4.0", result[2])
        @test result[3] == "MeijerG[{{1.5 + 2.5im}, {2.75 + 1.25im}}, {{3.0}, {4.0}}, 5.0]"
    end

    @testset "Numeric Types" begin
        # Float inputs
        plain, latex, mathematica = meijerg_to_latex_and_math((1.0, 2.0), (3.0, 4.0), 1, 1, 5.0, verbose=false)
        @test plain == "G_{2,2}^{1,1}(5.0 | 1.0, 2.0 ; 3.0, 4.0)"
        @test occursin("\\begin{matrix}", latex)
        @test occursin("\\end{matrix}", latex)
        @test occursin("\\middle|", latex)
        @test occursin("G_{2,2}^{1,1}", latex)
        @test occursin("5.0", latex)
        @test occursin("1.0, 2.0", latex)
        @test occursin("3.0, 4.0", latex)
        @test mathematica == "MeijerG[{{1.0}, {2.0}}, {{3.0}, {4.0}}, 5.0]"

        # Complex Float inputs
        plain, latex, mathematica = meijerg_to_latex_and_math((1.0+2.0im, 2.0+3.0im), (3.0, 4.0), 1, 1, 5, verbose=false)
        @test plain == "G_{2,2}^{1,1}(5 | 1.0 + 2.0im, 2.0 + 3.0im ; 3.0, 4.0)"
        @test occursin("\\begin{matrix}", latex)
        @test occursin("\\end{matrix}", latex)
        @test occursin("\\middle|", latex)
        @test occursin("G_{2,2}^{1,1}", latex)
        @test occursin("5", latex)
        @test occursin("1.0 + 2.0im", latex)
        @test occursin("2.0 + 3.0im", latex)
        @test occursin("3.0, 4.0", latex)
        @test mathematica == "MeijerG[{{1.0 + 2.0im}, {2.0 + 3.0im}}, {{3.0}, {4.0}}, 5]"

        # BigFloat inputs
        plain, latex, mathematica = meijerg_to_latex_and_math((BigFloat(1.5), BigFloat(2.25)), (BigFloat(3.0), BigFloat(4.0)), 1, 1, BigFloat(5.0), verbose=false)
        @test plain == "G_{2,2}^{1,1}(5.0 | 1.5, 2.25 ; 3.0, 4.0)"
        @test occursin("\\begin{matrix}", latex)
        @test occursin("\\end{matrix}", latex)
        @test occursin("\\middle|", latex)
        @test occursin("G_{2,2}^{1,1}", latex)
        @test occursin("5.0", latex)
        @test occursin("1.5, 2.25", latex)
        @test occursin("3.0, 4.0", latex)
        @test mathematica == "MeijerG[{{1.5}, {2.25}}, {{3.0}, {4.0}}, 5.0]"

        # Complex BigFloat inputs
        a1 = Complex(BigFloat(1.5), BigFloat(2.5))
        a2 = Complex(BigFloat(2.75), BigFloat(1.25))
        plain, latex, mathematica = meijerg_to_latex_and_math((a1, a2), (BigFloat(3.0), BigFloat(4.0)), 1, 1, BigFloat(5.0), verbose=false)
        @test plain == "G_{2,2}^{1,1}(5.0 | 1.5 + 2.5im, 2.75 + 1.25im ; 3.0, 4.0)"
        @test occursin("\\begin{matrix}", latex)
        @test occursin("\\end{matrix}", latex)
        @test occursin("\\middle|", latex)
        @test occursin("G_{2,2}^{1,1}", latex)
        @test occursin("5.0", latex)
        @test occursin("1.5 + 2.5im", latex)
        @test occursin("2.75 + 1.25im", latex)
        @test occursin("3.0, 4.0", latex)
        @test mathematica == "MeijerG[{{1.5 + 2.5im}, {2.75 + 1.25im}}, {{3.0}, {4.0}}, 5.0]"
    end
end;
