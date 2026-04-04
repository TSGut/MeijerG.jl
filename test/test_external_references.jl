@testset "External-reference checks" begin
    @testset "Reference values documented by mpmath" begin
        ref1 = 271.46290321152464592 - 703.03330399954820169im
        val1 = meijerg((), (1 + 1im,), (1,), (1,), 3 + 4im)
        # Complex reference 1 matches mpmath.
        @test val1 ≈ ref1 rtol=1e-14

        ref2 = -1.338096165935754898687431
        val2 = meijerg((-3.0,), (-0.5,), (-1.0,), (-2.5,), -0.5)
        # Real reference 2 matches mpmath.
        @test val2 ≈ ref2 rtol=1e-14

        ref3 = -(pi + 4) / (4pi)
        val3 = meijerg((-3.0,), (-0.5,), (-1.0,), (-2.5,), -1.0)
        # Real reference 3 matches closed form.
        @test val3 ≈ ref3 rtol=1e-14
    end

    @testset "Additional hard-coded references from mpmath" begin
        # Offline references generated with mpmath (high precision) and hard-coded
        # so tests do not depend on PyCall/mpmath at runtime.
        cases = [
            (a = (0.35, 1.1), b = (0.2, 0.7, 1.6), m = 1, n = 1, z = 0.6,
             ref = ComplexF64(-0.07241398474957796, 0.0)),
            (a = (0.45, 1.35), b = (0.12, 0.82, 1.74), m = 2, n = 1, z = 0.4,
             ref = ComplexF64(-0.32550031518980145, 0.0)),
            (a = (0.25, 0.8, 1.4), b = (0.15, 0.9), m = 1, n = 2, z = 1.8,
             ref = ComplexF64(0.17286201291248644, 0.0)),
            (a = (0.3, 0.95, 1.75), b = (0.2, 1.1), m = 1, n = 2, z = 2.4,
             ref = ComplexF64(-0.31897389835335743, 0.0)),
            (a = (0.45, 0.95, 1.55, 2.05), b = (0.12, 0.72, 1.22, 1.82), m = 1, n = 1, z = 0.9995,
             ref = ComplexF64(0.05214557679210809, 0.0)),
            (a = (0.4, 1.2), b = (0.15, 0.85, 1.75), m = 1, n = 1, z = 0.7 + 0.4im,
             ref = ComplexF64(-0.13263285540080822, -0.05383827730629095)),
            (a = (0.3, 0.9, 1.7), b = (0.2, 1.1), m = 1, n = 2, z = 1.4 + 0.6im,
             ref = ComplexF64(-0.27314274891696535, -0.10181497137331164)),
            (a = (0.45, 0.95, 1.55, 2.05), b = (0.12, 0.72, 1.22, 1.82), m = 1, n = 1, z = 0.92 + 0.08im,
             ref = ComplexF64(0.0383826526763178, 0.004354960240970144)),
            (a = (0.45, 0.95, 1.55, 2.05), b = (0.12, 0.72, 1.22, 1.82), m = 1, n = 1, z = 0.85 + 0.2im,
             ref = ComplexF64(0.03409281432680096, 0.00790153446472579)),
            (a = (0.45, 0.95, 1.55, 2.05), b = (0.12, 0.72, 1.22, 1.82), m = 1, n = 1, z = 0.98 + 0.02im,
             ref = ComplexF64(0.04374817674982247, 0.002495968608373242)),
        ]

        for case in cases
            value = meijerg(case.a, case.b, case.m, case.n, case.z)
            value_complex = value isa Real ? ComplexF64(value, 0.0) : ComplexF64(real(value), imag(value))
            @test value_complex ≈ case.ref rtol=1e-14 atol=1e-14
        end
    end

    @testset "BigFloat and Complex{BigFloat} hard-coded references from mpmath" begin
        # Offline references generated with mpmath at high precision and hard-coded
        # so this testset remains dependency-free at runtime.
        cases = [
            # Confluent
            (a = (big"0.35",),
             b = (big"0.10", big"1.10", big"1.85"),
             m = 2, n = 1, z = big"0.35",
             ref_re = "-0.1093138920999586896350449922954909715161896033882386207762090523816375124147",
             ref_im = "0.0"),
            (a = (big"1.10", big"2.10", big"3.00"),
             b = (big"0.20",),
             m = 1, n = 2, z = big"1.60",
             ref_re = "-82.05857615581717247370226927311796843433564488595643114064177906060640717272",
             ref_im = "0.0"),

            # Near confluent
            (a = (big"0.35",),
             b = (big"0.10", big"1.1000001", big"1.85"),
             m = 2, n = 1, z = big"0.35",
             ref_re = "-0.1093139037934360278622324154402682629361097440223195555130725179738909263458",
             ref_im = "0.0"),
            (a = (big"1.10", big"2.1000001", big"3.00"),
             b = (big"0.20",),
             m = 1, n = 2, z = big"1.60",
             ref_re = "-82.05865525827483849173693392241277805935735038813621069318681279645174470165",
             ref_im = "0.0"),

            # Near |z| = 1
            (a = (big"0.45", big"0.95", big"1.55", big"2.05"),
             b = (big"0.12", big"0.72", big"1.22", big"1.82"),
             m = 1, n = 1, z = big"0.9995",
             ref_re = "0.05214557679210786999980402035116356498445789864628528997202875992983155105482",
             ref_im = "0.0"),
            (a = (big"0.45", big"0.95", big"1.55", big"2.05"),
             b = (big"0.12", big"0.72", big"1.22", big"1.82"),
             m = 1, n = 1, z = Complex{BigFloat}(big"0.92", big"0.08"),
             ref_re = "0.03838265267631779050589642694654641427488561739468950322314982311904307020528",
             ref_im = "0.004354960240970137183535784693852029705695565487110648183868603908862938755014"),

            # Away from |z| = 1
            (a = (big"0.45", big"1.35"),
             b = (big"0.12", big"0.82", big"1.74"),
             m = 2, n = 1, z = big"0.40",
             ref_re = "-0.3255003151898014601890356241450776034565097395173994745470170907148008138558",
             ref_im = "0.0"),
            (a = (big"0.4", big"1.2"),
             b = (big"0.15", big"0.85", big"1.75"),
             m = 1, n = 1, z = Complex{BigFloat}(big"0.7", big"0.4"),
             ref_re = "-0.1326328554008082426306587393979191010061599004039250394697567147121110046468",
             ref_im = "-0.0538382773062909861453220750798390878588485494925070798279009961591681808435"),
            (a = (big"0.3", big"0.9", big"1.7"),
             b = (big"0.2", big"1.1"),
             m = 1, n = 2, z = Complex{BigFloat}(big"1.4", big"0.6"),
             ref_re = "-0.2731427489169650947274180635175955982893950168042290036921214065286151417604",
             ref_im = "-0.1018149713733116317874182081675873463740221372932547405581842600310849726563"),
        ]

        for case in cases
            reference = Complex{BigFloat}(parse(BigFloat, case.ref_re), parse(BigFloat, case.ref_im))
            value = meijerg(case.a, case.b, case.m, case.n, case.z)
            value_complex = value isa Real ? Complex{BigFloat}(BigFloat(value), big"0") : Complex{BigFloat}(BigFloat(real(value)), BigFloat(imag(value)))
            relerr = abs(value_complex - reference) / max(abs(reference), eps(BigFloat))
            @test relerr < big"1e-30"
        end
    end
end
