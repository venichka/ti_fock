using Test
using TI_quantum

canonical_path(p::AbstractString) = rstrip(normpath(p), ['/', '\\'])

@testset "TI_quantum smoke tests" begin
    @test heaviside(-1.0) == 0.0
    @test heaviside(0.0) == 0.5
    @test heaviside(2.0) == 1.0

    @test Factorial_1(0) == big(1)
    @test Factorial_1(5) == big(120)
    @test Binomial_1(6, 2) == big(15)

    r = 0.4
    @test isapprox(abs2(u(r)) - abs2(v(r)), 1.0; atol=1e-12)

    coeff = c_l_hyper(1, 1, 0, r)
    @test isfinite(abs2(convert(ComplexF64, coeff)))

    prob0 = p_s_fock(0, 0, 0, u(r), v(r))
    @test isfinite(prob0)
    @test 0.0 <= prob0 <= 1.0
end

@testset "Path management" begin
    mktempdir() do tmpdir
        figs, data = path(base_dir=tmpdir)
        @test isdir(figs)
        @test isdir(data)
        @test canonical_path(figs) == canonical_path(joinpath(tmpdir, "figs"))
        @test canonical_path(data) == canonical_path(joinpath(tmpdir, "data"))
    end

    mktempdir() do tmpdir
        withenv("TI_QUANTUM_OUTPUT_DIR" => tmpdir) do
            figs, data = path()
            @test isdir(figs)
            @test isdir(data)
            @test canonical_path(figs) == canonical_path(joinpath(tmpdir, "figs"))
            @test canonical_path(data) == canonical_path(joinpath(tmpdir, "data"))
        end
    end
end
