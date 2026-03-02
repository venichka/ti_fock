module TI_quantum

using BigCombinatorics
using Nemo

export path, heaviside, Factorial_1, Binomial_1
export u, v, c_l_hyper, p_ab, p_s_fock
export op_average, g2_function, e_field_average, e_field_k_average

# function definitions go here

"""
    TI_quantum.path(; base_dir=nothing, figs_dir="figs", data_dir="data", create=true)

Returns portable output paths for figures and data.

By default, paths are created inside the project root. You can override the
root with `base_dir` or by setting the `TI_QUANTUM_OUTPUT_DIR` environment
variable.

# Returns:
- PATH_FIGS: The path for the figures.
- PATH_DATA: The path for the data.
"""
function path(; base_dir::Union{Nothing,AbstractString}=nothing,
              figs_dir::AbstractString="figs",
              data_dir::AbstractString="data",
              create::Bool=true)
    output_root = if !isnothing(base_dir)
        abspath(expanduser(base_dir))
    elseif haskey(ENV, "TI_QUANTUM_OUTPUT_DIR") &&
           !isempty(strip(ENV["TI_QUANTUM_OUTPUT_DIR"]))
        abspath(expanduser(ENV["TI_QUANTUM_OUTPUT_DIR"]))
    else
        abspath(joinpath(@__DIR__, ".."))
    end

    path_figs = abspath(joinpath(output_root, figs_dir))
    path_data = abspath(joinpath(output_root, data_dir))

    if create
        mkpath(path_figs)
        mkpath(path_data)
    end

    return joinpath(path_figs, ""), joinpath(path_data, "")
end


"""
    TI_quantum.heaviside(x)

Computes the Heaviside step function of the input argument `x`, defined as 0.5*(sign(x) + 1).

# Arguments:
- x: The input argument.

# Returns:
- The value of the Heaviside step function of `x`.
"""
function heaviside(x)
    0.5*(sign.(x) .+ 1)
end

"""
    TI_quantum.u(r)

Computes the coefficient `u` in the Bogoliubov transformation for a single time interface, given the input argument `r`.

# Arguments:
- r: The input argument.

# Returns:
- The value of the coefficient `u`.
"""
function u(r::Real)
    0.5*(1.0+r)/sqrt(r)
end


"""
    TI_quantum.v(r)

Computes the coefficient `v` in the Bogoliubov transformation for a single time interface, given the input argument `r`.

# Arguments:
- r: The input argument.

# Returns:
- The value of the coefficient `v`.
"""
function v(r::Real)
    0.5*(1.0-r)/sqrt(r)
end


"""
    TI_quantum.Factorial_1(n::Integer)::BigInt

# Description:
This function calculates the factorial of a non-negative integer `n`, returning the result as a `BigInt`. It uses memoization to save previously calculated factorials, making subsequent calls to the function faster. It also checks that the argument `n` is non-negative and throws a `DomainError` if it is not.

# See also:

    https://oeis.org/A000142
    Julia's factorial function
    FallingFactorial and RisingFactorial functions in Julia

# Important
- Cache initialization is handled automatically on module load.

# Arguments:
- n

# Returns
- n!
"""
function Factorial_1(n::Integer)::BigInt
    _ensure_factorial_cache!()
    if n < 0
        throw(DomainError(n, "argument must be nonngative"))
    end
    if BigCombinatorics._has(Factorial_1, n)
        return BigCombinatorics._get(Factorial_1, n)
    end

    start = BigCombinatorics._max_arg(Factorial_1) + 1
    for m = start:n
        val = m * BigCombinatorics._get(Factorial_1, m - 1)
        BigCombinatorics._save(Factorial_1, m, val)
    end

    return BigCombinatorics._get(Factorial_1, n)
end


"""
    TI_quantum.Factorial_1()

# Description:
This function initializes the memoization table used by `Factorial_1`. It saves the values of `0!` and `1!` as `BigInts` to the table.
"""
function Factorial_1()
    BigCombinatorics._make(Factorial_1, Integer)
    BigCombinatorics._save(Factorial_1, 0, big(1))
    BigCombinatorics._save(Factorial_1, 1, big(1))
end

function _ensure_factorial_cache!()
    table = getfield(BigCombinatorics, :_master_table)
    if !haskey(table, Factorial_1)
        Factorial_1()
    end
    nothing
end

function __init__()
    _ensure_factorial_cache!()
end

Factorial_1()


"""
    TI_quantum.Binomial_1(n::Integer, k::Integer)::BigInt

# Description:
This function calculates the binomial coefficient "n choose k", returning the result as a `BigInt`. It uses the memoized `Factorial_1` function to calculate the factorials necessary for the calculation. It requires that `0 <= k <= n`.
"""
function Binomial_1(n::Integer, k::Integer)::BigInt
    return (Factorial_1(n) / (Factorial_1(n-k) * Factorial_1(k)))
end

"""
    TI_quantum.c_l_hyper(n_1, n_2, l, r)

Computes the coefficients in the expansion of the output wavefunction over the Fock states after the time interface, using the hypergeometric function.

# Arguments:
- n_1: The number of particles in the first mode.
- n_2: The number of particles in the second mode.
- l: The difference between the number of particles in the first and second modes.
- r: The input argument.

# Returns:
- The coefficients in the expansion of the output wavefunction over the Fock states after the time interface, using the hypergeometric function.
"""
function c_l_hyper(n_1::Int, n_2::Int, l::Int, r::Real; precision=4096)
    CC = AcbField(precision)
    u_0 = big(u(r))
    v_0 = big(v(r))
    (
    1 / abs(u_0) * conj(u_0)^n_1 * (- (v_0 / u_0))^l *
    sqrt((Factorial_1(n_1 + l) / Factorial_1(n_1)) *
         (Factorial_1(n_2 + l) / Factorial_1(n_2))) *
    u_0^(-n_2) *
        hypergeometric_2f1(CC(-n_1), CC(n_2 + l + 1), CC(l + 1),
            CC(abs2(v_0) / abs2(u_0)); flags=1)
    )
end

function c_l_hyper(n_1::Int, n_2::Int, l::Int, u::Complex, v::Complex; 
                   precision=4096)
    CC = AcbField(precision)
    (
    1 / abs(u) * conj(u)^n_1 * (- (v / u))^l * 
    sqrt((Factorial_1(n_1 + l) / Factorial_1(n_1)) * 
         (Factorial_1(n_2 + l) / Factorial_1(n_2))) *
        u^(-n_2) *
        convert(ComplexF64,
                hypergeometric_2f1(CC(-n_1), CC(n_2 + l + 1), CC(l + 1),
                                   CC(abs2(v) / abs2(u)); flags=1))
    )
end

function c_l_hyper(n_1::Int, n_2::Int, l::Int, u::Real, v::Real; 
                   precision=4096)
    CC = AcbField(precision)
    return convert(Float64,
                   1 / abs(u) * u^n_1 * (- (v / u))^l * 
                   sqrt((Factorial_1(n_1 + l) / Factorial_1(n_1)) * 
                        (Factorial_1(n_2 + l) / Factorial_1(n_2))) *
                   u^(-n_2) *
                   convert(ComplexF64,
                           hypergeometric_2f1(CC(-n_1), CC(n_2 + l + 1),
                                              CC(l + 1),
                                              CC(abs2(v) / abs2(u)); flags=1))
                   )
end


"""
    TI_quantum.p_s_fock(n_1::Int, n_2::Int, s::Int, u::Complex, v::Complex; precision=4096)

Computes the probability distribution for the number of photons in modes for the ininitial Fock state using the `c_l_hyper` function.

# Arguments:
- `n_1::Int`: The number of particles in the first mode.
- `n_2::Int`: The number of particles in the second mode.
- `s::Int`: The number of photon pairs produced.
- `u::Complex`: A complex number related to the system parameters.
- `v::Complex`: A complex number related to the system parameters.
- `precision::Int`: The precision used in the calculations (default is 4096).

# Returns:
- The probability distribution for the number of photons in modes.
"""
function p_s_fock(n_1::Int, n_2::Int, s::Int, u::Complex, v::Complex;
                  precision=4096)
    coeff = c_l_hyper(n_1, n_2, s, u, v; precision=precision)
    return abs2(coeff)
end
function p_s_fock(n_1::Int, n_2::Int, s::Int, u::Real, v::Real;
                  precision=4096)
    coeff = c_l_hyper(n_1, n_2, s, u, v; precision=precision)
    return abs2(coeff)
end



"""
    `TI_quantum.p_ab(alpha::Complex, beta::Complex, n_1::Int, 
              n_2::Int, r::Real, num::Int; precision=4096)`

Computes the probability to detect `n_1` and `n_2` photons in `+k` and `-k` modes, respectively, after the time interface.

# Arguments:
- `alpha`: initial value for `+k` mode.
- `beta`: initial value for `-k` mode.
- `n_1`: The number of particles in the `+k` mode.
- `n_2`: The number of particles in the '-k mode.
- `r`: The ratio of refractive indices: `before TI / after TI`.
- `num`: Number terms in sum.

# Returns:
- The probability to detect `n_1` and `n_2` photons in `+k` and `-k` modes.
"""
function p_ab(alpha::Real, beta::Real, n_1::Int, 
              n_2::Int, r::Real, num::Int; precision=4096, phase_fraction_a::Vector{Int}=[0, 1], 
              phase_fraction_b::Vector{Int}=[0, 1])
    C = CalciumField()
    CC = AcbField(precision)
    alpha = alpha * CC(exp(C(1im)*C(pi)*phase_fraction_a[1]/phase_fraction_a[2]))
    beta = beta * CC(exp(C(1im)*C(pi)*phase_fraction_b[1]/phase_fraction_b[2]))
    (abs(convert(Float64, exp(-abs(alpha)^2.0 - abs(beta)^2.0))) *
    abs(sum([
        convert(ComplexF64, (alpha)^(n) * (beta)^((n - n_1 + n_2)) *
        (1/sqrt(Factorial_1(n) * Factorial_1(n -  n_1 + n_2))) * 
        c_l_hyper(n, n-n_1+n_2, min(n_1, n_2) - min(n, n-n_1+n_2), r))
    for n = abs(n_1 - n_2) * Int(floor(heaviside(n_1 - n_2))):abs(n_1 - n_2) * Int(floor(heaviside(n_1 - n_2)))+num]))^2 )
end


"""
`op_average(alpha::Number, beta::Number, r::Real, num::Int, op::String)`

Compute the average value of the operator over the coherent states

# Operators
- `"a_k"`
- `"a^dag_k"`
- `"a_-k"`
- `"a^dag_-k"`
- `"n_1"`
- `"n_2"`
- `"adk_adk_ak_ak"`
- `"ad-k_ad-k_a-k_a-k"`
- `"adk_ad-k_a-k_ak"`
"""
function op_average(alpha::Number, beta::Number, r::Real, num::Int, op::String)
    if op == "a_k"
        n_1 = (n_01, n_02, s_0, n) -> n
        n_2 = (n_01, n_02, s_0, n)-> n_1(n_01, n_02, s_0, n) - 1 - n_01 + n_02
        s = (n_01, n_02, s_0, n) -> (n_01 + s_0 - min(n_01, n_02) - 
                                    (n_1(n_01, n_02, s_0, n) - 1 - 
                                    min(n_1(n_01, n_02, s_0, n), 
                                        n_2(n_01, n_02, s_0, n))))
        sqrt_fun = (n_01, n_02, s_0, n) -> sqrt(n_1(n_01, n_02, s_0, n) + 
                                                s(n_01, n_02, s_0, n) - 
                                                min(n_1(n_01, n_02, s_0, n),
                                                    n_2(n_01, n_02, s_0, n)))
        lim_fun = (n_01, n_02) -> (abs(n_01 - n_02 + 1) * 
                                    Int(floor(heaviside(n_01 - n_02 + 1))))
    elseif op == "a^dag_k"
        n_1 = (n_01, n_02, s_0, n) -> n
        n_2 = (n_01, n_02, s_0, n) -> n_1(n_01, n_02, s_0, n) + 1 - n_01 + n_02
        s = (n_01, n_02, s_0, n) -> (n_01 + s_0 - min(n_01, n_02) - 
                                    (n_1(n_01, n_02, s_0, n) + 1 - 
                                    min(n_1(n_01, n_02, s_0, n), 
                                        n_2(n_01, n_02, s_0, n))))
        sqrt_fun = (n_01, n_02, s_0, n) -> sqrt(n_1(n_01, n_02, s_0, n) + 
                                                s(n_01, n_02, s_0, n) + 1 - 
                                                min(n_1(n_01, n_02, s_0, n),
                                                    n_2(n_01, n_02, s_0, n)))
        lim_fun = (n_01, n_02) -> (abs(n_01 - n_02 - 1) * 
                                    Int(floor(heaviside(n_01 - n_02 - 1))))
    elseif op == "a_-k"
        n_1 = (n_01, n_02, s_0, n) -> n
        n_2 = (n_01, n_02, s_0, n) -> n_1(n_01, n_02, s_0, n) + 1 - n_01 + n_02
        s = (n_01, n_02, s_0, n) -> (n_01 + s_0 - min(n_01, n_02) - 
                                    (n_1(n_01, n_02, s_0, n) - 
                                    min(n_1(n_01, n_02, s_0, n), 
                                        n_2(n_01, n_02, s_0, n))))
        sqrt_fun = (n_01, n_02, s_0, n) -> sqrt(n_2(n_01, n_02, s_0, n) + 
                                                s(n_01, n_02, s_0, n) - 
                                                min(n_1(n_01, n_02, s_0, n),
                                                    n_2(n_01, n_02, s_0, n)))
        lim_fun = (n_01, n_02) -> (abs(n_01 - n_02 - 1) * 
                                    Int(floor(heaviside(n_01 - n_02 - 1))))
    elseif op == "a^dag_-k"
        n_1 = (n_01, n_02, s_0, n) -> n
        n_2 = (n_01, n_02, s_0, n) -> n_1(n_01, n_02, s_0, n) - 1 - n_01 + n_02
        s = (n_01, n_02, s_0, n) -> (n_01 + s_0 - min(n_01, n_02) - 
                                    (n_1(n_01, n_02, s_0, n) - 
                                    min(n_1(n_01, n_02, s_0, n), 
                                        n_2(n_01, n_02, s_0, n))))
        sqrt_fun = (n_01, n_02, s_0, n) -> sqrt(n_2(n_01, n_02, s_0, n) + 
                                                s(n_01, n_02, s_0, n) + 1 - 
                                                min(n_1(n_01, n_02, s_0, n),
                                                    n_2(n_01, n_02, s_0, n)))
        lim_fun = (n_01, n_02) -> (abs(n_01 - n_02 + 1) * 
                                    Int(floor(heaviside(n_01 - n_02 + 1))))
    elseif op == "n_k"
        n_1 = (n_01, n_02, s_0, n) -> n
        n_2 = (n_01, n_02, s_0, n) -> n_1(n_01, n_02, s_0, n) - n_01 + n_02
        s = (n_01, n_02, s_0, n) -> s_0
        sqrt_fun = (n_01, n_02, s_0, n) -> (n_1(n_01, n_02, s_0, n) + 
                                            s(n_01, n_02, s_0, n) - 
                                            min(n_1(n_01, n_02, s_0, n),
                                                n_2(n_01, n_02, s_0, n)))
        lim_fun = (n_01, n_02) -> (abs(n_01 - n_02) * 
                                    Int(floor(heaviside(n_01 - n_02))))
    elseif op == "n_-k"
        n_1 = (n_01, n_02, s_0, n) -> n
        n_2 = (n_01, n_02, s_0, n) -> n_1(n_01, n_02, s_0, n) - n_01 + n_02
        s = (n_01, n_02, s_0, n) -> s_0
        sqrt_fun = (n_01, n_02, s_0, n) -> (n_2(n_01, n_02, s_0, n) + 
                                            s(n_01, n_02, s_0, n) - 
                                            min(n_1(n_01, n_02, s_0, n),
                                                n_2(n_01, n_02, s_0, n)))
        lim_fun = (n_01, n_02) -> (abs(n_01 - n_02) * 
                                    Int(floor(heaviside(n_01 - n_02))))
    elseif op == "adk_adk_ak_ak"
        n_1 = (n_01, n_02, s_0, n) -> n
        n_2 = (n_01, n_02, s_0, n) -> n_1(n_01, n_02, s_0, n) - n_01 + n_02
        s = (n_01, n_02, s_0, n) -> s_0
        sqrt_fun = (n_01, n_02, s_0, n) -> ((n_1(n_01, n_02, s_0, n) + 
                                            s(n_01, n_02, s_0, n) - 
                                            min(n_1(n_01, n_02, s_0, n),
                                                n_2(n_01, n_02, s_0, n))) * 
                                            (n_1(n_01, n_02, s_0, n) + 
                                            s(n_01, n_02, s_0, n) - 
                                            min(n_1(n_01, n_02, s_0, n),
                                                n_2(n_01, n_02, s_0, n))-1))
        lim_fun = (n_01, n_02) -> (abs(n_01 - n_02) * 
                                    Int(floor(heaviside(n_01 - n_02))))
    elseif op == "ad-k_ad-k_a-k_a-k"
        n_1 = (n_01, n_02, s_0, n) -> n
        n_2 = (n_01, n_02, s_0, n) -> n_1(n_01, n_02, s_0, n) - n_01 + n_02
        s = (n_01, n_02, s_0, n) -> s_0
        sqrt_fun = (n_01, n_02, s_0, n) -> ((n_2(n_01, n_02, s_0, n) + 
                                            s(n_01, n_02, s_0, n) - 
                                            min(n_1(n_01, n_02, s_0, n),
                                                n_2(n_01, n_02, s_0, n))) * 
                                            (n_2(n_01, n_02, s_0, n) + 
                                            s(n_01, n_02, s_0, n) - 
                                            min(n_1(n_01, n_02, s_0, n),
                                                n_2(n_01, n_02, s_0, n))-1))
        lim_fun = (n_01, n_02) -> (abs(n_01 - n_02) * 
                                    Int(floor(heaviside(n_01 - n_02))))
    elseif op == "ad-k_adk_ak_a-k"
        n_1 = (n_01, n_02, s_0, n) -> n
        n_2 = (n_01, n_02, s_0, n) -> n_1(n_01, n_02, s_0, n) - n_01 + n_02
        s = (n_01, n_02, s_0, n) -> s_0
        sqrt_fun = (n_01, n_02, s_0, n) -> ((n_1(n_01, n_02, s_0, n) + 
                                            s(n_01, n_02, s_0, n) - 
                                            min(n_1(n_01, n_02, s_0, n),
                                                n_2(n_01, n_02, s_0, n))) * 
                                            (n_2(n_01, n_02, s_0, n) + 
                                            s(n_01, n_02, s_0, n) - 
                                            min(n_1(n_01, n_02, s_0, n),
                                                n_2(n_01, n_02, s_0, n))))
        lim_fun = (n_01, n_02) -> (abs(n_01 - n_02) * 
                                    Int(floor(heaviside(n_01 - n_02))))
    end

    result_re = Threads.Atomic{Float64}(0)
    result_im = Threads.Atomic{Float64}(0)
    lk = ReentrantLock()
    Threads.@threads for ijk in CartesianIndices((num, num, num))
        i, j, k = Tuple(ijk)[1], Tuple(ijk)[2], Tuple(ijk)[3]
        n_01 = i - 1
        n_02 = j - 1
        if op == "a^dag_k"
            if n_01 > n_02
                s_0 = k - 1
            else
                s_0 = k
            end
        elseif op == "a^dag_-k"
            if n_02 > n_01
                s_0 = k - 1
            else
                s_0 = k
            end
        else
            s_0 = k - 1
        end
            
        for n = lim_fun(n_01, n_02):lim_fun(n_01, n_02) + num
            result_01 = (conj(alpha)^n_01*conj(beta)^n_02*alpha^n_1(n_01, n_02, s_0, n)*beta^n_2(n_01, n_02, s_0, n)) 
            
            result_02 = convert(ComplexF64, (
            (1/sqrt(Factorial_1(n_01)*Factorial_1(n_02)*
                    Factorial_1(n_1(n_01, n_02, s_0, n))*
                    Factorial_1(n_2(n_01, n_02, s_0, n)))) * 
            conj(c_l_hyper(n_01,n_02,s_0-min(n_01, n_02), r))*
            c_l_hyper(n_1(n_01, n_02, s_0, n), n_2(n_01, n_02, s_0, n), 
                        s(n_01, n_02, s_0, n)-min(n_1(n_01, n_02, s_0, n),
                                                n_2(n_01, n_02, s_0, n)), r)
            )) * sqrt_fun(n_01, n_02, s_0, n)

            result = result_01 * result_02
            
            lock(lk) do
                Threads.atomic_add!(result_re, convert(Float64, real(result)))
                Threads.atomic_add!(result_im, convert(Float64, imag(result)))
            end
        end
    end
    exp(-abs(alpha)^2 - abs(beta)^2) * (result_re[] + im*result_im[])
end


"""
`e_field_average(r::Real, k::Real, t::Real, omega::Real, a_k, a_mk)`

Compute the average value of the E field operator over the coherent states. Note that the dimentional coefficient is omitted.

"""
function e_field_average(r::Real, k::Real, t::Real, omega::Real, a_k, a_mk)
    (im .* (
        a_k .* exp(im*(k*r - omega*t)) .- conj.(a_k) .* exp(-im*(k*r - omega*t)) .+
        a_mk.* exp(im*(-k*r -omega*t)) .- conj.(a_mk) .* exp(im*(k*r + omega*t))
    ))
end


"""
`e_field_k_average(r::Real, k::Real, t::Real, omega::Real, a_k)`

Compute the average value of the E field operator over the coherent states. Note that the dimentional coefficient is omitted.

"""
function e_field_k_average(r::Real, k::Real, t::Real, omega::Real, a_k)
    (im .* (
        a_k .* exp(im*(k*r - omega*t)) .- conj.(a_k) .* exp(-im*(k*r - omega*t)) ))
end


"""
`g2_function(alpha::Number, beta::Number, r::Real, num::Int, 
                        mode::String; num_1 = num)`

Compute the g^{(2)}(0) function for 
- `a_k` operators (`mode = "k"`)
- `a_-k` operators (`mode = "-k"`)
- `a_k` and a_-k operators (`mode = "k,-k"`)

"""
function g2_function(alpha::Number, beta::Number, r::Real, num::Int, 
                        mode::String; num_1 = num)
    if mode == "k"
        numerator_g2 = op_average(alpha, beta, r, num, "adk_adk_ak_ak")
        denominator_g2 = op_average(alpha, beta, r, num_1, "n_k")
        return numerator_g2 / denominator_g2^2, denominator_g2
    elseif mode == "-k"
        numerator_g2 = op_average(alpha, beta, r, num, "ad-k_ad-k_a-k_a-k")
        denominator_g2 = op_average(alpha, beta, r, num_1, "n_-k")
        return numerator_g2 / denominator_g2^2, denominator_g2
    elseif mode == "k,-k"
        numerator_g2 = op_average(alpha, beta, r, num, "ad-k_adk_ak_a-k")
        denominator_g2_1 = op_average(alpha, beta, r, num_1, "n_-k")
        denominator_g2_2 = op_average(alpha, beta, r, num_1, "n_k")
        return numerator_g2 / denominator_g2_1 / denominator_g2_2, denominator_g2_1, denominator_g2_2
    end
end


end # module
