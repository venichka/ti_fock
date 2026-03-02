### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 30e8847a-c88e-4f68-8e19-4a2044569ae3
begin
    using Pkg
    PROJECT_TOML = Base.current_project(@__DIR__)
    isnothing(PROJECT_TOML) && error("Project.toml not found from notebook location.")
    PROJECT_DIR = dirname(PROJECT_TOML)
    Pkg.activate(PROJECT_DIR)
    Pkg.instantiate()

    # Load local package code directly from the repository root.
    if !isdefined(Main, :TI_quantum)
        include(normpath(joinpath(PROJECT_DIR, "..", "src", "TI_quantum.jl")))
    end
end

# ╔═╡ 134acf44-bc13-4744-916b-52b1339d5936
begin
	using Plots, CairoMakie, WGLMakie, LaTeXStrings
	using BenchmarkTools
	using HDF5, FileIO, Printf, LaTeXStrings
	using Nemo, ClassicalOrthogonalPolynomials
	using StaticArrays, CoordinateTransformations
	using BigCombinatorics
	using QuantumOptics
	using Revise
	using .TI_quantum

	# initialize the table for factorials
	Factorial_1()

	NMAX = 20
	PATH_FIGS, PATH_DATA = path()
	mkpath(joinpath(PATH_FIGS, "Figs_fock"))
end

# ╔═╡ d5e9d066-8aac-4993-a8c1-b990998be10a
begin
	using PlutoUI
	slider_u = @bind u_percent PlutoUI.Slider(0:0.05:1, show_value=true, default=.7)
	slider_v = @bind v_percent PlutoUI.Slider(0:0.05:1, show_value=true, default=.7)
	slider_t = @bind t_val PlutoUI.Slider(0:0.2:20, show_value=true, default=6)
	slider_th = @bind thickness PlutoUI.Slider(1:1:10, show_value=true, default=4)
	slider_l = @bind limit PlutoUI.Slider(2.0:0.5:100.0, show_value=true, default=2.0)
	slider_r = @bind r_value PlutoUI.Slider(0.001:0.005:0.999, show_value=true, default=2.0)
end

# ╔═╡ f19b8814-28f3-11ef-37ed-cb82cc7cdffe
md"""
# Plots for Fock state paper
"""

# ╔═╡ 53ed7c53-2395-4a5b-b479-14e7ae7181cb
md"""
### Functions: ``c_s(n_1, n_2)``
"""

# ╔═╡ 0d1115cd-a481-4bff-a5ee-68a9cf4dde42
begin
	"""
	Doesn't work properly for s < 0
	"""
	function c_s_nn_jacobi(n::Int, s::Int, u::Complex, v::Complex)
	    u_abs2 = abs2(u)
	    v_abs2 = abs2(v)
	    jacobi_term = jacobip(n, s, 0, 1 - 2 * v_abs2 / u_abs2)
	    result = (-v/u)^s * (conj(u)/u)^n * jacobi_term / abs(u)
	    return result
	end
	function c_s_nn_jacobi(n::Int, s::Int, u::Real, v::Real)
	    u_abs2 = u^2
	    v_abs2 = v^2
	    jacobi_term = jacobip(n, s, 0, 1 - 2 * v_abs2 / u_abs2)
	    result = (-v/u)^s * jacobi_term / abs(u)
	    return result
	end
	"""
	Seems to work fine
	"""
	function c_0_nn_legendre(n::Int, u::Complex, v::Complex)
	    u_abs2 = u^2
	    v_abs2 = v^2
	    legendre_term = legendrep(n, 1 - 2 * v_abs2 / u_abs2)
	    result = (conj(u)/u)^n * legendre_term / abs(u)
	    return result
	end
	function c_0_nn_legendre(n::Int, u::Real, v::Real)
	    u_abs2 = u^2
	    v_abs2 = v^2
	    legendre_term = legendrep(n, 1 - 2 * v_abs2 / u_abs2)
	    result = legendre_term / abs(u)
	    return result
	end
	function p_fock(c_s::Number)::Float64
	    return abs2(c_s)
	end
end

# ╔═╡ d78dbee3-f040-4753-82e9-39362d4204a1
md"""
### Hyperbola functions
"""

# ╔═╡ 216474fc-a9fc-4629-b073-529c926e8e37
begin
	function f1_h(r::Float64)::SVector{2, Float64}
	    return SVector(0.0, coth((1 - r) * acoth(3 / 2)) - 1)
	end
	function f2_h(r::Float64)::SVector{2, Float64}
	    return SVector(1.0, -0.5)
	end
	function t_a(r::Float64)
	    f1 = f1_h(r)
	    f2 = f2_h(r)
	    return AffineMap(SMatrix{2,2,Float64}(hcat(f1, f2)), SVector(0.0, 0.0))
	end
	const t_emanuele = AffineMap(SMatrix{2,2,Float64}(0.0, 1.0, 1.0, 1.0), SVector(0.0, 0.0))
	function hyperbola_impl(f0::SVector{2,Float64}, f1::SVector{2,Float64}, f2::SVector{2,Float64}, x::Float64, y::Float64)::Float64
	    return det(SMatrix{2,2,Float64}(x - f0[1], y - f0[2], f2[1], f2[2]))^2 -
	           det(SMatrix{2,2,Float64}(f1[1], f1[2], x - f0[1], y - f0[2]))^2 -
	           det(SMatrix{2,2,Float64}(f1[1], f1[2], f2[1], f2[2]))^2
	end
	function hyperbola_vect(t::Float64, f0::SVector{2,Float64}, z::Float64, d1::Float64, r::Float64)::SVector{3,Float64}
	    return SVector(f0[1] + d1 * sinh(t),
	                   f0[2] - d1 * (coth(acoth(3/2)*(r - 1))+1)*cosh(t) - d1/2*sinh(t),
	                   z)
	end
	function t0_h(r::Float64)::Float64
	    f1 = f1_h(r)
	    f2 = f2_h(r)
	    return 0.25 * log(norm(f1 - f2)^2 / norm(f1 + f2)^2)
	end
	function vertex_h(r::Float64)::SVector{2,Float64}
	    t0 = t0_h(r)
	    f1 = f1_h(r)
	    f2 = f2_h(r)
	    return f1 * cosh(t0) + f2 * sinh(t0)
	end
end

# ╔═╡ 19713970-9020-4b20-9e1e-734dbfeb37e5
md"""
### Average value functions
"""

# ╔═╡ 156068dc-baa5-448a-9fe9-44d058650968
begin
function n_fock_sym(n::Int, v::Number)::Float64
    return n + (2 * n + 1) * abs2(v)
end
function n_fock_k_gen(n::Int, m::Int, v::Number)::Float64
    return n + (n + m + 1) * abs2(v)
end
function n_fock_mk_gen(n::Int, m::Int, v::Number)::Float64
    return m + (n + m + 1) * abs2(v)
end
function n_fock_n0(n::Int, v::Number)::Float64
    return n + (n + 1) * abs2(v)
end
function n_fock_k(n1::Int, n2::Int, u::Number, v::Number, num::Int)::Float64
    total = 0.0
    for s in 0:num
        total += (n1 + s - min(n1, n2)) * abs2(c_l_hyper(n1, n2, s - min(n1, n2), u, v))
    end
    return total
end
function n_fock_mk(n1::Int, n2::Int, u::Number, v::Number, num::Int)::Float64
    total = 0.0
    for s in 0:num
        total += (n2 + s - min(n1, n2)) * abs2(c_l_hyper(n1, n2, s - min(n1, n2), u, v))
    end
    return total
end
function entropy_vn_fock(n1::Int, n2::Int, u::Number, v::Number, num::Int; 
					 	 precision=4096)::Float64
    entropy = 0.0
    min_n1_n2 = min(n1, n2)
    for s in 0:num
        coef = c_l_hyper(n1, n2, s - min_n1_n2, u, v; precision=precision)
        abs2_coef = abs2(coef)
        if abs2_coef > 0
            entropy -= abs2_coef * log(abs2_coef)
        end
    end
    return entropy
end
function entropy_vn_00(u::Number, v::Number)::Float64
    ratio = abs2(v / u)
    if ratio == 0
        return 0.0
    elseif ratio == 1
        return Inf
    else
        term1 = -log(1 - ratio)
        term2 = (ratio * log(ratio)) / (1 - ratio)
        return term1 - term2
    end
end

"""
    beta(r::Real) -> Float64

Compute the inverse thermal energy ``\\beta = (k_B T)^{-1}`` for ''vacuum'' photons given the TI contrast `r`, assuming ``\\hbar = \\omega_0 = 1``.

# Arguments:
- `r`: TI contrast (``r = n_0 / n_1 = \\omega_1 / \\omega_0``)

# Returns:
- The inverse thermal energy ``\\beta`` as a Float64.
"""
function beta(r::Real)::Float64
    return (1 / (r)) * log((abs2(v(r)) + 1) / abs2(v(r)))
end

"""
    temp_vac(r::Real) -> Float64

Compute the thermal energy ``k_B T`` for ''vacuum'' photons given the TI contrast `r`, assuming ``\\hbar = \\omega_0 = 1``.

# Arguments:
- `r`: TI contrast (``r = n_0 / n_1 = \\omega_1 / \\omega_0``)

# Returns:
- The thermal energy ``k_B T`` as a Float64.
"""
function temp_vac(r::Real)::Float64
    return r / log((abs2(v(r)) + 1) / abs2(v(r)))
end

function g2_k(n_1::Int, n_2::Int, u::Number, v::Number, num::Int; 
			  precision=4096)::Float64
    numerator = 0.0
	denominator = 0.0
    min_n1_n2 = min(n_1, n_2)
    for s in 0:num
        coef = c_l_hyper(n_1, n_2, s - min_n1_n2, u, v; precision=precision)
        abs2_coef = abs2(coef)
		term = (n_1 + s - min_n1_n2) * (n_1 + s - min_n1_n2 - 1)
		numerator += abs2_coef * term
		denominator += (n_1 + s - min_n1_n2) * abs2_coef
    end
    return numerator / abs2(denominator)
end

function g2_cross(n_1::Int, n_2::Int, u::Number, v::Number, num::Int; 
			  	  precision=4096)::Float64
    numerator = 0.0
	denominator_k = 0.0
	denominator_mk = 0.0
    min_n1_n2 = min(n_1, n_2)
    for s in 0:num
        coef = c_l_hyper(n_1, n_2, s - min_n1_n2, u, v; precision=precision)
        abs2_coef = abs2(coef)
		term = (n_1 + s - min_n1_n2) * (n_2 + s - min_n1_n2)
		numerator += abs2_coef * term
		denominator_k += (n_1 + s - min_n1_n2) * abs2_coef
		denominator_mk += (n_2 + s - min_n1_n2) * abs2_coef
    end
    return numerator / (denominator_k * denominator_mk)
end
	
"""
    g2_sym(n::Int, v::Number) -> Float64

Compute the second-order correlation function ``g^{(2)}(0)`` for a symmetric case (equal initial number of photons in modes).
"""
function g2_sym(n::Int, v::Number)::Float64
    n_vac = abs2(v)
    numerator = (6n^2 + 6n + 2) * n_vac^2 + (6n^2 + 2n) * n_vac + n^2 - n
    denominator = ((1 + 2n) * n_vac + n)^2
    return numerator / denominator
end

"""
    g2_sym_cross(n::Int, v::Number) -> Float64

Compute the cross second-order correlation function  ``g^{(2)}_{\\mathrm{cross}}(0)`` for a symmetric case (equal initial number of photons in modes).

"""
function g2_sym_cross(n::Int, v::Number)::Float64
    n_vac = abs2(v)
    numerator = (6n^2 + 6n + 2) * n_vac^2 + (6n^2 + 4n + 1) * n_vac + n^2
    denominator = ((1 + 2n) * n_vac + n)^2
    return numerator / denominator
end

"""
    g2_n0(n::Int, v::Number) -> Float64

Compute the second-order correlation function ``g^{(2)}(0)`` for a |n,0> initial state.
"""
function g2_n0(n::Int, v::Number)::Float64
    n_vac = abs2(v)
    numerator = (n+2) * (n+1) * n_vac^2 + 2*n*(n+1) * n_vac + n*(n-1)
    denominator = ((1 + n) * n_vac + n)^2
    return numerator / denominator
end

"""
    g2_nm(n::Int, m::Int, v::Number) -> Float64

Compute the second-order correlation function ``g^{(2)}_k(0)`` for a |n,m> initial state.
"""
function g2_nm(n::Int, m::Int, v::Number)::Float64
    n_vac = abs2(v)
    numerator = (n^2+m^2+3*(n+m)+4*n*m+2)*n_vac^2 + (2*n^2+4*n*m+2*n)*n_vac + n*(n-1)
    denominator = ((1 + n + m) * n_vac + n)^2
    return numerator / denominator
end

"""
	var_n_nm_k(n::Int, m::Int, v::Number) -> Float64

Compute variance of the number of photons in k-mode for a |n,m> initial state.
"""
function var_n_nm(n::Int, m::Int, v::Number)::Float64
    return n_fock_k_gen(n, m, v)^2 * (g2_nm(n, m, v) - 1) + n_fock_k_gen(n, m, v)
end

	
end

# ╔═╡ a2e10c75-44ff-442e-b9aa-30996ee25be7
md"""
### Asymptote functions
"""

# ╔═╡ d6a16977-c42e-417f-8829-8fb17e200422
begin

"""
    c0_nn_inf(n::Int, u::Number, v::Number)

Compute the asymptotic form of the coefficient for symmetric modes in the limit when n approaches infinity and s = 0.

# Arguments:
- n: The number photons in modes.
- u: Bogoliubov coefficient.
- v: Bogoliubov coefficient.

# Returns:
- The computed coefficient.
"""
function c0_nn_inf(n::Int, u::Number, v::Number)
    pi = Float64(π)
    coeff = -(1 / (16 * n^(3/2) * sqrt(pi) * abs(v)^(3/2))) * (conj(u)/u)^n
    term1 = abs(u)^2 * sin(1/4 * (pi + (2 - 4n) * acos(1 - 2abs(v)^2 / abs(u)^2)))
    term2 = 2 * (1 - 8n) * abs(v) * sin(1/4 * (pi + (2 + 4n) * acos(1 - 2abs(v)^2 / abs(u)^2)))
    return coeff * (term1 + term2)
end

"""
    cs_nn_inf(n::Int, s::Int, u::Number, v::Number)

Compute the asymptotic form of the coefficient for symmetric modes as n approaches infinity and ``s \neq 0``.

# Arguments:
- n: The number of photons in modes.
- s: The number of photon pairs created.
- u: Bogoliubov coefficient.
- v: Bogoliubov coefficient.

# Returns:
- The computed coefficient.
"""
function cs_nn_inf(n::Int, s::Int, u::Number, v::Number)
    pi = Float64(π)
    coeff = 1 / (16 * n^(3/2) * sqrt(pi)) * abs(v)^(-(3/2) - s) * (conj(u)/u)^n * (-(v/sign(u)))^s
    term1 = 2 * (-3 + 8n) * abs(v) * sin(1/4 * (pi - 2pi*s + 2 * (1 + 2n + s) * acos(1 - 2abs(v)^2 / abs(u)^2)))
    term2 = u * conj(u) * ((-1 + s) * s * sin(1/4 * (pi + 2pi*s - 2 * (-1 + 2n + s) * acos(1 - 2abs(v)^2 / abs(u)^2)))
                           + 2 * s^2 * sin(1/4 * (pi + 2pi*s - 2 * (1 + 2n + s) * acos(1 - 2abs(v)^2 / abs(u)^2)))
                           + (-1 + s + s^2) * sin(1/4 * (pi + 2pi*s - 2 * (3 + 2n + s) * acos(1 - 2abs(v)^2 / abs(u)^2))))
    return coeff * (term1 + term2)
end

"""
    cmn_nn(n::Int, u::Number, v::Number)

Compute the coefficient for n_1 = n_2 = -s = n.
"""
function cmn_nn(n::Int, u::Number, v::Number)
    return v^(-n) * abs(u)^(-1 - 2n) * abs(v)^(2n) * conj(u)^n
end

"""
    cs_00(s::Int, u::Number, v::Number)

Compute the coefficient for a vacuum initial state.
"""
function cs_00(s::Int, u::Number, v::Number)
    return (-(v/u))^s / abs(u)
end

"""
    p_0_nn_approx(n::Int, u::Complex, v::Complex) -> Float64

Compute the approximate probability distribution of the photon numbers for s=0, n1=n2=n as n->∞.

"""
function p_0_approx(n::Int, u::Number, v::Number)::Float64
    abs_u2 = abs2(u)
    abs_v2 = abs2(v)
    cos_inverse_arg = 1 - 2 * abs_v2 / abs_u2
    angle = acos(cos_inverse_arg)
    sin_term = sin((2 * n + 1) * angle)
    return (1 + sin_term) / (2 * π * n * abs(v))
end

end

# ╔═╡ ecae7206-cc04-4a03-89c4-9c21672ad0b1
md"""
### Misc functions
"""

# ╔═╡ 2da5b512-8cb1-4fae-a4af-371361ccb0f6
begin
	function plot_coordinate_mesh(x::Function, y::Function, ax::Makie.Axis; 
								  t_val::Number = 10, u_length::Int = 10, 
								  v_length::Int = 10, t_length::Int = 100, 
								  u_lim::Number = 5, v_lim::Number = 5, 
								  t_val_v::Number = 10, kwargs...)
	
	    u_range = range(-u_lim, u_lim, length=u_length)
	    v_range = range(-v_lim, v_lim, length=v_length)
	    
	    t_range = range(-t_val, t_val, length=t_length)
		t_range_v = range(-t_val_v, t_val_v, length=t_length)
	
	    # Transform to xy-plane
	    
	    lines!(ax, [(x(u_range[1], v), y(u_range[1], v)) for v in t_range]; kwargs...)
	    for i in u_range[2:end]
	        lines!(ax, [(x(i, v), y(i, v)) for v in t_range_v]; kwargs...)
	    end
	    for i in v_range
	        lines!(ax, [(x(u, i), y(u, i)) for u in t_range]; kwargs...)
	    end
	end
end

# ╔═╡ a900ceb7-5e4f-4ddf-a7fc-a0d1e7f52819
let
	CairoMakie.activate!()

	r_range = 10.0.^(range(-0.8, 0.8, 1000))
	f = Figure()
	ax = Axis(f[1, 1],
	    xlabel = L"r = \omega_1 / \omega_0",
	    ylabel = L"k_B T / \hbar \omega_0",
		xscale=log10,
		yscale=log10,
		xlabelsize=20,
		ylabelsize=20,
	)
	Makie.lines!(ax, r_range, map(r -> 1 / beta(r), r_range), color = :tomato, linestyle = :solid, linewidth=2.0)
	Makie.lines!(ax, r_range[1:16*length(r_range) ÷ 30], map(r -> 0.25, r_range[1:16*length(r_range) ÷ 30]), color = :tomato, linestyle = :dash, linewidth=1)
	Makie.lines!(ax, r_range[14*length(r_range) ÷ 30:end], map(r -> r^2 / 4, r_range[14*length(r_range) ÷ 30:end]), color = :tomato, linestyle = :dash, linewidth=1)
	Makie.text!(L"\left(k_B T / \hbar \omega_0 \right)_{r \rightarrow 0} = 1/4", 
		position = (0.2, 0.26),
		fontsize=18)
	Makie.text!(L"\left(k_B T / \hbar \omega_0 \right)_{r \rightarrow \infty} = r^2/4", 
		position = (2.0, 1.3),
		fontsize=18,
		rotation=0.25*pi
	)
	# save(PATH_FIGS*"Figs_fock/Fig_1.pdf", f)
	f
end

# ╔═╡ ed403f68-031e-4164-af23-23f8b4632691
let
	WGLMakie.activate!()
	nmax = 50
	r = 0.02
	x = y = 0:nmax
	z = -nmax:1:nmax
	probability = [(s >= -min(n_1, n_2)) ? p_fock(c_l_hyper(n_1, n_2, s, u(r), v(r))) : 0.0 for n_1 = x, n_2 = y, s = z]

	# colormap with transparency at the beginning
	cmap = :plasma
	cmap_hm = :inferno
	n = 101
	g(x) = 0.5*(tanh(15*(x-0.1)) + 1)
	g_hm(x) = 0.5*(tanh(10*(x-0.1)) + 1)
	alphas = [(x == 0) ? 0.0 : g(x) for x in range(0, 1, length = n)]
	alphas_hm = [(x == 0) ? 0.0 : g_hm(x) for x in range(0, 1, length = n)]
	cmap_alpha = resample_cmap(cmap, n; alpha = alphas)
	cmap_alpha_hm = resample_cmap(cmap_hm, n; alpha = alphas_hm)
	# cmap = cgrad([:blue, :green, :yellow, :red], alpha=[0.0, 1.0, 1.0, 1.0], rev=false)


	prob_nz = probability[probability .> 0.0]
	n_levels = 8
	levels = 10.0.^range(log10.(0.002), log10.(maximum(prob_nz)), n_levels)

	fig, ax, plt = Makie.volume((0, nmax), (0, nmax), (-nmax, nmax), probability; 
		interpolate=true,
		algorithm = :mip, 
		absorption=4f0,
		colormap=cmap_alpha,
    	nan_color=:transparent,
    	transparency=true,
		figure = (;
        	size = (600, 600)
        ),
    	axis = (;
	        type = Axis3,
	        perspectiveness = 0.5,
			aspect = :data,
	        azimuth = 0.32*pi,
	        elevation = 0.2*pi,
			# viewmode = :fit,
	        xlabel = L"n_k",
	        ylabel = L"n_{-k}",
	        zlabel = L"s \textrm{, photon pairs}",
			xlabelsize=24,
			ylabelsize=24,
			zlabelsize=24,
			limits=(0, 50, 0, 50, nothing, nothing),
	        )
	)
	Makie.heatmap!((0, nmax), (0, nmax), probability[:,:, 50];
		colormap=cmap_alpha_hm,
		interpolate=true,
		# levels = 50,
		transformation=(:xy, minimum(z)),
		)
	Makie.heatmap!((0, 50), (-50, 50), [probability[i, i, j] for i = 1:51, j = 1:101];
		colormap=cmap_alpha_hm,
		interpolate=true,
		# levels = levels,
		transformation=(:xz, minimum(y)),
		)
	Colorbar(fig[1, 2], plt; 
		label = L"p_s(n_k, n_{-k})", 
		height = Relative(0.5), 
		labelsize=24,
	)
	# save(PATH_FIGS*"Figs_fock/Fig_2.png", fig, px_per_unit=4)
	fig
end

# ╔═╡ 05bf53c8-22c5-43ea-90cd-338f4cb230e7
Plots.plot(range(0,1,1000), [map(x -> 0.5*(tanh(15*(x-0.1)) + 1), range(0,1,1000)),
							 map(x -> x^(1/3), range(0,1,1000))])

# ╔═╡ be504ed6-e88b-4150-a7ba-64788021d5b7
let
	CairoMakie.activate!()

	nmax = 50
	r = 0.02
	x = 0:nmax
	y = -nmax:1:nmax
	# z = vec([(n, n + s, (s >= 0) ? p_fock(c_l_hyper(n, 0, s, u(r), v(r))) : NaN) for n = x, s = y])
	z = vec([(n, n + s, (s >= -n) ? p_fock(c_l_hyper(n, n, s, u(r), v(r))) : NaN) for n = x, s = y])

	xs = [i[1] for i in z]
	ys = [i[2] for i in z]
	zs = [i[3] for i in z]

	# mesh functions
	r_value=0.008
	f0 = SVector(0.0, 0.0)
	x(u, v) = v*t_emanuele(hyperbola_vect(u, f0, 0.0, 1.0, r_value)[1:2])[1]
	y(u, v) = v*t_emanuele(hyperbola_vect(u, f0, 0.0, 1.0, r_value)[1:2])[2]

	# colormap with transparency at the beginning
	cmap = :plasma
	n = 101
	g(x) = 0.5*(tanh(15*(x-0.05)) + 1)
	alphas = [(x == 0) ? 0.0 : g(x) for x in range(0, 1, length = n)]
	cmap_alpha = resample_cmap(cmap, n; alpha = alphas)

	prob_nz = zs[zs .> 0.0]
	n_levels = 3
	levels = 10.0.^range(log10.(0.002), log10.(maximum(prob_nz)), n_levels)

	f = Figure()
	ax = Axis(f[1, 1],
	    xlabel = L"n, \textrm{ initial photons}",
		ylabel = L"n + s, \textrm{ total photons in } \mathbf{k}",
		xgridvisible=false,
		ygridvisible=false,
		xlabelsize=28,
		ylabelsize=28,
		xticklabelsize=24,
		yticklabelsize=24,
		# title=L"|n,0\rangle",
		title=L"|n,n\rangle",
		titlesize=32,
		limits=(-0.5, nmax, -0.5, nmax),
		aspect=DataAspect(),
	)
	p1 = Makie.heatmap!(ax, xs, ys, zs;
		colormap=cmap_alpha,
		interpolate=false,
		colorscale=identity,
	)

	# For asymmetric state
	# Makie.band!(ax, 0:nmax,
	# 	map(x -> n_fock_k_gen(x, 0, v(r)) - sqrt(var_n_nm(x, 0, v(r))), 0:nmax), 
	# 	map(x -> n_fock_k_gen(x, 0, v(r)) + sqrt(var_n_nm(x, 0, v(r))), 0:nmax), 
	# 	color = (:aliceblue, 0.8))
	# Makie.lines!(ax, 0:nmax, map(x ->n_fock_k_gen(x, 0, v(r)), 0:nmax), color=:white, linestyle=:dash, linewidth=1.5)

	# For symmetric state
	plot_coordinate_mesh(x, y, ax; t_val=2.0, t_val_v=95.5, u_length=80, v_length=18, 
		u_lim = 10,
		v_lim = 95.5,
		t_length=200, 
		color=:royalblue4, 
		linewidth=1,
		linestyle=:dash
	)
	Makie.lines!(ax, 0:50, 0:50, color=:white, linestyle=:dash, linewidth=1.5)
	Makie.band!(ax, 0:nmax,
		map(x -> n_fock_k_gen(x, x, v(r)) - sqrt(var_n_nm(x, x, v(r))), 0:nmax), 
		map(x -> n_fock_k_gen(x, x, v(r)) + sqrt(var_n_nm(x, x, v(r))), 0:nmax), 
		color = (:aliceblue, 0.8))
	Makie.lines!(ax, 0:nmax, map(x ->n_fock_k_gen(x, x, v(r)), 0:nmax), color=:red, linestyle=:solid, linewidth=1.5)

	
	Colorbar(f[1, 2], p1; 
		label = L"p_s(n, 0)", 
		height = Relative(0.5), 
		ticklabelsize=24,
		labelsize=28,
	)
	# save(PATH_FIGS*"Figs_fock/Fig_1_1_new_var.pdf", f)
	f
end

# ╔═╡ 0f469be7-1653-49e4-bbf8-cb0f228d1a91
md"""
$slider_u

$slider_v

$slider_t

$slider_th

$slider_l

$slider_r
"""

# ╔═╡ 5b365b15-f017-4e53-985a-cb9ec962b125
let

	CairoMakie.activate!()
	
	u_length = 30
	v_length = 10
	t_length = 200
	# Initial setup
	f0 = SVector(0.0, 0.0)
	x(u, v) = v*t_emanuele(hyperbola_vect(u, f0, 0.0, 1.0, r_value)[1:2])[1]
	y(u, v) = v*t_emanuele(hyperbola_vect(u, f0, 0.0, 1.0, r_value)[1:2])[2]
	# x(u, v) = u^2 - v^2
	# y(u, v) = u * v + cos(u * v)
	u_range = range(-t_val, t_val, length=u_length)
	v_range = range(-t_val, t_val, length=v_length)
	
	u_val = minimum(u_range) + (maximum(u_range) - minimum(u_range)) * u_percent
	v_val = minimum(v_range) + (maximum(v_range) - minimum(v_range)) * v_percent
	t_min = -t_val
	t_max = t_val
	t_range = range(t_min, t_max, t_length)
	
	f = Figure()
	ax = Axis(f[1, 1],
		# xlabel = L"n, \textrm{ initial photons}",
		# ylabel = L"n + s, \textrm{ total photons}",
		# xlabelsize=20,
		# ylabelsize=20,
		# limits=(0, nmax, 0, 2*nmax),
		aspect=DataAspect(),
	)
	ax2 = Axis(f[1,2],
		aspect=DataAspect(),
		limits=(-limit, limit, -limit, limit),
	)
		
	# Plot in uv-plane
	uv_plot = lines!(ax, [(u_range[1], v) for v in t_min:t_max], color = :red, label = "u lines")
	for i in u_range[2:end]
		lines!(ax, [(i, v) for v in t_min:t_max], color = :red, label = "u lines")
	end
	for i in v_range
	    lines!(ax, [(u, i) for u in t_min:t_max], color = :blue, label = "v lines")
	end
	
	lines!(ax, [(u, v_val) for u in t_min:t_max], color = :blue, linewidth = thickness)
	lines!(ax, [(u_val, v) for v in t_min:t_max], color = :red, linewidth = thickness)
	
	
	# Transform to xy-plane
	T(u, v) = (x(u, v), y(u, v))
	x_vals(u, v) = x(u, v)
	y_vals(u, v) = y(u, v)
	
	xy_plot = lines!(ax2, [(x_vals(u_range[1], v), y_vals(u_range[1], v)) for v in t_range], color = :red, label = "u lines")
	for i in u_range[2:end]
		lines!(ax2, [(x_vals(i, v), y_vals(i, v)) for v in t_range], color = :red, label = "u lines")
	end
	for i in v_range
	    lines!(ax2, [(x_vals(u, i), y_vals(u, i)) for u in t_range], color = :blue, label = "v lines")
	end
	
	lines!(ax2, [(x_vals(u, v_val), y_vals(u, v_val)) for u in t_range], color = :red, linewidth = thickness)
	lines!(ax2, [(x_vals(u_val, v), y_vals(u_val, v)) for v in t_range], color = :blue, linewidth = thickness)
	
	f

end

# ╔═╡ 748b46be-2955-41b8-aedb-7b0383d07c9b
let
	f = Figure()
	ax = Axis(f[1,1], 
		limits=(0, 50, 0, 100)
	)
	
	# mesh functions
	r = 0.02
	f0 = SVector(0.0, 0.0)
	x(u, v) = v*t_emanuele(hyperbola_vect(u, f0, 0.0, 1.0, r)[1:2])[1]
	y(u, v) = v*t_emanuele(hyperbola_vect(u, f0, 0.0, 1.0, r)[1:2])[2]
	plot_coordinate_mesh(x, y, ax; t_val=3, t_val_v=140, u_length=30, v_length=10, 
		u_lim = 10,
		v_lim = 100,
		t_length=200, color=:black, linestyle=:dash
	)
	f
end

# ╔═╡ b37d5f91-98de-4acf-a645-a62aeb4754f7
let
	CairoMakie.activate!()

	nmax = 100
	r = 0.1
	x = 0:nmax
	y = -nmax:1:nmax
	m_range = [0, 1, 2, 5, 10, 50]
	z = [vec([(n, n + s, (s >= -min(m, n)) ? p_fock(c_l_hyper(n, m, s, u(r), v(r))) : NaN) for n = x, s = y]) for m in m_range]
	# z = vec([(n, n + s, (s >= -min(n,50)) ? p_fock(c_l_hyper(n, 50, s, u(r), v(r))) : NaN) for n = x, s = y])

	# colormap with transparency at the beginning
	cmap = :plasma
	n = 101
	g(x) = 0.5*(tanh(15*(x - 0.05)) + 1)
	alphas = [(x == 0) ? 0.0 : g(x) for x in range(0, 1, length = n)]
	cmap_alpha = resample_cmap(cmap, n; alpha = alphas)


	f = Figure(size = (1100, 600))
	num_col = length(m_range) ÷ 2
	for i in eachindex(m_range)
		f_idx = ((i - 1) ÷ num_col + 1, mod(i - 1, num_col) + 1)
		m = m_range[i]
		xs = [j[1] for j in z[i]]
		ys = [j[2] for j in z[i]]
		zs = [j[3] for j in z[i]]
		ax = Axis(f[f_idx...],
		    xlabel = (f_idx[1] == 2) ? L"n" : "",
		    ylabel = (f_idx[2] == 1) ? L"n + s" : "",
			title = L"|n, %$m \rangle",
			titlesize=22,
			xticklabelsvisible=(f_idx[1] == 2) ? true : false,
			yticklabelsvisible=(f_idx[2] == 1) ? true : false,
			xgridvisible=false,
			ygridvisible=false,
			xlabelsize=20,
			ylabelsize=20,
			limits=(0, nmax, 0, nmax),
			aspect=DataAspect(),
		)
		p1 = Makie.heatmap!(ax, xs, ys, zs;
			colormap=cmap_alpha,
			interpolate=false,
			colorscale=identity,
		)
		Colorbar(f[f_idx[1], f_idx[2]+3], p1; 
			label = (f_idx[2]+3 == 6) ? L"p_s(n, n)" : "", 
			height = Relative(0.8), 
			labelsize=20,
			ticklabelsize=10,
		)
	end
	
	
	# save(PATH_FIGS*"Figs_fock/Fig_4_cb_r.pdf", f)
	f
end

# ╔═╡ e8c18003-6717-450c-8e0c-78f09b939717
let
	CairoMakie.activate!()

	nmax = 100
	r = 0.1
	x = 0:nmax
	y = -nmax:1:nmax
	m_range = [0, 1, 2, 50]
	z = [vec([(n, n + s, (s >= -min(m, n)) ? p_fock(c_l_hyper(n, m, s, u(r), v(r))) : NaN) for n = x, s = y]) for m in m_range]
	# append!(z, [vec([(n, n + s, (s >= -min(n, n)) ? p_fock(c_l_hyper(n, n, s, u(r), v(r))) : NaN) for n = x, s = y])])
	# append!(m_range, [10])

	# colormap with transparency at the beginning
	cmap = :plasma
	n = 101
	g(x) = 0.5*(tanh(15*(x - 0.05)) + 1)
	alphas = [(x == 0) ? 0.0 : g(x) for x in range(0, 1, length = n)]
	cmap_alpha = resample_cmap(cmap, n; alpha = alphas)


	f = Figure(size = (800, 350))
	num_col = 4
	for i in eachindex(m_range)
		f_idx = ((i - 1) ÷ num_col + 1, mod(i - 1, num_col) + 1)
		m = m_range[i]
		xs = [j[1] for j in z[i]]
		ys = [j[2] for j in z[i]]
		zs = [j[3] for j in z[i]]
		ax = Axis(f[f_idx...],
		    xlabel = (f_idx[1] == 1) ? L"n, \textrm{ initial}" : "",
		    ylabel = (f_idx[2] == 1) ? L"n + s, \textrm{ total in } \mathbf{k}" : "",
			# title = (i < length(m_range)) ? L"|n, %$m \rangle" : L"|n, n \rangle",
			title = L"|n, %$m \rangle",
			titlesize=22,
			xticklabelsvisible=(f_idx[1] == 1) ? true : false,
			yticklabelsvisible=(f_idx[2] == 1) ? true : false,
			xgridvisible=false,
			ygridvisible=false,
			xlabelsize=20,
			ylabelsize=20,
			limits=(0, nmax, 0, nmax),
			aspect=DataAspect(),
		)
		p1 = Makie.heatmap!(ax, xs, ys, zs;
			colormap=cmap_alpha,
			interpolate=false,
			colorscale=identity,
		)
		Colorbar(f[f_idx[1]+1, f_idx[2]], p1; 
			# label = (i < length(m_range)) ? L"p_s(n, %$m)" : L"p_s(n, n)", 
			label = L"p_s(n, %$m)",
			width = Relative(1), 
			vertical = false,
			labelsize=20,
			ticklabelsize=10,
		)
	end
	rowsize!(f.layout, 1, Auto(0.5))
	
	# save(PATH_FIGS*"Figs_fock/Fig_4_new.pdf", f)
	f
end

# ╔═╡ 35338683-6dab-44dc-9787-121ce2c183cd
let
	CairoMakie.activate!()

	nmax = 100
	r = 0.1
	x = 0:nmax
	y = -nmax:1:nmax
	m_range = [0, 1, 2, 50]
	z = [vec([(n, n + s, (s >= -min(m, n)) ? p_fock(c_l_hyper(n, m, s, u(r), v(r))) : NaN) for n = x, s = y]) for m in m_range]
	# append!(z, [vec([(n, n + s, (s >= -min(n, n)) ? p_fock(c_l_hyper(n, n, s, u(r), v(r))) : NaN) for n = x, s = y])])
	# append!(m_range, [10])

	# colormap with transparency at the beginning
	cmap = :plasma
	n = 101
	g(x) = 0.5*(tanh(15*(x - 0.05)) + 1)
	alphas = [(x == 0) ? 0.0 : g(x) for x in range(0, 1, length = n)]
	cmap_alpha = resample_cmap(cmap, n; alpha = alphas)


	f = Figure(size = (400, 450))
	num_col = 2
	for i in eachindex(m_range)
		f_idx = ((i - 1) ÷ num_col + 1, mod(i - 1, num_col) + 1)
		m = m_range[i]
		xs = [j[1] for j in z[i]]
		ys = [j[2] for j in z[i]]
		zs = [j[3] for j in z[i]]
		ax = Axis(f[f_idx...],
		    xlabel = (f_idx[1] == 2) ? L"n, \textrm{ initial}" : "",
		    ylabel = (f_idx[2] == 1) ? L"n + s, \textrm{ total in } \mathbf{k}" : "",
			# title = (i < length(m_range)) ? L"|n, %$m \rangle" : L"|n, n \rangle",
			title = L"|n, %$m \rangle",
			titlesize=22,
			xticklabelsvisible=(f_idx[1] == 2) ? true : false,
			yticklabelsvisible=(f_idx[2] == 1) ? true : false,
			xgridvisible=false,
			ygridvisible=false,
			xlabelsize=20,
			ylabelsize=20,
			limits=(0, nmax, 0, nmax),
			aspect=DataAspect(),
		)
		p1 = Makie.heatmap!(ax, xs, ys, zs;
			colormap=cmap_alpha,
			interpolate=false,
			colorscale=identity,
		)
	end
	# rowsize!(f.layout, 1, Auto(0.5))
	
	# save(PATH_FIGS*"Figs_fock/Fig_2_new.pdf", f)
	f
end

# ╔═╡ d4159cb2-8269-4ea4-aadb-29c426a10e85
let
	CairoMakie.activate!()

	nmax = 100
	r = 0.01
	x = 0:nmax
	y = -nmax:1:nmax
	n_plus_s = 100
	z = vec([(n_s, n_plus_s - n_s, (n_plus_s - 2*n_s >= -min(n_s, n_s)) ? p_fock(c_l_hyper(n_s, n_s, n_plus_s - 2*n_s, u(r), v(r))) : NaN) for n_s = 0:n_plus_s])

	# colormap with transparency at the beginning
	cmap = :plasma
	n = 101
	g(x) = 0.5*(tanh(15*(x - 0.05)) + 1)
	alphas = [(x == 0) ? 0.0 : g(x) for x in range(0, 1, length = n)]
	cmap_alpha = resample_cmap(cmap, n; alpha = alphas)


	f = Figure()
	xs = [j[1] for j in z]
	ys = [j[3] for j in z]
	ax = Axis(f[1,1],
		xlabel = L"n_k, n_{-k}",
		ylabel = L"\textrm{Probability}",
		title = L"|n, n \rangle",
		titlesize=22,
		# xgridvisible=false,
		# ygridvisible=false,
		xlabelsize=20,
		ylabelsize=20,
		# limits=(0, nmax, 0, nmax),
		# aspect=DataAspect(),
	)
	p1 = Makie.lines!(ax, xs, ys;
	
	)
	
	# save(PATH_FIGS*"Figs_fock/Fig_4_new.pdf", f)
	f
end

# ╔═╡ 3a5ac21b-67eb-4aa8-959f-6c080e5219e5
let
	CairoMakie.activate!()

	nmax = 40
	n = 0
	r_range = 1.0 .+ 10.0.^range(-2, 2, 5*nmax)
	# r_range = 10.0.^range(-0.05, 0.01, 5*nmax)
	T = map(r -> temp_vac(r), r_range)
	x = r_range
	y = -n:1:nmax
	z = vec([(temp_vac(r), n + s, (s >= -n) ? p_fock(c_l_hyper(n, n, s, u(r), v(r))) : NaN) for r = x, s = y])

	xs = [i[1] for i in z]
	ys = [i[2] for i in z]
	zs = [i[3] for i in z]

	# colormap with transparency at the beginning
	cmap = :plasma
	n_grad = 101
	g(x) = 0.5*(tanh(15*(x-0.05)) + 1)
	alphas = [(x == 0) ? 0.0 : g(x) for x in range(0, 1, length = n_grad)]
	cmap_alpha = resample_cmap(cmap, n_grad; alpha = alphas)

	f = Figure()
	ax = Axis(f[1, 1],
	    xlabel = L"k_B T / \hbar \omega_j",
	    ylabel = L"n + s, \textrm{ total photons}",
		xgridvisible=false,
		ygridvisible=false,
		xlabelsize=20,
		ylabelsize=20,
		xscale=log10,
		limits=(T[1], T[end], -0.5, nmax+n+0.5),
		aspect=1,
	)
	p1 = Makie.heatmap!(ax, xs, ys, zs;
		colormap=cmap_alpha,
		interpolate=false,
		colorscale=identity,
	)
	Makie.band!(ax, T, zeros(length(T)), map(r -> n_fock_sym(0, v(r)) + 3*sqrt(var_n_nm(0,0,v(r))), r_range), color = (:red, 0.05))
	Makie.lines!(ax, T, map(r -> n_fock_sym(0, v(r)), r_range), color=:black, linestyle=:dash, linewidth=1.5)
	Makie.lines!(ax, T, map(r -> n_fock_sym(0, v(r)) + 3*sqrt(var_n_nm(0,0,v(r))), r_range), color=(:red, 0.5), linewidth=0.5)

	text!(800, 30, text = L"p_s(0,0)", 
		align = (:center, :center),
		color=:firebrick4,
		fontsize = 20,
	)
	text!(400, 14, text = L"n_\mathrm{vac}", 
		align = (:center, :center),
		color=:black,
		fontsize = 20,
	)
	Colorbar(f[1, 2], p1; 
		label = L"p_s(n, n)", 
		height = Relative(0.5), 
		labelsize=24,
	)
	# save(PATH_FIGS*"Figs_fock/Fig_4_2_rev.pdf", f)
	f
end

# ╔═╡ 953b86c7-65a5-4bf0-93fb-2003d5b3d555
let
 	f = Figure()
	ax = Axis(f[(1, 1)...])
	x(u,v) = u^2 + v^2
	y(u, v) = u*v
	plot_coordinate_mesh(x, y, ax)
	f
end

# ╔═╡ d6ec3e69-e2e0-471b-aa84-e16019efbee6
md"""
## Statistics
"""

# ╔═╡ 5978f543-f872-441f-8f7b-71a4f110105b
let
	CairoMakie.activate!()

	r = 10.0 .^ range(-3.0, 0.0, 200)
	n_vac = map(x -> v(x)^2, r)
	
	f = Figure()
	ax = Axis(f[1, 1],
		xlabel = L"r = \omega_1 / \omega_0",
		ylabel = L"g^{(2)}_k(0)",
		xlabelsize=20,
		ylabelsize=20,
		xscale=log10,
		# yscale=log10,
		# limits=(0.009, 1.1, -0.1, 8.0),
		# aspect=DataAspect(),
	)

	for n in 0:100
		lines!(ax, n_vac, map(x -> g2_n0(n, v(x)), r),
			color = :blue, 
			linewidth = 2,
			alpha=0.2,
			label=(n == 0) ? L"|n,0\rangle_0" : "",
		)
		lines!(ax, n_vac, map(x -> g2_sym(n, v(x)), r),
			color = :red, 
			linewidth = 2,
			alpha=0.2,
			label=(n == 0) ? L"|n,n\rangle_0" : "",
		)
		# if n > 0
		# 	lines!(ax, r, map(x -> g2_sym_cross(n, v(x)), r),
		# 		color = :green, 
		# 		linewidth = 2,
		# 		alpha=0.2,
		# 	)
		# end
		if n == 0
			Legend(f[1, 1], ax, "";
				nbanks = 1, 
				orientation=:vertical,
				tellheight = false,
        		tellwidth = false,
				halign = :right,
				valign = :bottom,
				margin = (10, 10, 10, 10),
			)
		end
	end
	Makie.hlines!(1.5; color=:red, linestyle=:dash, linewidth=1)
	Makie.hlines!(1.0; color=:blue, linestyle=:dash, linewidth=1)
	Makie.text!(30, 1.2;
		text=L"\downarrow \; n",
		fontsize=28,
		color=:blue,
	)
	Makie.text!(0.7, 0.3;
		text=L"\uparrow \; n",
		fontsize=28,
		color=:red,
	)
	Makie.text!(0.01, 1.85;
		text=L"n = 0 \textrm{ (vacuum)}",
		fontsize=22,
		color=:plum4,
	)
	Makie.scatter!(1, 2; color=:white, strokecolor=:black, strokewidth=1, marker=:circle, markersize=10)

	# save(PATH_FIGS*"Figs_fock/Fig_4_a.pdf", f)
	f
end

# ╔═╡ ec51735a-38f8-4bf8-8954-91353dbe8b9f
g2_k(10, 0, u(0.1), v(0.1), 1000)

# ╔═╡ 390f2ad3-1203-4c0f-81d7-55bf1ded25d7
# ╠═╡ disabled = true
#=╠═╡
let
	CairoMakie.activate!()

	n_1 = 6
	n_2 = 1
	r = 10.0 .^ range(-1.0, 1.0, 100)
	g2_gen = map(x -> g2_k(n_1, n_2, u(x), v(x), 1000), r)
	
	f = Figure()
	ax = Axis(f[1, 1],
		xlabel = L"r = \omega_1 / \omega_0",
		ylabel = L"g^{(2)}(0)",
		# xlabelsize=20,
		# ylabelsize=20,
		xscale=log10,
		# yscale=log10,
		# limits=(0.009, 1.1, -0.1, 8.0),
		# aspect=DataAspect(),
	)

	lines!(ax, r, g2_gen,
		color = :black, 
		linewidth = 2,
		# alpha=0.2,
	)
	lines!(ax, r, map(x -> g2_nm(n_1, n_2, v(x)), r),
		color = :red, 
		linewidth = 2,
		# alpha=0.2,
	)
	f
end
  ╠═╡ =#

# ╔═╡ 0292ac23-70e3-416c-b0a6-31e74bf459d7
# ╠═╡ disabled = true
#=╠═╡
let
	CairoMakie.activate!()

	n_1 = 10
	n_2 = 4
	r = 10.0 .^ range(-1.0, 1.0, 100)
	n_gen = map(x -> n_fock_k(n_1, n_2, u(x), v(x), 1000), r)
	
	f = Figure()
	ax = Axis(f[1, 1],
		xlabel = L"r = \omega_1 / \omega_0",
		ylabel = L"\langle n \rangle",
		# xlabelsize=20,
		# ylabelsize=20,
		xscale=log10,
		# yscale=log10,
		# limits=(0.009, 1.1, -0.1, 8.0),
		# aspect=DataAspect(),
	)

	lines!(ax, r, n_gen,
		color = :black, 
		linewidth = 2,
		# alpha=0.2,
	)
	lines!(ax, r, map(x -> n_fock_k_gen(n_1, n_2, v(x)), r),#map(x -> n_fock_n0(n_1, v(x)) + 0.0*n_fock_n0(n_2, v(x)), r),
		color = :red, 
		linewidth = 2,
		# alpha=0.2,
	)
	f
end
  ╠═╡ =#

# ╔═╡ ace0e2f7-cbf1-4767-a2cc-5f69b1a47f19
let
	CairoMakie.activate!()

	nmax = 100
	r = 0.02
	x = y = 0:nmax
	z = vec([(n, m, g2_nm(n, m, v(r))) for n = x, m = y])

	xs = [i[1] for i in z]
	ys = [i[2] for i in z]
	zs = [i[3] for i in z]

	# colormap with transparency at the beginning
	cmap = :plasma
	n = 101
	g(x) = 0.5*(tanh(16*(x-0.1)) + 1)
	alphas = [(x == 0) ? 0.0 : g(x) for x in range(0, 1, length = n)]
	cmap_alpha = resample_cmap(cmap, n; alpha = alphas)

	f = Figure()
	ax = Axis(f[1, 1],
	    xlabel = L"n_k",
	    ylabel = L"n_{-k}",
		xgridvisible=false,
		ygridvisible=false,
		xlabelsize=20,
		ylabelsize=20,
		xscale=Makie.pseudolog10,
		yscale=Makie.pseudolog10,
		# limits=(0, nmax, 0, nmax),
		aspect=DataAspect(),
	)
	p1 = Makie.heatmap!(ax, xs, ys, zs;
		colormap=cmap_alpha,
		interpolate=false,
		colorscale=identity,
	)
	Makie.contour!(ax, xs, ys, zs; levels=[2.0, 1.498, 1.4, 1.3, 1.2, 1.1], labels=true, color=:azure)
	Colorbar(f[1, 2], p1; 
		label = L"g^{(2)}_k(0)", 
		height = Relative(0.5), 
		labelsize=24,
	)

	# save(PATH_FIGS*"Figs_fock/Fig_4_c_log.pdf", f)
	f
end

# ╔═╡ ff13b396-afdb-4155-82f9-88330ca51fd1
let
	CairoMakie.activate!()

	# r = 10.0 .^ range(-3.5, 3.5, 1000)
	r = 10.0 .^ range(0.02, 3.0, 1000)
	n_vac = map(x -> v(x)^2, r)
	n = 8
	m_range = 0:2:20
	cmap = cgrad(:Cassatt1, length(m_range), categorical = true)
	
	f = Figure()
	ax = Axis(f[1, 1],
		xlabel = L"r = \omega_1 / \omega_0",
		ylabel = L"g^{(2)}_k(0)",
		xlabelsize=20,
		ylabelsize=20,
		xscale=log10,
		# yscale=log10,
		# limits=(0.009, 1.1, -0.1, 8.0),
		# aspect=DataAspect(),
	)

	idx = 0
	for m in m_range
		idx += 1
		lines!(ax, n_vac, map(x -> g2_nm(n, m, v(x)), r),
			color = cmap[idx], 
			linewidth = 2,
			label=L"|%$n, %$m \rangle_0",
			# alpha=0.2,
		)
	end
	Legend(f[1, 2], ax, L"|n,m\rangle_0"; 
		nbanks = 1,
		orientation=:vertical,
	)
	# save(PATH_FIGS*"Figs_fock/Fig_4_b.pdf", f)

	f
end

# ╔═╡ c34aca81-de76-4b6b-b389-5f8f73129c29
let
	r = 10.0 .^ range(0.0, 3.0, 1000)
	n_vac = map(x -> v(x)^2, r)
	g2 = map(x -> g2_nm(100, 0, v(x)), r)
	lines(g2)
end

# ╔═╡ Cell order:
# ╠═f19b8814-28f3-11ef-37ed-cb82cc7cdffe
# ╠═30e8847a-c88e-4f68-8e19-4a2044569ae3
# ╠═134acf44-bc13-4744-916b-52b1339d5936
# ╟─53ed7c53-2395-4a5b-b479-14e7ae7181cb
# ╠═0d1115cd-a481-4bff-a5ee-68a9cf4dde42
# ╟─d78dbee3-f040-4753-82e9-39362d4204a1
# ╟─216474fc-a9fc-4629-b073-529c926e8e37
# ╟─19713970-9020-4b20-9e1e-734dbfeb37e5
# ╟─156068dc-baa5-448a-9fe9-44d058650968
# ╟─a2e10c75-44ff-442e-b9aa-30996ee25be7
# ╟─d6a16977-c42e-417f-8829-8fb17e200422
# ╟─ecae7206-cc04-4a03-89c4-9c21672ad0b1
# ╠═2da5b512-8cb1-4fae-a4af-371361ccb0f6
# ╠═a900ceb7-5e4f-4ddf-a7fc-a0d1e7f52819
# ╠═ed403f68-031e-4164-af23-23f8b4632691
# ╠═05bf53c8-22c5-43ea-90cd-338f4cb230e7
# ╠═be504ed6-e88b-4150-a7ba-64788021d5b7
# ╠═d5e9d066-8aac-4993-a8c1-b990998be10a
# ╟─0f469be7-1653-49e4-bbf8-cb0f228d1a91
# ╠═5b365b15-f017-4e53-985a-cb9ec962b125
# ╠═748b46be-2955-41b8-aedb-7b0383d07c9b
# ╠═b37d5f91-98de-4acf-a645-a62aeb4754f7
# ╠═e8c18003-6717-450c-8e0c-78f09b939717
# ╠═35338683-6dab-44dc-9787-121ce2c183cd
# ╠═d4159cb2-8269-4ea4-aadb-29c426a10e85
# ╠═3a5ac21b-67eb-4aa8-959f-6c080e5219e5
# ╠═953b86c7-65a5-4bf0-93fb-2003d5b3d555
# ╟─d6ec3e69-e2e0-471b-aa84-e16019efbee6
# ╠═5978f543-f872-441f-8f7b-71a4f110105b
# ╠═ec51735a-38f8-4bf8-8954-91353dbe8b9f
# ╟─390f2ad3-1203-4c0f-81d7-55bf1ded25d7
# ╟─0292ac23-70e3-416c-b0a6-31e74bf459d7
# ╠═ace0e2f7-cbf1-4767-a2cc-5f69b1a47f19
# ╠═ff13b396-afdb-4155-82f9-88330ca51fd1
# ╠═c34aca81-de76-4b6b-b389-5f8f73129c29
