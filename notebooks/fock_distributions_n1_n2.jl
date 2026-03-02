### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 5186070a-c993-41f9-a34a-d18c8fd06e28
begin
	using Plots
	using BenchmarkTools
	using PlutoUI
	using HDF5, FileIO, Printf, LaTeXStrings
	using Nemo
	using BigCombinatorics
	using Revise

	NMAX = 50
end

# ╔═╡ e2a6ec34-defa-11ed-36e9-01e8a5006b3a
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
    using .TI_quantum

	# initialize the table for factorials
	Factorial_1()

    PATH_FIGS, PATH_DATA = path()
end

# ╔═╡ b6ce4930-5d86-421e-bfcd-09addd6f5f9f
md"""
# Compute $P_{fock}(n_1, n_2)$
"""

# ╔═╡ 138a0225-abb6-4b2f-accd-8b53563c9711
begin
	r_list = range(0.001, 0.999, 1000)
end

# ╔═╡ dab45c67-5af9-405b-aeb4-4d349a0bb5e9
# ╠═╡ disabled = true
#=╠═╡
sum([abs.(real(convert(ComplexF64,c_l_hyper(50, 50, l, 0.001)))).^2 for l=-50:100000])
  ╠═╡ =#

# ╔═╡ 38a10c4e-9f03-4ffc-bb19-c51f8cf13371
@bind r_slider PlutoUI.Slider(r_list, default=r_list[end])

# ╔═╡ e6ba2d26-a8aa-46fd-b31f-cd3f45987875
let
	n = [string(i) for i = 0:NMAX-1]
	l = [string(i) for i = -NMAX-1:NMAX]
	prob = zeros(length(n), length(l))
	for i in eachindex(n), j in eachindex(l)
		n_i = i - 1
		l_j = (-NMAX-1:NMAX)[j]
		if l_j < -n_i
			prob[i, j] = 0.0
		else
			prob[i, j ]	= abs.(real(convert(ComplexF64,c_l_hyper(n_i, n_i, l_j, r_slider)))).^2		
		end
	end
	#prob = [abs.(real(convert(ComplexF64,c_l_hyper(n_i, n_i, l_j, r_slider)))).^2 for n_i=0:NMAX-1, l_j=0:NMAX-1]
	heatmap(n, l, (prob.^(1/3))', aspect_ratio = 1/3, 
		ylabel = L"$l\;$(number of pairs)",
		xlabel = L"n_1, n_2 = n",
	)
	annotate!([(NMAX ÷ 5 * 1, NMAX ÷ 10 * 9, (L"r = "*string(round(r_slider, digits=3)), 12, :blue, :center))])
end

# ╔═╡ Cell order:
# ╟─b6ce4930-5d86-421e-bfcd-09addd6f5f9f
# ╠═e2a6ec34-defa-11ed-36e9-01e8a5006b3a
# ╠═5186070a-c993-41f9-a34a-d18c8fd06e28
# ╠═138a0225-abb6-4b2f-accd-8b53563c9711
# ╠═dab45c67-5af9-405b-aeb4-4d349a0bb5e9
# ╠═38a10c4e-9f03-4ffc-bb19-c51f8cf13371
# ╠═e6ba2d26-a8aa-46fd-b31f-cd3f45987875
