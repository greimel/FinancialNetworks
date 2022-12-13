### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# ╔═╡ 7d711677-35dd-46d2-ba87-f2994d483dc3
using NetworksUtils: FinancialNetwork, CompleteNetwork, figure, γNetwork, IslandNetwork, InterbankMarket, is_regular, label, numbered_graphplot, numbered_graphplot!, nested_circles, componentwise_circle

# ╔═╡ 952825f3-7f66-4435-a850-cbe160106f6f
using DataFrames

# ╔═╡ c59a6cc9-f9bd-4831-9bd8-7bc331d5790f
using Chain: @chain

# ╔═╡ 63ed5af7-2c4d-4184-8559-ec6977f04709
using Roots

# ╔═╡ d5ac4a10-48f0-4e59-a186-eff9dd059717
using Statistics: mean

# ╔═╡ 162e1da7-0d14-464f-8f64-66804009bc5d
using Graphs

# ╔═╡ c46f9af5-6bd5-49fe-8086-44b238d5a8a1
using LinearAlgebra: I

# ╔═╡ fc4fa10d-2b39-4f10-bac8-7a1e5467669d
using Makie

# ╔═╡ fcb09a02-70ca-4da0-bab9-6f2d9da6a195
using SimpleWeightedGraphs

# ╔═╡ 3d037a99-7978-4fce-b603-535a90362f40
# ╠═╡ skip_as_script = true
#=╠═╡
using PlutoUI: PlutoUI, as_svg, Slider, TableOfContents
  ╠═╡ =#

# ╔═╡ 92c4efe9-6ba1-4895-a3c8-05ec312a3820
using Colors

# ╔═╡ d004095e-2910-4a7c-8ee4-b5d13c34419e
using PlutoUI: CheckBox

# ╔═╡ 7e12dd72-2d66-4d94-8af6-084e80a63885
using StructArrays

# ╔═╡ fde5439e-d1ed-43e7-be89-d8f41eff3817
# ╠═╡ skip_as_script = true
#=╠═╡
using AlgebraOfGraphics
  ╠═╡ =#

# ╔═╡ 1b470abb-bf1e-4ec6-a5c8-5f05e86c87dc
using DataFrameMacros

# ╔═╡ 57945142-4153-4c89-bef6-bb8bb0b6a7d8
using DocStringExtensions

# ╔═╡ ec9e91ef-fff8-4d7f-8df3-86df19ced7ed
md"""
`FinancialNetworks-nb.jl` | _last updated on December 13_
"""

# ╔═╡ 213aa0db-1aac-4471-b7eb-42f14d23a5c2
md"""
# Useful functions for Financial networks
"""

# ╔═╡ b90bc5c2-5d7c-4e0c-bf71-08b1af48e754
initial_network(interbank_market) = adjacency_matrix(interbank_market.network)

# ╔═╡ b6238421-76d2-4da0-a1e0-d0e94bed1b0d
updated_network(interbank_market) = interbank_market.payments

# ╔═╡ 2e88a6a7-79a5-4ef2-b3f4-ac916e09edf6
begin
	payables(IM) = sum(initial_network(IM), dims = 2)
	receivables(IM) = sum(initial_network(IM), dims = 1)
	
	paid(IM) = sum(updated_network(IM), dims = 2)
	received(IM) = sum(updated_network(IM), dims = 1)
end

# ╔═╡ 8eb2aa7c-2dc1-11ed-24ab-e9cce1f5013d
md"""
## Firms (Projects)
"""

# ╔═╡ 8fd6c40f-17b1-461b-82ec-100c6c987701
firms_doc = raw"""
Firms start with a working capital of ``1``. The working capital is needed for production. If the firm stays fully capitalized it will pay out dividends ``z_i = a - \varepsilon_i`` in period ``t = 1`` and ``A`` in period ``t=2``. The payments are 
```math
(z_i, A).
```

Owner ``i \in \{1, \ldots, N\}`` of the firm can withdraw (liquidate) a fraction ``\ell_i \in [0, 1]`` of their share ``\omega_i`` of the working capital of the firm. The total fraction liquidated is
```math
\bar \ell = \sum_{i=1}^n = \omega_i \ell_i.
```

In that case the firm will withold part of the dividends to replace the capital. The firms payoff will be
```math
(z_i - \alpha \bar\ell, A)
```
Each individual owner ``i`` will have some costs and revenues from liquidating the fraction ``\ell_i`` if their share.
```math
(\alpha \ell_i + \ell_i \zeta A, - \ell_i A)
```
""";

# ╔═╡ 3c867227-aa93-4ade-9434-a1a61904656a
Markdown.parse(firms_doc)

# ╔═╡ eae2ef96-deeb-4e17-85fe-00ef53d13a46
"""
$firms_doc

$(FIELDS)
"""
Base.@kwdef struct Firm{T}
	"Payoff in ``t=2``"
	A::T = 3.0
	"Payoff in ``t=1`` (absent shock; think _dividend_)"
	a::T = 1.0
	"share of capital that needs to be replaced in case of liquidation - check relationship to _alienability_ of capital"
	α = 0.5
	"fraction of non-replacable capital that can be recovered in case of liquidation"
	ζ::T = 0.1
	"fraction of replacable capital that can be recovered in case of liquidation"
	ζ_α::T = 1.0
	#function Firm(A, a, α, ζ, ζ_α)		
	#	new{Float64}(float(A), float(a), float(α), float(ζ), float(ζ_α))
	#end
end

# ╔═╡ 3d741a24-60c5-4462-b048-796efe3ef0e4
"""
### Dividend payment in ``t = 1 ``

The baseline dividend ``a`` is reduced if there is a shock to the project (``ε > 0``) or if the firm is partially liquidated (``\\bar ℓ > 0``). In the latter case the firm will first use their "profit" to replace the liquidated replacable capital ``\\bar ℓ α`` and pay out the rest as dividends. That is,

```math
\\max\\{a - ε - \\bar ℓ ⋅ α, 0 \\}
```

In the model of Acemoglu, Ozdaglar and Tahbaz-Salehi (``α = 0``) this reduces to ``\\max\\{a - ε, 0\\}``. 

Dividends are paid out by share ``ω_i``.
"""
dividend((; a, α), ε, ℓ̄) = max(a - ε - ℓ̄ * α, 0.0)

# ╔═╡ 24fbe0a5-f5dd-49b3-a667-963b25168191
"""
### Proceeds from liquidation (in ``t=1``)

Use the _internal_ function `_recovery_` which computes the proceeds from liquidating
* replacable capital (machines; ``ℓ ⋅ ζ_α ⋅ α``) and
* non-replacable capital (field; ``ℓ ⋅ ζ ⋅ A``)
* total (``ℓ ⋅ (ζ_α ⋅ α + ζ ⋅ A)``)

In the model of Acemoglu, Ozdaglar and Tahbaz-Salehi (``α = 0``) this reduces to (``ℓ ⋅ ζ ⋅ A``).
"""
function recovery((; ζ, A, α, ζ_α), ℓ)
	rec_machines = ℓ * ζ_α * α
	rec_field    = ℓ * ζ   * A
	rec = rec_machines + rec_field
	(; rec, rec_machines, rec_field)
end

# ╔═╡ bfed935f-0b08-491a-a29a-c118c35134d6
"""
### Cash-flow from firm to bank in ``t=1``

Sum of `dividend` and `recovery`.

"""
function cashflow₁(firm, ε, ℓ, ℓ̄)
	div = dividend(firm, ε, ℓ̄)
	(; rec, rec_machines, rec_field) = recovery(firm, ℓ)
	cf = div + rec
	(; cf, div, rec, rec_machines, rec_field)
end

# ╔═╡ f5f971e9-bb3c-44f2-b269-bdaf17bff943
let
	f = Firm()
	ε = 0.0
	ℓ = 0.0
	ℓ̄ = 0.1
	dividend(f, ε, ℓ̄)
	recovery(f, ℓ)
	cashflow₁(f, ε, ℓ, ℓ̄)
end

# ╔═╡ 8709bbc9-ee03-42ea-ac7c-e2b1c123eaed
"""
### Cash-flow from firm to bank in ``t=2``

**TODO** Handle the case where the replacable capital cannot be fully replaced

"""
cashflow₂((; A), ℓ) = (1 - ℓ) * A

# ╔═╡ 57247692-b967-4880-9574-af5a45fffc0b
md"""
## Banks
"""

# ╔═╡ 75410a0d-a352-435a-88d8-7b520bf2005f
"""
This is a bank

$(FIELDS)
"""
Base.@kwdef struct Bank
	"Outside (senior) obligations (liability)."
	ν = 4.5
	"Outside (senior) assets (cash)."
	c = 0.0
	"shares in firms"
	shares = [0.5]
end

# ╔═╡ 5e0ce818-4ac9-4efb-ad24-7223cac5cb12
md"""
## Payments and Equilibrium
"""

# ╔═╡ 2da46c3d-dd64-424f-8fe9-b6ba1156f423
function _repay_((; c, ν), (; x̄, ȳ), firm, (; ωᵢ, ε, ℓ̄₋ᵢ,))
	ℓ̄(ℓᵢ) = (1-ωᵢ) * ℓ̄₋ᵢ + ωᵢ*ℓᵢ
	assets(ℓᵢ) = x̄ + c + ωᵢ * cashflow₁(firm, ε, ℓᵢ, ℓ̄(ℓᵢ)).cf

	if ν + ȳ ≤ assets(0.0)
		out = (; ℓ=0, y_pc=1.0, ν_pc=1.0)
	elseif assets(0.0) < ν + ȳ ≤ assets(1.0)
		ℓᵢ = find_zero(ℓ -> ν + ȳ - assets(ℓ), (0.0, 1.0))
		out = (; ℓ=ℓᵢ,   y_pc=1.0, ν_pc=1.0)
	else
		ℓᵢ = 1.0
		ass = assets(ℓᵢ)
		if ν ≤ ass < ν + ȳ
			out = (; ℓ=ℓᵢ, y_pc = (ass - ν)/ȳ, ν_pc = 1.0)
		else # assets < ν
			out = (; ℓ=ℓᵢ, y_pc = 0.0, ν_pc = ass/ν)
		end
	end		

	(; out..., assets = assets(out.ℓ), cashflow₁(firm, ε, out.ℓ, ℓ̄(out.ℓ))..., firm_pc = 1 - out.ℓ, x̄, ȳ = ȳ * out.y_pc)
	
end

# ╔═╡ 62a50fdf-4c81-4adf-9366-40362cbefb83
function repay((; c, ν, shares), (; x̄, ȳ), firms, εs, ℓ̄₋ᵢ)
	i_firm = only(findall(shares .> 0))
	ωᵢ = shares[i_firm]
	firm = firms[i_firm]
	ε = εs[i_firm]

	_repay_((; c, ν), (; x̄, ȳ), firm, (; ωᵢ, ε, ℓ̄₋ᵢ,))
	
end

# ╔═╡ 1cad0937-5ecd-4b0a-871c-de8c2f75d7b7
function iterate_payments(banks, IM, firms, shares, εs, ℓs = zeros(length(banks)))
	x = updated_network(IM)
	y = initial_network(IM)
	
	x_new = copy(y)
	out = map(enumerate(banks)) do (i, bank)
		# compute repayment
		x̄ = received(IM)[i]
		ȳ = payables(IM)[i]

		i_firm = findall(shares[:,i] .> 0) |> only
		
		co_lenders = @chain begin
			shares[i_firm,:]
			findall(_ .> 0)
			filter(≠(i), _)
		end
		# FIXME! This assumes equal shares in firms
		ℓ₋ᵢ = length(co_lenders) > 0 ? mean(ℓs[co_lenders]) : 0.0
		rpy = repay(bank, (; x̄, ȳ), firms, εs, ℓ₋ᵢ)
		(; y_pc, ℓ) = rpy
	
		ℓs[i] = ℓ
		# update payment matrix
		x_new[i, :] .*= y_pc
		# return bank's choices
		(; bank = i, rpy...)
	end

	(; x_new, out, ℓs)
end

# ╔═╡ 58481389-ded5-4b52-9760-ece0d381451e
function equilibrium(banks, IM₀, firms, shares, εs; maxit = 100)
	IM = deepcopy(IM₀)
	x = updated_network(IM)
	y = initial_network(IM)

	ℓs = zeros(length(banks))
	
	x_new = copy(x)
	for it ∈ 1:maxit
		ℓs_old = copy(ℓs)
		(; x_new, out, ℓs) = iterate_payments(banks, IM, firms, shares, εs, copy(ℓs_old))
		converged = x_new ≈ x && ℓs_old ≈ ℓs
		x .= x_new
		
		if converged || it == maxit
			bank_df = DataFrame(out)
			
			ℓ̄ = shares * bank_df.ℓ
			pc_div((; a, α), ℓ) = (a - α * ℓ)/a
			δ_pc = pc_div.(firms, ℓ̄)
			
			firm_df = DataFrame(; ℓ̄, δ_pc)
			
			# interbank_df = DataFrame(; x, y)

			return (; bank_df, firm_df, it, success = it != maxit, banks, IM)
		end
	end
	
end

# ╔═╡ f49f7e4d-20bc-4dd6-9586-207dbfd19cef
md"""
## Setting up the environment
"""

# ╔═╡ 6b31b940-0a2e-4bbf-9534-bdfa5a02c3c4
n_banks = 6

# ╔═╡ afa36073-424e-4ff8-b584-554ae66cee06
md"""
### When will there be transmission?

1. ``a ≤ α`` so that there is a chance that dividends shrink towards zero
2. ``0 ≤ ν + y - c - x ≤ a``, so that there is default when dividends shrink to zero
"""

# ╔═╡ 51f5c242-54af-4511-9838-58326756197d
#=╠═╡
md"""
* ``a``: $(@bind a Slider(0:0.01:1.0, default = 0.6, show_value = true))
* ``ε``: $(@bind ε Slider(0:0.01:2.0, default = 0.0, show_value = true))
* ``γ``: $(@bind γ Slider(0:0.05:1.0, default = 0.5, show_value = true))
* ``α``: $(@bind α Slider(0:0.05:0.99, default = 0.5, show_value = true))
* ``ζ``: $(@bind ζ Slider(0:0.005:1.0, default = 0.0, show_value = true))
* ``ζ_α``: $(@bind ζ_α Slider(0:0.01:1.0, default = 1.0, show_value = true))
* ``\bar y``:  $(@bind ȳ Slider(0:0.005:0.99, default = 0.25, show_value = true))
* ``c``:  $(@bind c Slider(0:0.005:0.99, default = 0.25, show_value = true))
* ``ν``:  $(@bind ν Slider(0:0.005:0.99, default = 0.50, show_value = true))

"""
  ╠═╡ =#

# ╔═╡ eabcafbd-5c35-4ee9-a50c-d129f8d3a34a
#=╠═╡
begin
	n_islands = 2
	#ȳ = 0.25
	n = IslandNetwork(n_islands, n_banks ÷ n_islands, ȳ; γ)
	IM = InterbankMarket(n)
end;
  ╠═╡ =#

# ╔═╡ 0a2dfc44-67f8-4e9f-991f-906509116d66
#=╠═╡
RES = let
	shocked_bank = 2
	show_firms = true
	common_lenders = true
	
	if !common_lenders #AOTS
		n_firms = n_banks
		shares = I(n_banks)
	else
		n_firms = n_banks ÷ 2
		shares = zeros(n_firms, n_banks)
		for i_bank ∈ 1:n_banks
			i_firm = (i_bank+1) ÷ 2
			shares[i_firm, i_bank] = 0.5
		end
	end

	firms = [Firm(; α, ζ, ζ_α, a) for _ ∈ 1:n_firms]
	εs = zeros(n_banks); εs[shocked_bank] = ε
	
	banks = [Bank(; ν, c=c - (i==shocked_bank)*ε, shares=shares[:,i]) for i ∈ 1:n_banks]

	IMx = deepcopy(IM)
	
	(; bank_df, firm_df) = equilibrium(banks, IMx, firms, shares, εs)

	(; bank_df, firm_df, IM=IMx, show_firms, shares, firms, banks)
end
  ╠═╡ =#

# ╔═╡ 1d00ca98-985a-489c-9ef4-d74f558d9b71
#=╠═╡
RES.firm_df
  ╠═╡ =#

# ╔═╡ d5dc8002-0d8e-4681-a6a5-4d5fc49f9288
# ╠═╡ skip_as_script = true
#=╠═╡
import CairoMakie
  ╠═╡ =#

# ╔═╡ 72de5ce2-6bab-4f9a-b45a-4559fa3affea
function mix_with_white(color, α)
	(; r, g, b) = RGB(color)
	(r, g, b) .* (1-α) .+ α |> x -> RGB(x...)
end

# ╔═╡ d816a40a-67c1-4a1a-ae0c-d26502ba3918
mix_with_white(colorant"limegreen", 0.5)

# ╔═╡ cbcdef2e-d854-47c9-bfc4-595c90ea50fb
md"""
### Visualize bank-firm-network
"""

# ╔═╡ 8028be99-7559-4974-91a6-36ccab2d4a7f
function visualize_bank_firm_network!(ax, IM, shares, bank_df, firm_df; r = 1.4, start = Makie.automatic, layout=Makie.automatic, kwargs...)

	n = IM.network
	g₀ = SimpleWeightedDiGraph(n.Y)
	A = adjacency_matrix(g₀)

	big_A = [A       shares';
	         0 * shares 0 * I] 	

	g = big_A |> SimpleWeightedDiGraph

	n_banks = size(A, 1)
	n_firms = size(shares, 1)

	edge_attr_df = map(edges(g)) do (; src, dst, weight)
		if src ≤ n_banks && dst ≤ n_banks
			linewidth = weight
			linestyle = Makie.automatic
			arrow_size = 7
		elseif dst > n_banks
			linewidth = 0.5
			arrow_size = 0.0
			linestyle = "--"
		end
		(; linewidth, arrow_size, linestyle)
	end |> DataFrame

	arrow_attr = (; markersize = edge_attr_df.arrow_size)
	edge_attr = (; edge_attr_df.linewidth, edge_attr_df.linestyle)

	start = start === Makie.automatic ? 1/n_banks : start
	layout = layout === Makie.automatic ? nested_circles(n_banks; start, r) : layout
	numbered_graphplot!(ax, g;
		layout,
		#figure = figure(220),
		nlabels = [string.(1:n_banks); ["F$i" for i ∈ 1:n_firms]],
		node_marker = [fill(:circle, n_banks); fill(:rect, n_firms)],
		node_color = [ifelse.(bank_df.y_pc .< 1.0, :red, ifelse.(bank_df.ℓ .> 0.0, :orange, :lightgray)); [mix_with_white(colorant"limegreen", 1-α) for α ∈ firm_df.δ_pc]],
		edge_attr,
		arrow_attr,
		extend_limits = 0.1,
		edge_plottype = :beziersegments,
		kwargs...
		#axis = (; title = label(n))
	)

	nothing
end

# ╔═╡ af1b2d6d-0ac8-4bbf-b0da-adf65f27330e
function add_legend!(figpos; kwargs...)
	dict = [:lightgray => "solvent", :orange => "insolvent", :red => "bankrupt"]
	
	elements = [
		MarkerElement(marker = :circle, strokewidth=1, markersize = 20, color = c) for c in first.(dict)
	]

	Legend(figpos, elements, last.(dict); kwargs...)

	figpos
end

# ╔═╡ ba5c477d-ad8c-44f1-b815-529bd247274b
function visualize_bank_firm_network(IM, shares, bank_df, firm_df; figure = figure(320, 220), add_legend=false, kwargs...)
	fig = Figure(; figure...)
	visualize_bank_firm_network!(Axis(fig[1,1]), IM, shares, bank_df, firm_df; kwargs...)

	if add_legend
		add_legend!(fig[0,:], orientation=:horizontal, framevisible=false)
	end
	fig # |> as_svg
end

# ╔═╡ 517b4257-4097-4ab4-8ddb-34e877f6ad15
#=╠═╡
let
	(; IM, shares, bank_df, firm_df, show_firms) = RES
	
	if show_firms
		kws = (;)
	else
		kws = (; layout = componentwise_circle)
		shares = fill(0.0, (0, size(shares, 2)))
	end

	visualize_bank_firm_network(IM, shares, bank_df, firm_df; add_legend=true, kws...) |> as_svg

end
  ╠═╡ =#

# ╔═╡ 8d9c009d-7434-4d8b-884f-2725639b2748
let
	fig = Figure()
	add_legend!(fig[1,1], orientation = :horizontal, framevisible=false)
	fig
end

# ╔═╡ ed67cfcb-b83d-4b8b-b698-00e91056523a
md"""
# Visualize bank balance sheet
"""

# ╔═╡ 2a5f0419-3cf5-483a-921d-f53a49e3c096
function balance_sheet_palette()
	color_df = DataFrame(
		color = ["external", "interbank", "firm", "liquidated", "shortfall"],
		i_color = [1, 2, 4, 3, 5]
	)
	@transform!(color_df, :wong = Makie.wong_colors()[:i_color])
	color = color_df.color .=> color_df.wong
end

# ╔═╡ 0f81552e-c4cf-495e-b1c7-782fe2d7b204
#=╠═╡
md"""
* size of shock $(@bind ε_cash Slider(0.0:0.05:1.0, show_value = true, default = 0.0))
* show illiquid firm value ``A`` $(@bind show_illiquid CheckBox(default = false))
* recovery rate ``ζ`` $(@bind recovery_rate Slider(0:0.1:0.5, show_value = true, default = 0.0))
""" # |> aside
  ╠═╡ =#

# ╔═╡ 789efa32-9af5-4325-ae94-4bae9e9ecbfc
function _balance_sheet_df_((; c, x, div, ill, rec, ν_paid, y_paid, shortfall))
	DataFrame([
		(color = "external",  side = "receivable", ill=false, val=c, lab="cash 𝑐"),
		(color = "interbank", side = "receivable", ill=false, val=x, lab="IB deposit 𝑥"),
		(color = "firm",   side = "receivable", ill=false, val=div, lab="dividend δ"),
		
		(color = "firm",   side = "receivable", ill=true,  val=ill, lab="illiquid"),
		(color = "liquidated", side = "receivable", ill=false, val=rec, lab=""),
		
		(color = "external",  side = "payable",    ill=false, val=ν_paid, lab="deposits ν"),
		(color = "interbank", side = "payable",    ill=false, val=y_paid, lab="IB debt 𝑦"),
		(color = "shortfall", side = "payable",    ill=false,  val=shortfall, lab=""),
	])
end

# ╔═╡ 82ecd66b-d854-455e-b57e-42a318cd66cd
#=╠═╡
function balance_sheet_df_new(bank, firm, (; x, y); ωᵢ=1.0, ε=0.0, ℓ̄₋ᵢ=0.0)

	(; y_pc, ν_pc, ℓ, rec_field, rec_machines, div) = 
		_repay_(bank, (; x̄=x, ȳ=y), firm, (; ωᵢ, ε, ℓ̄₋ᵢ))

	y_paid = y_pc * y
	ν_paid = ν_pc * ν
	shortfall = (1-ν_pc) * ν + (1-y_pc) * y
	rec =  ωᵢ * (rec_field + rec_machines)
	ill = (1-ℓ) * ωᵢ * firm.A
	div = ωᵢ * div
	
	df = _balance_sheet_df_((; c, x, div, ill, rec, ν_paid, y_paid, shortfall))

	(; df, ℓ, shortfall)
end
  ╠═╡ =#

# ╔═╡ 05688f6e-afa8-407b-867b-884a753c4ea9
#=╠═╡
function visualize_simple_balance_sheet((; c, ν), firm, (; x, y))

	(; df, ℓ, shortfall) = balance_sheet_df_new((; c, ν), firm, (; x, y))

	if !show_illiquid || ℓ ≈ 1
		@subset!(df, :lab ≠ "illiquid")
	end
	if ℓ * firm.ζ ≈ 0
		@subset!(df, :color ≠ "liquidated")
	end
	@subset!(df, :val > 0)
	
	plt = data(df) * mapping(
		:side => sorter("receivable", "payable") => "", 
		:val,
		stack=:color, color=:color => "",
		bar_labels=:lab => verbatim
	) * visual(BarPlot, flip_labels_at = 0)
	
	draw(
		plt, 
		palettes = (; color=balance_sheet_palette()),
		legend = (; position = :top, titleposition=:left, framevisible=false),
	)
end
  ╠═╡ =#

# ╔═╡ 6e9977c1-b74e-4092-a858-1c1516bd3fe5
#=╠═╡
let
	c = max(0.7 - ε_cash, 0.0)
	ν = 1.2
	x = 1.0
	a = 0.1
	y = 1.0
	A = 1.0
	ζ = recovery_rate

	firm = (; a, A, ζ, α=0.0, ζ_α = 1.0)
	visualize_simple_balance_sheet((; c, ν), firm, (; x, y))
	
end
  ╠═╡ =#

# ╔═╡ 1a7cbefb-2e2b-4cc5-81e7-caccc970e29a
#=╠═╡
let
	c = 0.7
	ν = 1.2
	x = 1.0
	a = 0.5
	y = 1.0
	A = 1.0
	ζ = recovery_rate

	c = max(c - ε_cash, zero(c))

	firm = (; a, A, ζ, α = 0.0, ζ_α = 1.0)

	(; df, ℓ) = balance_sheet_df_new((; c, ν), firm, (; x, y))
	
	if !show_illiquid || (1-ℓ) * A ≈ 0
		@subset!(df, :lab ≠ "illiquid")
	end
	if ℓ * ζ ≈ 0
		@subset!(df, :color ≠ "liquidated")
	end
	@subset!(df, :val > 0)
	
	fig = Figure(; font="CMU", resolution=(350, 300))
	ax = Axis(fig[2,1], xticks = (1:2, ["receivables", "payables"]), limits = (nothing, nothing, 0, 3.5))

	color_df = DataFrame(
		color = ["external", "interbank", "firm", "liquidated", "shortfall"],
		i_color = [1, 2, 4, 3, 5]
	)
	
	@chain df begin
		leftjoin(_, color_df, on=:color)
		disallowmissing!
		@transform!(:i_side = :side == "receivable" ? 1 : 2)
		@transform!(:fill = (Makie.wong_colors()[:i_color], 0.5 + 0.5 * (:ill == false)))
		barplot!(ax, _.i_side, _.val, stack = _.i_color, color = _.fill, bar_labels = _.lab, flip_labels_at = 0
		)
	end

	df = @chain color_df begin
		@subset(:color ∈ ["external", "interbank", "firm"])
		@transform!(:legend = PolyElement(color = Makie.wong_colors()[:i_color]))
	end

	Legend(fig[1,1], df.legend, df.color, orientation = :horizontal)

	fig |> as_svg
end
  ╠═╡ =#

# ╔═╡ ca99d6ff-1ab1-4de1-9009-2ebdbb49140e
function balance_sheet_df((; ν, c), (; a, A, ζ), (; ℓᵢ, ωᵢ, y_pc, ν_pc, x̄, ȳ, div, rec_field, rec_machines))

	x = x̄
	y = ȳ / y_pc
	y_paid = y * y_pc
	ν_paid = ν_pc * ν
	shortfall = (1-ν_pc) * ν + (1-y_pc) * y
	rec =  ωᵢ * (rec_field + rec_machines)
	ill = (1-ℓᵢ) * ωᵢ * A
	div = ωᵢ * div
	
	df = _balance_sheet_df_((; c, x, div, ill, rec, ν_paid, y_paid, shortfall))
	
	bs_max = max(ωᵢ * (A + a) + x̄ + c, ȳ + ν) * 1.05
	(; df, bs_max)
end

# ╔═╡ 5e65f0a8-93c9-4698-8cf3-30df9992334a
#=╠═╡
function visualize_balance_sheets!(figpos, bank_df, banks, firms, shares)

	function firm(i)
		i_firm = findall(shares[:,i] .> 0) |> only
		ωᵢ = shares[i_firm,i]
		firm = firms[i_firm]
		(; ωᵢ, firm)
	end

	f_sh = firm.(1:length(banks)) |> StructArray

	bank_df = @chain bank_df begin
		rename(:ℓ => :ℓᵢ)
		@transform!(:ωᵢ = @bycol f_sh.ωᵢ)
	end
	
	bank_dfs = [
		balance_sheet_df(banks[i], f_sh.firm[i], bank_df[i,:]) for i ∈ 1:length(banks)
	] |> StructArray

	combined_bank_df = vcat(
		bank_dfs.df..., source = :bank
	)

	@transform!(combined_bank_df, :bank = "bank $(:bank)")
	
	bs_max = maximum(bank_dfs.bs_max)

	plt = @chain combined_bank_df begin
		data(_) * mapping(
			:side => sorter("receivable", "payable") => "",
			:val => "",
			color = :color => "", stack = :color,
			layout = :bank
		) * visual(BarPlot)
	end

	fg = draw!(figpos[1,1], plt, 
		axis = (limits = (nothing, nothing, 0, 1.05 * bs_max),),
		palettes=(color=balance_sheet_palette(),)
	)
	legend!(figpos[1,2], fg)

	nothing
end
  ╠═╡ =#

# ╔═╡ 7e46af06-58f4-44d3-98c0-f9b41f4d8d84
#=╠═╡
let
	(; IM, shares, firms, banks, bank_df, firm_df) = RES

	fig = Figure(resolution = (800, 300))
		
	visualize_bank_firm_network!(Axis(fig[1,1]), IM, shares, bank_df, firm_df; start = 1/n_banks)

	visualize_balance_sheets!(fig[1,2:4], bank_df, banks, firms, shares)

	fig |> as_svg
end
  ╠═╡ =#

# ╔═╡ a15f8d3c-2bb2-40d0-9439-701fbd8dd47c
md"""
# Appendix
"""

# ╔═╡ 5e26f5e4-0af9-4c2c-bba7-cfa53547e5e9
#=╠═╡
TableOfContents()
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
AlgebraOfGraphics = "cbdf2221-f076-402e-a563-3d30da359d67"
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
Chain = "8be319e6-bccf-4806-a6f7-6fae938471bc"
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
DataFrameMacros = "75880514-38bc-4a95-a458-c2aea5a3a702"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DocStringExtensions = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
Graphs = "86223c79-3864-5bf0-83f7-82e725a168b6"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
NetworksUtils = "4943429a-ba68-4c19-ade3-7332adbb3997"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Roots = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
SimpleWeightedGraphs = "47aef6b3-ad0c-573a-a1e2-d07658019622"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"

[compat]
AlgebraOfGraphics = "~0.6.12"
CairoMakie = "~0.9.3"
Chain = "~0.5.0"
Colors = "~0.12.8"
DataFrameMacros = "~0.4.0"
DataFrames = "~1.4.3"
DocStringExtensions = "~0.9.2"
Graphs = "~1.7.4"
Makie = "~0.18.3"
NetworksUtils = "~0.1.2"
PlutoUI = "~0.7.48"
Roots = "~2.0.8"
SimpleWeightedGraphs = "~1.2.1"
StructArrays = "~0.6.13"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.3"
manifest_format = "2.0"
project_hash = "d677988e84a9aa754c3f9a22c8986268d4aea6a5"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "69f7020bd72f069c219b5e8c236c1fa90d2cb409"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.2.1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "52b3b436f8f73133d7bc3a6c71ee7ed6ab2ab754"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.3"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "195c5505521008abea5aee4f96930717958eac6f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.4.0"

[[deps.AlgebraOfGraphics]]
deps = ["Colors", "Dates", "Dictionaries", "FileIO", "GLM", "GeoInterface", "GeometryBasics", "GridLayoutBase", "KernelDensity", "Loess", "Makie", "PlotUtils", "PooledArrays", "RelocatableFolders", "StatsBase", "StructArrays", "Tables"]
git-tree-sha1 = "f4d6d0f2fbc6b2c4a8eb9c4d47d14b9bf9c43d23"
uuid = "cbdf2221-f076-402e-a563-3d30da359d67"
version = "0.6.12"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e81c509d2c8e49592413bfb0bb3b08150056c79d"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Automa]]
deps = ["Printf", "ScanByte", "TranscodingStreams"]
git-tree-sha1 = "d50976f217489ce799e366d9561d56a98a30d7fe"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "0.8.2"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.CairoMakie]]
deps = ["Base64", "Cairo", "Colors", "FFTW", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "SHA", "SnoopPrecompile"]
git-tree-sha1 = "20bd6ace08bb83bf5579e8dfb0b1e23e33518b04"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.9.3"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.Chain]]
git-tree-sha1 = "8c4920235f6c561e401dfe569beb8b924adad003"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.5.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e7ff6cadf743c098e08fca25c91103ee4303c9bb"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.6"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON", "Test"]
git-tree-sha1 = "61c5334f33d91e570e1d0c3eb5465835242582c4"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "1fd869cc3875b57347f7027521f561cf46d1fcd8"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.19.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.CommonSolve]]
git-tree-sha1 = "9441451ee712d1aec22edad62db1a9af3dc8d852"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.3"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "aaabba4ce1b7f8a9b34c015053d3b1edf60fa49c"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.4.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fb21ddd70a051d882a1686a5a550990bbe371a95"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.4.1"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "e08915633fcb3ea83bf9d6126292e5bc5c739922"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.13.0"

[[deps.DataFrameMacros]]
deps = ["DataFrames", "MacroTools"]
git-tree-sha1 = "92ae44e8d08667be722ca197c97e60bcff1db968"
uuid = "75880514-38bc-4a95-a458-c2aea5a3a702"
version = "0.4.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SnoopPrecompile", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "0f44494fe4271cc966ac4fea524111bef63ba86c"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.4.3"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.Dictionaries]]
deps = ["Indexing", "Random", "Serialization"]
git-tree-sha1 = "e82c3c97b5b4ec111f3c1b55228cebc7510525a2"
uuid = "85a47980-9c8c-11e8-2b9f-f7ca1fa99fb4"
version = "0.3.25"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "7fe1eff48e18a91946ff753baf834ff4d5c03744"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.78"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "c36550cb29cbe373e95b3f40486b9a4148f89ffd"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.2"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.Extents]]
git-tree-sha1 = "5e1e4c53fa39afe63a7d356e30452249365fba99"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.1"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "90630efff0894f8142308e334473eba54c433549"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.5.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "7be5f99f7d15578798f338f5433b6c432ea8037b"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "802bfc139833d2ba893dd9e62ba1767c88d708ae"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.5"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "cabd77ab6a6fdff49bfd24af2ebe76e6e018a2b4"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.0.0"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics"]
git-tree-sha1 = "38a92e40157100e796690421e34a11c107205c86"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLM]]
deps = ["Distributions", "LinearAlgebra", "Printf", "Reexport", "SparseArrays", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "StatsModels"]
git-tree-sha1 = "884477b9886a52a84378275737e2823a5c98e349"
uuid = "38e38edf-8417-5370-95a0-9cbb8c7f171a"
version = "1.8.1"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "fb28b5dc239d0174d7297310ef7b84a11804dfab"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.0.1"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "fe9aea4ed3ec6afdfbeb5a4f39a2208909b162a6"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.5"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "fb83fbe02fe57f2c068013aa94bcdf6760d3a7a7"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+1"

[[deps.GraphMakie]]
deps = ["GeometryBasics", "Graphs", "LinearAlgebra", "Makie", "NetworkLayout", "StaticArrays"]
git-tree-sha1 = "693642d05cc34f336c3b01e7b722024ab308ac78"
uuid = "1ecd5474-83a3-4783-bb4f-06765db800d2"
version = "0.4.3"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "ba2d094a88b6b287bd25cfa86f301e7693ffae2f"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.7.4"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "678d136003ed5bceaab05cf64519e3f956ffa4ba"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.9.1"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions", "Test"]
git-tree-sha1 = "709d864e3ed6e3545230601f94e11ebc65994641"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.11"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "acf614720ef026d38400b3817614c45882d75500"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.4"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "342f789fd041a55166764c351da1710db97ce0e0"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.6"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "87f7662e03a649cffa2e05bf19c303e168732d3e"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.2+0"

[[deps.Indexing]]
git-tree-sha1 = "ce1566720fd6b19ff3411404d4b977acd4814f9f"
uuid = "313cdc1a-70c2-5d6a-ae34-0150d3930a38"
version = "1.1.1"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "842dd89a6cb75e02e85fdd75c760cdc43f5d6863"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.6"

[[deps.IntervalSets]]
deps = ["Dates", "Random", "Statistics"]
git-tree-sha1 = "3f91cd3f56ea48d4d2a75c2a65455c5fc74fa347"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.3"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "a77b273f1ddec645d1b7c4fd5fb98c8f90ad10a5"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "9816b296736292a80b9a3200eb7fbb57aaa3917a"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.5"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Loess]]
deps = ["Distances", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "46efcea75c890e5d820e670516dc156689851722"
uuid = "4345ca2d-374a-55d4-8d30-97f9976e7612"
version = "0.5.4"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "94d9c52ca447e23eac0c0f074effbcd38830deb5"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.18"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "2ce8695e1e699b68702c03402672a69f54b8aca9"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.2.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Makie]]
deps = ["Animations", "Base64", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "Distributions", "DocStringExtensions", "FFMPEG", "FileIO", "FixedPointNumbers", "Formatting", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageIO", "InteractiveUtils", "IntervalSets", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MakieCore", "Markdown", "Match", "MathTeXEngine", "MiniQhull", "Observables", "OffsetArrays", "Packing", "PlotUtils", "PolygonOps", "Printf", "Random", "RelocatableFolders", "Serialization", "Showoff", "SignedDistanceFields", "SnoopPrecompile", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun"]
git-tree-sha1 = "d3b9553c2f5e0ca588e4395a9508cef024bd9e8a"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.18.3"

[[deps.MakieCore]]
deps = ["Observables"]
git-tree-sha1 = "c1885d865632e7f37e5a1489a164f44c54fb80c9"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.5.2"

[[deps.MappedArrays]]
git-tree-sha1 = "e8b359ef06ec72e8c030463fe02efe5527ee5142"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.1"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.Match]]
git-tree-sha1 = "1d9bc5c1a6e7ee24effb93f175c9342f9154d97f"
uuid = "7eb4fadd-790c-5f42-8a69-bfa0b872bfbf"
version = "1.2.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "Test", "UnicodeFun"]
git-tree-sha1 = "f04120d9adf4f49be242db0b905bea0be32198d1"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.5.4"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.MiniQhull]]
deps = ["QhullMiniWrapper_jll"]
git-tree-sha1 = "9dc837d180ee49eeb7c8b77bb1c860452634b0d1"
uuid = "978d7f02-9e05-4691-894f-ae31a51d76ca"
version = "0.4.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "b34e3bc3ca7c94914418637cb10cc4d1d80d877d"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.3"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore"]
git-tree-sha1 = "18efc06f6ec36a8b801b23f076e3c6ac7c3bf153"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.0.2"

[[deps.NetworkLayout]]
deps = ["GeometryBasics", "LinearAlgebra", "Random", "Requires", "SparseArrays"]
git-tree-sha1 = "cac8fc7ba64b699c678094fa630f49b80618f625"
uuid = "46757867-2c16-5918-afeb-47bfcb05e46a"
version = "0.4.4"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NetworksUtils]]
deps = ["GraphMakie", "Graphs", "InteractiveUtils", "Makie", "Markdown", "NetworkLayout", "SimpleWeightedGraphs", "Statistics"]
git-tree-sha1 = "ffdf14f92843b4a125a7351494f88e3f504a7f1f"
uuid = "4943429a-ba68-4c19-ade3-7332adbb3997"
version = "0.1.2"

[[deps.Observables]]
git-tree-sha1 = "6862738f9796b3edc1c09d0890afce4eca9e7e93"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.4"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "f71d8950b724e9ff6110fc948dff5a329f901d64"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.8"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "923319661e9a22712f24596ce81c54fc0366f304"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.1+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6e9dba33f9f2c44e08a020b0caf6903be540004"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.19+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.40.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "cf494dca75a69712a72b80bc48f59dcf3dea63ec"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.16"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "f809158b27eba0c18c269cf2a2be6ed751d3e81d"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.3.17"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "1155f6f937fa2b94104162f01fa400e192e4272f"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.4.2"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "03a7a85b76381a3d04c7a1656039197e70eda03d"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.11"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "84a314e3926ba9ec66ac097e3635e270986b0f10"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.50.9+0"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "b64719e8b4504983c7fca6cc9db3ebc8acc2a4d6"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.1"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f6cf8e7944e50901594838951729a1861e668cb8"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.2"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "SnoopPrecompile", "Statistics"]
git-tree-sha1 = "21303256d239f6b484977314674aef4bb1fe4420"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.1"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "efc140104e6d0ae3e7e30d56c98c4a927154d684"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.48"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "LaTeXStrings", "Markdown", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "d8ed354439950b34ab04ff8f3dfd49e11bc6c94b"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "d7a7aef8f8f2d537104f170139553b14dfe39fe9"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.2"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.QhullMiniWrapper_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Qhull_jll"]
git-tree-sha1 = "607cf73c03f8a9f83b36db0b86a3a9c14179621f"
uuid = "460c41e3-6112-5d7f-b78c-b6823adb3f2d"
version = "1.0.0+1"

[[deps.Qhull_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "695c3049ad94fa38b7f1e8243cdcee27ecad0867"
uuid = "784f63db-0788-585a-bace-daefebcd302b"
version = "8.0.1000+0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "97aa253e65b784fd13e83774cadc95b38011d734"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.6.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.Roots]]
deps = ["ChainRulesCore", "CommonSolve", "Printf", "Setfield"]
git-tree-sha1 = "a3db467ce768343235032a1ca0830fc64158dadf"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.0.8"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
git-tree-sha1 = "bc12e315740f3a36a6db85fa2c0212a848bd239e"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.4.2"

[[deps.ScanByte]]
deps = ["Libdl", "SIMD"]
git-tree-sha1 = "2436b15f376005e8790e318329560dcc67188e84"
uuid = "7b38b023-a4d7-4c5e-8d43-3f3097f304eb"
version = "0.3.3"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.ShiftedArrays]]
git-tree-sha1 = "503688b59397b3307443af35cd953a13e8005c16"
uuid = "1277b4bf-5013-50f5-be3d-901d8477a67a"
version = "2.0.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleWeightedGraphs]]
deps = ["Graphs", "LinearAlgebra", "Markdown", "SparseArrays", "Test"]
git-tree-sha1 = "a6f404cc44d3d3b28c793ec0eb59af709d827e4e"
uuid = "47aef6b3-ad0c-573a-a1e2-d07658019622"
version = "1.2.1"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "8fb59825be681d451c246a795117f317ecbcaa28"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.2"

[[deps.SnoopPrecompile]]
git-tree-sha1 = "f604441450a3c0569830946e5b33b78c928e1a85"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "4e051b85454b4e4f66e6a6b7bdc452ad9da3dcf6"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.10"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "5783b877201a82fc0014cbf381e7e6eb130473a4"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.0.1"

[[deps.StatsModels]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Printf", "REPL", "ShiftedArrays", "SparseArrays", "StatsBase", "StatsFuns", "Tables"]
git-tree-sha1 = "a5e15f27abd2692ccb61a99e0854dfb7d48017db"
uuid = "3eaba693-59b7-5ba5-a881-562e759f1c8d"
version = "0.6.33"

[[deps.StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArraysCore", "Tables"]
git-tree-sha1 = "13237798b407150a6d2e2bce5d793d7d9576e99e"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.13"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "c79322d36826aa2f4fd8ecfa96ddb47b174ac78d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "f8cd5b95aae14d3d88da725414bdde342457366f"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.6.2"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "8a75929dcd3c38611db2f8d08546decb514fcadf"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.9"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.URIs]]
git-tree-sha1 = "e59ecc5a41b000fa94423a578d29290c7266fc10"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "58443b63fb7e465a8a7210828c91c08b92132dff"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.14+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "libpng_jll"]
git-tree-sha1 = "d4f63314c8aa1e48cd22aa0c17ed76cd1ae48c3c"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.3+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"
"""

# ╔═╡ Cell order:
# ╟─ec9e91ef-fff8-4d7f-8df3-86df19ced7ed
# ╟─213aa0db-1aac-4471-b7eb-42f14d23a5c2
# ╠═7d711677-35dd-46d2-ba87-f2994d483dc3
# ╠═b90bc5c2-5d7c-4e0c-bf71-08b1af48e754
# ╠═b6238421-76d2-4da0-a1e0-d0e94bed1b0d
# ╠═2e88a6a7-79a5-4ef2-b3f4-ac916e09edf6
# ╟─8eb2aa7c-2dc1-11ed-24ab-e9cce1f5013d
# ╟─3c867227-aa93-4ade-9434-a1a61904656a
# ╟─8fd6c40f-17b1-461b-82ec-100c6c987701
# ╠═eae2ef96-deeb-4e17-85fe-00ef53d13a46
# ╠═f5f971e9-bb3c-44f2-b269-bdaf17bff943
# ╟─3d741a24-60c5-4462-b048-796efe3ef0e4
# ╟─24fbe0a5-f5dd-49b3-a667-963b25168191
# ╟─bfed935f-0b08-491a-a29a-c118c35134d6
# ╟─8709bbc9-ee03-42ea-ac7c-e2b1c123eaed
# ╟─57247692-b967-4880-9574-af5a45fffc0b
# ╠═75410a0d-a352-435a-88d8-7b520bf2005f
# ╟─5e0ce818-4ac9-4efb-ad24-7223cac5cb12
# ╠═952825f3-7f66-4435-a850-cbe160106f6f
# ╠═58481389-ded5-4b52-9760-ece0d381451e
# ╠═1cad0937-5ecd-4b0a-871c-de8c2f75d7b7
# ╠═c59a6cc9-f9bd-4831-9bd8-7bc331d5790f
# ╠═2da46c3d-dd64-424f-8fe9-b6ba1156f423
# ╠═62a50fdf-4c81-4adf-9366-40362cbefb83
# ╠═63ed5af7-2c4d-4184-8559-ec6977f04709
# ╟─f49f7e4d-20bc-4dd6-9586-207dbfd19cef
# ╠═6b31b940-0a2e-4bbf-9534-bdfa5a02c3c4
# ╠═d5ac4a10-48f0-4e59-a186-eff9dd059717
# ╠═162e1da7-0d14-464f-8f64-66804009bc5d
# ╠═c46f9af5-6bd5-49fe-8086-44b238d5a8a1
# ╟─afa36073-424e-4ff8-b584-554ae66cee06
# ╠═0a2dfc44-67f8-4e9f-991f-906509116d66
# ╟─51f5c242-54af-4511-9838-58326756197d
# ╠═1d00ca98-985a-489c-9ef4-d74f558d9b71
# ╠═517b4257-4097-4ab4-8ddb-34e877f6ad15
# ╠═eabcafbd-5c35-4ee9-a50c-d129f8d3a34a
# ╠═d5dc8002-0d8e-4681-a6a5-4d5fc49f9288
# ╠═fc4fa10d-2b39-4f10-bac8-7a1e5467669d
# ╠═fcb09a02-70ca-4da0-bab9-6f2d9da6a195
# ╠═3d037a99-7978-4fce-b603-535a90362f40
# ╠═d816a40a-67c1-4a1a-ae0c-d26502ba3918
# ╠═72de5ce2-6bab-4f9a-b45a-4559fa3affea
# ╟─cbcdef2e-d854-47c9-bfc4-595c90ea50fb
# ╠═ba5c477d-ad8c-44f1-b815-529bd247274b
# ╠═92c4efe9-6ba1-4895-a3c8-05ec312a3820
# ╠═8028be99-7559-4974-91a6-36ccab2d4a7f
# ╠═af1b2d6d-0ac8-4bbf-b0da-adf65f27330e
# ╠═8d9c009d-7434-4d8b-884f-2725639b2748
# ╟─ed67cfcb-b83d-4b8b-b698-00e91056523a
# ╠═d004095e-2910-4a7c-8ee4-b5d13c34419e
# ╠═6e9977c1-b74e-4092-a858-1c1516bd3fe5
# ╠═05688f6e-afa8-407b-867b-884a753c4ea9
# ╠═2a5f0419-3cf5-483a-921d-f53a49e3c096
# ╟─0f81552e-c4cf-495e-b1c7-782fe2d7b204
# ╠═1a7cbefb-2e2b-4cc5-81e7-caccc970e29a
# ╠═789efa32-9af5-4325-ae94-4bae9e9ecbfc
# ╠═82ecd66b-d854-455e-b57e-42a318cd66cd
# ╠═ca99d6ff-1ab1-4de1-9009-2ebdbb49140e
# ╠═7e12dd72-2d66-4d94-8af6-084e80a63885
# ╠═fde5439e-d1ed-43e7-be89-d8f41eff3817
# ╠═7e46af06-58f4-44d3-98c0-f9b41f4d8d84
# ╠═5e65f0a8-93c9-4698-8cf3-30df9992334a
# ╠═1b470abb-bf1e-4ec6-a5c8-5f05e86c87dc
# ╟─a15f8d3c-2bb2-40d0-9439-701fbd8dd47c
# ╠═57945142-4153-4c89-bef6-bb8bb0b6a7d8
# ╠═5e26f5e4-0af9-4c2c-bba7-cfa53547e5e9
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
