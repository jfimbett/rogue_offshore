#********************************************************************************************
#* This file belief_formation.jl illustrates how agents in the model form beliefs.
#*
#* 21/04/2025
#********************************************************************************************

# =============================================================================
#  BELIEF‑UPDATING ILLUSTRATION FOR THE SECRET ACCOUNT  (21 Apr 2025)
# =============================================================================

using Parameters, Distributions, Plots, LaTeXStrings
include("ROFI01_parameters.jl")        # <- your parameter file
mp = ModelParameters()                 # instantiate

#───────────────────────── 1.  Prior g₀(B) ───────────────────────────────────
grid = range(0.0, 4.0; length = 800)
Δb   = grid[2] - grid[1]
g0   = pdf.(LogNormal(log(1.5), 0.5), grid)

left_mask  = grid .<= 1.0              # 𝔅(+η) = [0,1]
right_mask = .!left_mask

#───────────────────────── 2.  I‑projection posterior g₁ ─────────────────────
κ = mp.κ
mass_L = sum(g0[left_mask]) * Δb
mass_R = sum(g0[right_mask]) * Δb

g1 = zeros(length(grid))
g1[left_mask]  .= g0[left_mask]  * κ       / mass_L
g1[right_mask] .= g0[right_mask] * (1-κ)   / mass_R

#───────────────────────── 3.  Convolution to get g₂ ─────────────────────────
γs    = mp.γ_s                 # payout from the secret account, assume z = 2 times eta
s_pos = 2*mp.η + mp.η
s_neg = 2*mp.η - mp.η                

function shift_scale(pdf_vec, grid, Δb, shift, scale)
    N = length(grid); g = zeros(N)
    for (i,b) in enumerate(grid)
        x = (b - shift)/scale
        j = Int(floor(x/Δb)) + 1          # nearest index
        1 ≤ j ≤ N && (g[i] = pdf_vec[j] / abs(scale))
    end
    g
end

g2 = κ      .* shift_scale(g1, grid, Δb, s_pos, 1-γs) .+
     (1 - κ) .* shift_scale(g1, grid, Δb, s_neg, 1-γs)
g2 ./= sum(g2) * Δb            # renormalise

#───────────────────────── 4.  Three‑panel plot ──────────────────────────────
default(; size = (1000, 300), linewidth = 2)

p1 = plot(grid, g0, label = L"g_{0}", color = :blue,
          title = L"\textbf{Prior}")
vline!(p1, [1.0], label = "Boundary", color = :gray, linestyle = :dot)
plot!(p1, xlabel = L"B_t^{s}", ylabel = L"\text{density}")

p2 = plot(grid, g0, label = L"g_0", color = :blue, alpha = 0.25)
plot!(p2, grid, g1, label = L"g_{1}", color = :red)
vline!(p2, [1.0], label = "", color = :gray, linestyle = :dot)
annotate!(p2, 0.35, maximum(g1)*0.1, text(L"\kappa", 10))
annotate!(p2, 2.5,  maximum(g1)*0.1, text(L"1-\kappa", 10))
plot!(p2, title = L"\textbf{I-projection}",
          xlabel = L"B_t^{s}")

p3 = plot(grid, g2, label = L"g_{2}", color = :black,
          linestyle = :dash,
          title = L"\textbf{Posterior}")
plot!(p3, xlabel = L"B_{t+1}^{s}")

plot(p1, p2, p3; layout = (1,3))