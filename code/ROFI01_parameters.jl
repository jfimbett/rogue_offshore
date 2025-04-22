# ********************************************************************************************
# * Initial Parameters of the Model
# * Last updated: 21‑Apr‑2025
# ********************************************************************************************

using Parameters

@with_kw struct ModelParameters
    # ── Technology & shocks ────────────────────────────────────────────────────────────
    α::Float64   = 0.60   # Capital share (Cobb‑Douglas) – mid‑point of 0.55‑0.65 range
    ρ::Float64   = 0.90   # AR(1) persistence of the observable productivity shock z
    σ::Float64   = 0.12   # Std. dev. of ε (so quarterly σ≈0.06 in logs)
    η::Float64   = 0.20   # Magnitude of private shock (≈ 20 % of avg. TFPR)
    κ::Float64   = 0.50   # Prob{η = +η}  (symmetric for baseline)
    δ::Float64   = 0.10   # Depreciation rate of physical capital
    θ::Float64   = 2.00   # Quadratic investment‐adjustment cost coefficient

    # ── Financing block ────────────────────────────────────────────────────────────────
    r::Float64   = 0.04   # Risk‑free rate (annual, real terms)
    ϕ::Float64   = 0.05   # External‑finance premium (5 ¢ per €1 raised)

    # ── Compensation & governance ──────────────────────────────────────────────────────
    ν::Float64   = 0.02   # Bonus share of accounting earnings
    β::Float64   = 0.01   # Manager’s equity stake
    λ_x::Float64 = 0.25   # Convex cost coeff. for (possibly illegal) transfer‑pricing
    λ_s::Float64 = 0.10   # Convex cost coeff. for outright expropriation
    l::Float64   = 2.00   # Litigation/penalty multiple if expropriation is exposed

    # ── Taxation & repatriation ───────────────────────────────────────────────────────
    τ_c::Float64 = 0.25   # Statutory corporate income‑tax rate in the home country
    τ_x::Float64 = 0.10   # Effective tax on *voluntarily* repatriated offshore profits
    τ_d::Float64 = 0.40   # Punitive tax rate if the data‑leak forces disclosure
    ω::Float64   = 0.30   # Arm’s‑length ceiling: x_t / k_t ≤ ω (transfer‑pricing cap)

    γ_x::Float64 = 0.25   # Fraction of “legal” offshore balance repatriated each period
    γ_s::Float64 = 0.05   # Leakage rate from hidden expropriation account
    ξ::Float64   = 0.70   # Share of secretly repatriated funds that actually reaches SHs

    # ── Information risk ──────────────────────────────────────────────────────────────
    p::Float64   = 0.02   # Per‑period prob. that *any* hidden subsidiary leaks (2 %)
end
