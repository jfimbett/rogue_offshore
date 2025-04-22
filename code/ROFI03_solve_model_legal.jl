#********************************************************************************************
#* This file solve_model_legal.jl solves the model in the paper when tax avoidance is legal.
#*
#* 21/04/2025
#********************************************************************************************

using Parameters
using Distributions
using Plots
using LaTeXStrings
using Interpolations
include("ROFI01_parameters.jl")       
    

"""
    tauchen(nz, ρ, σ; m = 3)

Return (grid, P) where `grid` has `nz` points covering ±m·σ_z
around the mean, and `P` is the nz×nz transition matrix.
"""
function tauchen(nz, ρ, σ; m = 3)
    σ_z = σ / sqrt(1 - ρ^2)                 # unconditional s.d.
    z_max =  m * σ_z
    z_min = -m * σ_z
    grid  = LinRange(z_min, z_max, nz)
    step  = grid[2] - grid[1]

    P = zeros(nz, nz)
    for i in 1:nz
        for j in 1:nz
            if j == 1
                P[i,j] = cdf(Normal(),
                              (grid[1] - ρ*grid[i] + step/2)/σ)
            elseif j == nz
                P[i,j] = 1 -
                         cdf(Normal(),
                             (grid[end] - ρ*grid[i] - step/2)/σ)
            else
                z_low  = (grid[j]   - ρ*grid[i] - step/2)/σ
                z_high = (grid[j]   - ρ*grid[i] + step/2)/σ
                P[i,j] = cdf(Normal(), z_high) - cdf(Normal(), z_low)
            end
        end
    end
    return collect(grid), P
end

function steady_state(mp::ModelParameters)
    # Unpack parameters
    @unpack α, δ, τ_c, θ, r = mp

    # Steady state equations
    k_ss   = (((1-τ_c)*α)/(r+δ+0.5*θ*δ^2))^(1/(1-α)) 

    return k_ss
end

function d_1(k, kp, c, cp, x, s, Bx, z, η_t, mp::ModelParameters)
    # Unpack all parameters
    @unpack α, ρ, σ, η, κ, δ, θ, r, ϕ, c̄, ν, β, λ_x, λ_s, l, ν, τ_c, τ_x, τ_d, ω, γ_x, γ_s, ξ, p = mp

    i = kp-(1-δ)*k  # Investment
    Ψ = 0.5*θ*k*(i/k)^2

    # Compute dividends
    d1 = (1-τ_c)*(1-ν)*((z+η_t)*k^α-x-s)+τ_c*δ*k - i - Ψ - cp -c*(1+r*(1-τ_c)) + (1-τ_x)*γ_x*Bx 

    # if dividends are negative they must be financed with equity at a cost of ϕ
    if d1 < 0
        d1 = d1 - ϕ*abs(d1)
    end

    return d1
end

function u_1(k, kp, c, cp, x, s, Bx, Bs, z, η_t, mp::ModelParameters)
    @unpack α, ρ, σ, η, κ, δ, θ, r, ϕ, c̄, ν, β, λ_x, λ_s, l, ν, τ_c, τ_x, τ_d, ω, γ_x, γ_s, ξ, p = mp

    u_1 = ν*((z+η_t)*k^α-x-s) + ξ*γ_s*Bs - 0.5*λ_s*k*(s/k)^2 + β*d_1(k, kp, c, cp, x, s, Bx, z, η_t, mp) 

    return u_1
end

@with_kw struct GridParameters
    #-------- State Variables -----------
    nk::Int = 10   # Number of grid points for capital
    nc::Int = 10   # Number of grid points for cash
    nz::Int = 5    # Number of grid points for z
    nη::Int = 2    # Number of grid points for η
    nbx::Int = 5    # Number of grid points for Bx
    nbs::Int = 5    # Number of grid points for Bs

    #-------- Control Variables -----------
    nkp::Int = 10  # Number of grid points for kp
    ncp::Int = 10  # Number of grid points for cp
    nxs::Int = 10  # Number of grid points for x
    nss::Int = 10  # Number of grid points for s

    #-------- Other Parameters -----------
    lb_k::Float64 = 0.5  # Factor lower bound for k (lb_k*k_ss)
    ub_k::Float64 = 2.0  # Factor upper bound for k (ub_k*k_ss)
    c_bar::Float64 = 5.0  # Upper bound for cash 
end

@with_kw struct Grids 

    #-------- State Variables -----------
    k_grid::Vector{Float64} = []  # Capital grid
    c_grid::Vector{Float64} = []  # Cash grid
    z_grid::Vector{Float64} = []  # z grid
    η_grid::Vector{Float64} = []  # η grid
    Bx_grid::Vector{Float64} = [] # Bx grid
    Bs_grid::Vector{Float64} = [] # Bs grid

    #-------- Control Variables -----------
    kp_grid::Vector{Float64} = [] # kp grid
    cp_grid::Vector{Float64} = [] # cp grid
    x_grid::Vector{Float64} = []  # x grid
    s_grid::Vector{Float64} = []  # s grid

    #-------- Probability Matrices -----------
    Q_z::Matrix{Float64} = []  # Transition matrix for z
    Q_η::Matrix{Float64} = []  # Transition matrix for η

end




function create_grids(mp::ModelParameters, gp::GridParameters)
    @unpack nk, nc, nz, nη, nbx, nbs, nkp, ncp, nxs, nss, lb_k, ub_k, c_bar = gp

    #------------------------State Variables------------------------

    # Capital grid
    k_ss = steady_state(mp)
    k_min = 0.5*k_ss
    k_max = 2.0*k_ss
    k_grid = LinRange(k_min, k_max, nk)

    # cash grid
    c_min = 0.0
    c_max = c_bar
    c_grid = LinRange(c_min, c_max, nc)

    # z grid (use tauchen method)
    @unpack ρ, σ = mp
    z_grid, Q_z = tauchen(nz, ρ, σ; m = 3)
    
    # η grid, easy, η_t = {+η, -η} with prob. κ and 1-κ
    @unpack η, κ = mp
    η_grid = [-η, η]
    Q_η = [κ 1.0 - κ; 1.0 - κ κ]

    # Bx grid and Bs use same grid as cash
    Bx_grid = LinRange(c_min, c_max, nbx)
    Bs_grid = LinRange(c_min, c_max, nbs)

    #------------------------Control Variables------------------------
    kp_grid = LinRange(k_min, k_max, nkp)
    cp_grid = LinRange(c_min, c_max, ncp)
    x_grid  = LinRange(c_min, c_max, nxs)
    s_grid  = LinRange(c_min, c_max, nss)

    return Grids(k_grid=k_grid, c_grid=c_grid, z_grid=z_grid, η_grid=η_grid, 
                 Bx_grid=Bx_grid, Bs_grid=Bs_grid, kp_grid=kp_grid, 
                 cp_grid=cp_grid, x_grid=x_grid, s_grid=s_grid, 
                 Q_z=Q_z, Q_η=Q_η)

end


# solve the model

function bellman_step!(U, Kp, Cp, X, S, mp::ModelParameters, gp::GridParameters, grids::Grids)

    @unpack k_grid, c_grid, z_grid, η_grid, Bx_grid, Bs_grid, kp_grid, cp_grid, x_grid, s_grid = grids

    # interpolate U
    Uf = LinearInterpolation((k_grid, c_grid, z_grid, η_grid, Bx_grid, Bs_grid), U, extrapolation_bc=Flat()) 

    # nested loops plus enumeration of the state variables
    # loop only on state variables 
    for (i_k, k) in enumerate(k_grid)
        for (i_c, c) in enumerate(c_grid)
            for (i_z, z) in enumerate(z_grid)
                for (i_η, η_t) in enumerate(η_grid)
                    for (i_Bx, Bx) in enumerate(Bx_grid)
                        for (i_Bs, Bs) in enumerate(Bs_grid)
                           
                            # In vectorize form, we need to compute the value function for all control variables to choose the one that maximizes
                            cont = TU(k, c, x, η_t, Bx, Bs, kp_grid, cp_grid, x_grid, s_grid, Uf, mp, grids)
                            
                            # max and argmax
                            
                            max_index = argmax(cont)
                            U[i_k,i_c,i_z,i_η,i_Bx,i_Bs,i_Bx,i_Bs] = cont[max_index]
                            
                            Kp[i_k,i_c,i_z,i_η,i_Bx,i_Bs,i_Bx,i_Bs] = kp_grid[max_index]
                            Cp[i_k,i_c,i_z,i_η,i_Bx,i_Bs,i_Bx,i_Bs] = cp_grid[max_index]
                            X[i_k,i_c,i_z,i_η,i_Bx,i_Bs,i_Bx,i_Bs]  = x_grid[max_index]
                            S[i_k,i_c,i_z,i_η,i_Bx,i_Bs,i_Bx,i_Bs]  = s_grid[max_index]

                        end
                    end
                end
            end
        end
    end

end


    
# Bellman operator
function TU(k, c, x, η_t, Bx, Bs, kp_grid, cp_grid, x_grid, s_grid, Uf, mp::ModelParameters, gr::Grids)
    # kp_grid, cp_grid, x_grid, s_grid are 1‑D vectors of any numeric type
    nk, nc, nx, ns = length.( (kp_grid, cp_grid, x_grid, s_grid) )
    N = nk * nc * nx * ns                     # total number of combinations

    # ── 1. Create an iterator over the Cartesian product ─────────────────────────
    prod_iter = Iterators.product(kp_grid, cp_grid, x_grid, s_grid)

    # ── 2. Collect it into a matrix  (rows = combos, cols = variables) ───────────
    #       Step a: materialise the tuples as one long vector
    flat = vcat(collect(prod_iter)...)        # 4N‑element vector
    #       Step b: reshape and transpose to get (N × 4)
    combo_mat = reshape(flat, 4, N)'          # N rows, 4 columns

    # for every single row, apply the function u_1
    # and return a vector of the same size as the number of combinations
    u_1_vec = zeros(N)                       # preallocate
    for i in 1:N
        kp, cp, x, s = combo_mat[i, :]
        u_1_vec[i] = u_1(k, kp, c, cp, x, s, Bx, Bs, z, η_t, mp)
    end

    # now same for the continuation value, Uf is a function of the state variables that has been interpolated

    ev = zeros(N)  # preallocate
    @unpack Q_z, Q_η = gr
    for i in 1:N
        kp, cp, x, s = combo_mat[i, :]
        # expected value based on the transition matrix
        ev[i] = 0.0
        for (i_zp, zp) in enumerate(z_grid)
            for (i_ηp, ηp) in enumerate(η_grid)
                Bxp = (1-γ_x)*Bx+x
                Bsp = (1-γ_s)*Bs+s
                ev[i] += (1-p)*Q_z[i_zp, i_z] * Q_η[i_ηp, i_η] * Uf(kp, cp, zp, ηp, Bxp, Bsp) - p*Bsp
            end
        end
    end

    return u_1_vec + (1/(1+r))*ev

end

function solve_model_legal(mp::ModelParameters, gp::GridParameters)
    # Unpack parameters
    @unpack nk, nc, nz, nη, nbx, nbs, nkp, ncp, nxs, nss = gp
    # Unpack model parameters
    @unpack α, ρ, σ, η, κ, δ, θ, r, ϕ = mp

    # display the number of total states and controls
    println("Number of states: ", nk*nc*nz*nη*nbx*nbs)
    println("Number of controls: ", nkp*ncp*nxs*nss)

    # Initialize value function and policy functions, we call it U since it is different than the value of the firm
    U = zeros(nk,nc,nz,nη,nbs,nxs,nkp,ncp,nxs, nss);
    # similar for the four policy functions (use similar matrices)
    Kp = similar(U);
    Cp = similar(U);
    X  = similar(U);
    S  = similar(U);

    @time bellman_step!(U, Kp, Cp, X, S, mp::ModelParameters, gp::GridParameters, grids::Grids);
end
    
mp = ModelParameters()
gp = GridParameters()
grids = create_grids(mp, gp)

