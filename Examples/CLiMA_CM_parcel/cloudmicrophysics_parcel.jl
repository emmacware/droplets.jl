import OrdinaryDiffEq as ODE
import CairoMakie as MK
import Thermodynamics as TD
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import ClimaParams as CP
using Droplets
using ComponentArrays

include(joinpath(pkgdir(CM), "parcel", "Parcel.jl"))
include("helperfunctions.jl")

FT = Float32
Ns::Int = 2^14

# Get free parameters
tps = TD.Parameters.ThermodynamicsParameters(FT)
wps = CMP.WaterProperties(FT)
# Constants
ρₗ = wps.ρw
ρᵢ = wps.ρi
R_v = TD.Parameters.R_v(tps)
R_d = TD.Parameters.R_d(tps)

# Common Initial conditions
Nₐ = FT(0)
Nₗ = FT(200 * 1e6)
Nᵢ = FT(0)
r₀ = FT(8e-6)
p₀ = FT(800 * 1e2)
T₀ = FT(273.15 + 7.0)
ln_INPC = FT(0)
e = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
Sₗ = FT(1)
md_v = (p₀ - e) / R_d / T₀
mv_v = e / R_v / T₀

# Simulation parameters passed into ODE solver
w = FT(10)                                 # updraft speed
const_dt = FT(0.5)                         # model timestep
t_max = FT(20)

run_parameters = [
    ["Monodisperse", "Condensation"],
    ["Gamma", "Condensation"],
    ["Monodisperse", "Superdroplets"],
    ["Gamma", "Superdroplets"]
]

# Data from Rogers(1975) Figure 1
# https://www.tandfonline.com/doi/abs/10.1080/00046973.1975.9648397
#! format: off
Rogers_time_supersat = [0.0645, 0.511, 0.883, 1.4, 2.07, 2.72, 3.24, 3.89, 4.53, 5.87, 7.16, 9.79, 16.0, 19.8]
Rogers_supersat = [0.0268, 0.255, 0.393, 0.546, 0.707, 0.805, 0.863, 0.905, 0.938, 0.971, 0.978, 0.963, 0.910, 0.885]
Rogers_time_radius = [0.561, 2, 3.99, 10.7, 14.9, 19.9]
Rogers_radius = [8.0, 8.08, 8.26, 8.91, 9.26, 9.68]

# Setup the plots
fig = MK.Figure(size = (800, 400))
ax1 = MK.Axis(fig[2, 1], xlabel = "Time [s]", ylabel = "Supersaturation [%]")
ax2 = MK.Axis(fig[2, 2], xlabel = "Time [s]", ylabel = "radius [μm]")
ax3 = MK.Axis(fig[2, 3], xlabel = "Time [s]", ylabel = "q_liq [g/kg]")

#! format: on

for runs in run_parameters

    local pp = parcel_params{FT}(
        liq_size_distribution = runs[1],
        condensation_growth = runs[2],
        const_dt = const_dt,
        w = w,
    )
    # solve ODE
    coagsettings = coag_settings{FT}(Ns=Ns,Δt=1,ΔV = 1,n0=Nₗ,R0=r₀)
    if runs[1] == "Monodisperse"
        drops = init_monodisperse(coagsettings)
    else
        drops = init_ξ_const(coagsettings)
    end

    if runs[2] == "Superdroplets"
        ml_v = sum(drops.X.*drops.ξ.*ρₗ)
    else
        ml_v = Nₗ * 4 / 3 * FT(π) * ρₗ * r₀^3
    end
    
    qᵥ = mv_v / (md_v + mv_v + ml_v)
    qₗ = ml_v / (md_v + mv_v + ml_v)
    qᵢ = FT(0)
    IC = [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, ln_INPC]
    d_Y = ComponentVector{FT}(IC=IC,R =drops.R,ξ =drops.ξ, X =drops.X)

    sol = run_parcel_sd(d_Y, FT(0), t_max, pp)

    DSD = runs[1]
    if runs[2] == "Condensation"
        model = "CLiMA"
    elseif runs[2] == "Superdroplets"
        model = "Droplets"
        if runs[1] == "Gamma"
            DSD = "Exponential"
        end
    end

    # Plot results
    MK.lines!(ax1, sol.t, (sol[1, :] .- 1) * 100.0, label = model * " " * DSD *"Dist")
    MK.lines!(ax3, sol.t, sol[5, :] * 1e3)

    sol_Nₗ = sol[8, :]
    sol_Nᵢ = sol[9, :]
    # Compute the current air density
    sol_T = sol[3, :]
    sol_p = sol[2, :]
    sol_qᵥ = sol[4, :]
    sol_qₗ = sol[5, :]
    sol_qᵢ = sol[6, :]
    local q = TD.PhasePartition.(sol_qᵥ + sol_qₗ + sol_qᵢ, sol_qₗ, sol_qᵢ)
    local ts = TD.PhaseNonEquil_pTq.(tps, sol_p, sol_T, q)
    local ρₐ = TD.air_density.(tps, ts)
    # Compute the mean particle size based on the distribution
    distr = sol.prob.p.liq_distr
    moms = distribution_moments.(distr, sol_qₗ, sol_Nₗ, ρₗ, ρₐ)
    local r = similar(sol_T)

    if pp.condensation_growth == "Superdroplets"
        for (i,output) in enumerate(sol.u)
            r[i] = (3/(4pi)*sum(output.X)/Ns)^(1/3)
        end
    else 
        for it in range(1, length(sol_T))
            r[it] = moms[it].r
        end
    end
    MK.lines!(ax2, sol.t, r * 1e6)
end

MK.lines!(ax1, Rogers_time_supersat, Rogers_supersat, label = "Rogers_1975", color = :black)
MK.lines!(ax2, Rogers_time_radius, Rogers_radius, color = :black)

legend = MK.axislegend(
    ax1,
    framevisible = false,
    labelsize = 12,
    orientation = :horizontal,
)
fig[1,:] = legend

MK.display(fig)
MK.save("liquid_only_parcel.svg", fig)
nothing
