include(joinpath("/Users/emmaware/.julia/dev/Droplets/src/SDfunc/constants.jl"))
include(joinpath("/Users/emmaware/.julia/dev/Droplets/src/SDfunc/coalescence.jl"))
include(joinpath("/Users/emmaware/.julia/dev/Droplets/src/SDfunc/condensation.jl"))
include(joinpath("/Users/emmaware/.julia/dev/Droplets/src/SDfunc/setup.jl"))
include(joinpath("/Users/emmaware/.julia/dev/Droplets/src/SDfunc/updateposition.jl"))
include(joinpath("/Users/emmaware/.julia/dev/Droplets/src/SDfunc/density.jl"))
using RecursiveArrayTools
import OrdinaryDiffEq as ODE
import CairoMakie as MK
import Thermodynamics as TD
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import ClimaParams as CP

include(joinpath(pkgdir(CM), "parcel", "Parcel.jl"))

FT = Float64

# Get free parameters
tps = TD.Parameters.ThermodynamicsParameters(FT)
wps = CMP.WaterProperties(FT)
# Constants
ρₗ = wps.ρw
ρl = wps.ρw
ρᵢ = wps.ρi
R_v = TD.Parameters.R_v(tps)
R_d = TD.Parameters.R_d(tps)

# Initial conditions
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
ml_v = Nₗ * 4 / 3 * FT(π) * ρₗ * r₀^3
qᵥ = mv_v / (md_v + mv_v + ml_v)
qₗ = ml_v / (md_v + mv_v + ml_v)
qᵢ = FT(0)
IC = [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, ln_INPC]

V = FT(1)
ξstart, Rstart, Xstart,Mstart = init_ξ_const(250,V,Nₗ,r₀,1e-16)
# Rstart .= r₀
# Xstart = Rstart.^3 .* 4 / 3 * π


sd_ml_v = calc_ρw(ξstart,Xstart,V)
sd_qᵥ = mv_v / (md_v + mv_v + sd_ml_v)
sd_qₗ = sd_ml_v / (md_v + mv_v + sd_ml_v)
sd_IC = ArrayPartition([Sₗ, p₀, T₀, sd_qᵥ, sd_qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, ln_INPC],ξstart, Rstart, Xstart)
# sd_IC = ArrayPartition([Sₗ, p₀, T₀, sd_qᵥ, sd_qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, ln_INPC],crystals)



# Simulation parameters passed into ODE solver
w = FT(10)                                 # updraft speed
const_dt = FT(0.5)                         # model timestep
t_max = FT(20)
size_distribution_list = ["Monodisperse"]#, "Gamma"]
condensation_growth = "Condensation"

# Data from Rogers(1975) Figure 1
# https://www.tandfonline.com/doi/abs/10.1080/00046973.1975.9648397
#! format: off
Rogers_time_supersat = [0.0645, 0.511, 0.883, 1.4, 2.07, 2.72, 3.24, 3.89, 4.53, 5.87, 7.16, 9.79, 16.0, 19.8]
Rogers_supersat = [0.0268, 0.255, 0.393, 0.546, 0.707, 0.805, 0.863, 0.905, 0.938, 0.971, 0.978, 0.963, 0.910, 0.885]
Rogers_time_radius = [0.561, 2, 3.99, 10.7, 14.9, 19.9]
Rogers_radius = [8.0, 8.08, 8.26, 8.91, 9.26, 9.68]
#! format: on

# Setup the plots
fig = MK.Figure(size = (800, 600))
ax1 = MK.Axis(fig[1, 1], ylabel = "Supersaturation [%]",title = "Ns = 250 superdroplets")
ax2 = MK.Axis(fig[3, 1], xlabel = "Time [s]", ylabel = "Temperature [K]")
ax3 = MK.Axis(fig[2, 1], ylabel = "q_vap [g/kg]")
ax4 = MK.Axis(fig[2, 2], xlabel = "Time [s]", ylabel = "q_liq [g/kg]")
ax5 = MK.Axis(fig[1, 2], ylabel = "radius [μm]")
ax6 = MK.Axis(fig[3, 2], xlabel = "Time [s]", ylabel = "air pressure")
MK.lines!(ax1, Rogers_time_supersat, Rogers_supersat, label = "Rogers_1975", color = :black)
MK.lines!(ax5, Rogers_time_radius, Rogers_radius, label = "Rogers_1975", color = :black)


DSD = "Monodisperse"
#for DSD in size_distribution_list
    #local
    params_ = parcel_params{FT}(
                size_distribution = DSD,
                condensation_growth = condensation_growth,
                const_dt = const_dt/20,
                w = w,
            )


            # pp = params_
            # FT = eltype(sd_IC.x[1][1])
        
            # info = "\nSize distribution: $(pp.size_distribution)\n"
            # if pp.size_distribution == "Monodisperse"
            #     distr = Monodisperse{FT}()
            # elseif pp.size_distribution == "Gamma"
            #     distr = Gamma{FT}()
            # else
            #     throw("Unrecognized size distribution")
            # end
        
            # info *= "Aerosol: $(chop( string(typeof(pp.aerosol)), head = 29, tail = 9))\n"
        
            # info *= "Deposition: $(pp.deposition)\n"
            # if pp.deposition == "None"
            #     dep_params = Empty{FT}()
            # elseif pp.deposition == "MohlerAF"
            #     dep_params = MohlerAF{FT}(pp.ips, pp.aerosol, pp.tps, pp.const_dt)
            # elseif pp.deposition == "MohlerRate"
            #     dep_params = MohlerRate{FT}(pp.ips, pp.aerosol, pp.tps, pp.const_dt)
            # elseif pp.deposition == "ABDINM"
            #     dep_params = ABDINM{FT}(pp.tps, pp.aerosol, pp.r_nuc, pp.const_dt)
            # elseif pp.deposition == "P3_dep"
            #     dep_params = P3_dep{FT}(pp.ips, pp.const_dt)
            # else
            #     throw("Unrecognized deposition mode")
            # end
        
            # info *= "Heterogeneous: $(pp.heterogeneous)\n"
            # if pp.heterogeneous == "None"
            #     imm_params = Empty{FT}()
            # elseif pp.heterogeneous == "ABIFM"
            #     imm_params = ABIFM{FT}(pp.tps, pp.aerosol, pp.A_aer, pp.const_dt)
            # elseif pp.heterogeneous == "P3_het"
            #     imm_params = P3_het{FT}(pp.ips, pp.const_dt)
            # elseif pp.heterogeneous == "Frostenberg_random"
            #     imm_params =
            #         Frostenberg_random{FT}(pp.ip, pp.sampling_interval, pp.const_dt)
            # elseif pp.heterogeneous == "Frostenberg_mean"
            #     imm_params = Frostenberg_mean{FT}(pp.ip, pp.const_dt)
            # elseif pp.heterogeneous == "Frostenberg_stochastic"
            #     imm_params = Frostenberg_stochastic{FT}(pp.ip, pp.γ, pp.const_dt)
            # else
            #     throw("Unrecognized heterogeneous mode")
            # end
        
            # info *= "Homogeneous: $(pp.homogeneous)\n"
            # if pp.homogeneous == "None"
            #     hom_params = Empty{FT}()
            # elseif pp.homogeneous == "ABHOM"
            #     hom_params = ABHOM{FT}(pp.tps, pp.ips, pp.const_dt)
            # elseif pp.homogeneous == "P3_hom"
            #     hom_params = P3_hom{FT}(pp.const_dt)
            # else
            #     throw("Unrecognized homogeneous mode")
            # end
        
            # info *= "Condensation growth: $(pp.condensation_growth)\n"
            # if pp.condensation_growth == "None"
            #     ce_params = Empty{FT}()
            # elseif pp.condensation_growth == "Condensation"
            #     ce_params = CondParams{FT}(pp.aps, pp.tps)
            # else
            #     throw("Unrecognized condensation growth mode")
            # end
        
            # info *= "Deposition growth: $(pp.deposition_growth)\n"
            # if pp.deposition_growth == "None"
            #     ds_params = Empty{FT}()
            # elseif pp.deposition_growth == "Deposition"
            #     ds_params = DepParams{FT}(pp.aps, pp.tps)
            # else
            #     throw("Unrecognized deposition growth mode")
            # end
            # @info info
        
            # # Parameters for the ODE solver
            # p = (
            #     distr = distr,
            #     dep_params = dep_params,
            #     imm_params = imm_params,
            #     hom_params = hom_params,
            #     ce_params = ce_params,
            #     ds_params = ds_params,
            #     wps = pp.wps,
            #     tps = pp.tps,
            #     r_nuc = pp.r_nuc,
            #     w = pp.w,
            # )
    params_m = parcel_params{FT}(
                size_distribution = DSD,
                condensation_growth = condensation_growth,
                const_dt = const_dt,
                w = w,
            )
    params_g = parcel_params{FT}(
                size_distribution = "Gamma",
                condensation_growth = condensation_growth,
                const_dt = const_dt,
                w = w,
            )           
    # solve ODE
    sol = sd_run_parcel(sd_IC, FT(0), FT(t_max), params_)
    also = run_parcel(IC,FT(0), FT(t_max), params_m)
    gamma = run_parcel(IC,FT(0), FT(t_max), params_g)
    # Initialize empty arrays for each variable
Sₗ_series = []
p₀_series = []
T₀_series = []
sd_qᵥ_series = []
sd_qₗ_series = []
qᵢ_series = []
Nₐ_series = []
Nₗ_series = []
Nᵢ_series = []
ln_INPC_series = []
ξstart_series = []
Rstart_series = []
Xstart_series = []

# Iterate over the solution
for state in sol.u
    # Unpack the ArrayPartition
    (Sₗ, p0, T0, sd_qv, sd_ql, qi, Na, Nl, Ni, lnINPC), ξs, R, X = state.x
    # Append the values to the corresponding arrays
    push!(Sₗ_series, Sₗ)
    push!(p₀_series, p0)
    push!(T₀_series, T0)
    push!(sd_qᵥ_series, sd_qv)
    push!(sd_qₗ_series, sd_ql)
    push!(qᵢ_series, qi)
    push!(Nₐ_series, Na)
    push!(Nₗ_series, Nl)
    push!(Nᵢ_series, Ni)
    push!(ln_INPC_series, lnINPC)
    push!(ξstart_series, ξs)
    push!(Rstart_series, R)
    push!(Xstart_series, X)
end

    # Plot results
    MK.lines!(ax1, sol.t, (Sₗ_series  .- 1) * 100.0, label = "Droplets.jl", color = :blue)
    MK.lines!(ax2, sol.t, T₀_series.+0, label = "Droplets.jl", color = :blue)
    MK.lines!(ax3, sol.t, sd_qᵥ_series * 1e3, label = "Droplets.jl", color = :blue)
    MK.lines!(ax4, sol.t, sd_qₗ_series * 1e3, label = "Droplets.jl", color = :blue)
    rₗ = [(3*mean(X)/4pi)^(1/3) for X in Xstart_series]
    MK.lines!(ax5, sol.t, rₗ * 1e6,label= "SD effective rad",color = :blue)
    MK.lines!(ax6, sol.t, p₀_series.+0, label = "Droplets.jl", color = :blue)


    MK.lines!(ax1, also.t, (also[1, :] .- 1) * 100.0, label = DSD, color = :red)
    MK.lines!(ax2, also.t, also[3, :], label = DSD, color = :red)
    MK.lines!(ax3, also.t, also[4, :] * 1e3, label = DSD, color = :red)
    MK.lines!(ax4, also.t, also[5, :] * 1e3, label = DSD, color = :red)
    MK.lines!(ax6, also.t, also[2, :], label = DSD, color = :red)
    sol_Nₗ = also[8, :]
    sol_Nᵢ = also[9, :]
    # Compute the current air density
    sol_T = also[3, :]
    sol_p = also[2, :]
    sol_qᵥ = also[4, :]
    sol_qₗ = also[5, :]
    sol_qᵢ = also[6, :]
    q = TD.PhasePartition.(sol_qᵥ + sol_qₗ + sol_qᵢ, sol_qₗ, sol_qᵢ)
    ts = TD.PhaseNonEquil_pTq.(tps, sol_p, sol_T, q)
    ρₐ = TD.air_density.(tps, ts)
    # Compute the mean particle size based on the distribution
    distr = also.prob.p.distr
    moms = distribution_moments.(distr, sol_qₗ, sol_Nₗ, ρₗ, ρₐ, sol_qᵢ, Nᵢ, ρᵢ)
    rₗ = similar(sol_T)
    for it in range(1, length(sol_T))
        rₗ[it] = moms[it].rₗ
    end
    MK.lines!(ax5, also.t, rₗ * 1e6, label = DSD, color = :red)

    MK.lines!(ax1, gamma.t, (gamma[1, :] .- 1) * 100.0, label = "Gamma", color = :green)
    MK.lines!(ax2, gamma.t, gamma[3, :], label = "Gamma", color = :green)
    MK.lines!(ax3, gamma.t, gamma[4, :] * 1e3, label = "Gamma", color = :green)
    MK.lines!(ax4, gamma.t, gamma[5, :] * 1e3, label = "Gamma", color = :green)
    sol_Nₗ = gamma[8, :]
    sol_Nᵢ = gamma[9, :]
    # Compute the current air density
    sol_T = gamma[3, :]
    sol_p = gamma[2, :]
    sol_qᵥ = gamma[4, :]
    sol_qₗ = gamma[5, :]
    sol_qᵢ = gamma[6, :]
    q = TD.PhasePartition.(sol_qᵥ + sol_qₗ + sol_qᵢ, sol_qₗ, sol_qᵢ)
    ts = TD.PhaseNonEquil_pTq.(tps, sol_p, sol_T, q)
    ρₐ = TD.air_density.(tps, ts)
    # Compute the mean particle size based on the distribution
    distr = gamma.prob.p.distr
    moms = distribution_moments.(distr, sol_qₗ, sol_Nₗ, ρₗ, ρₐ, sol_qᵢ, Nᵢ, ρᵢ)
    rₗ = similar(sol_T)
    for it in range(1, length(sol_T))
        rₗ[it] = moms[it].rₗ
    end
    MK.lines!(ax5, gamma.t, rₗ * 1e6, label = "Gamma", color = :green)
    MK.lines!(ax6, gamma.t, gamma[2, :], label = "Gamma", color = :green)


    display(fig)


MK.axislegend(
    ax1,
    framevisible = false,
    labelsize = 12,
    orientation = :horizontal,
    nbanks = 2,
    position = :rb,
)
MK.axislegend(
    ax2)
display(fig)
MK.save("liquid_only_parcel.png", fig)

