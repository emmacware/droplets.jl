# #Crystals!


# mutable struct Crystals{FT <: Real, T <: Union{Tuple{FT,FT}, Tuple{FT,FT,FT}}}
#     a::FT # equatorial radius (m)
#     c::FT # polar radius (m)
#     # r::FT # wet radius (m)
#     loc::T #location (m)
#     ξ::Int #multiplicity
#     M::FT # mass (g)
#     ρi::FT # ice density (kg/m^3)
#     ρd::FT # deposition density (kg/m^3)

# end

# function create_crystals(Ns,Nx,Ny,Δx,Δy,astart,cstart,ξstart,Mstart,ρistart,ρdstart)
#     crystals = []
#     for i in 1:Ns
#         a = astart[i]
#         c = cstart[i]
#         loc = [rand(0.003:Nx*Δx), rand(0.003:Ny*Δy)]
#         ξ = ξstart[i]
#         M = Mstart[i]
#         ρi = ρistart[i]
#         ρd = ρdstart[i]
#         crystal = Crystals(a,c,loc,ξ,M,ρi,ρd)
#         push!(crystals, crystal)
#     end
#     return crystals 
# end

# ################################################################################
# # Imagine a setup like this!!:

# Base.@kwdef struct deposition_params{FT} <: {FT} #not sure how to do the type here
#     deposition = "Capacitance" # "None", "Capacitance", "DepositionCoeff"
#     ρ_dep_p = "Constant" # "Constant", "Optimized"
#     constants = constants()
#     const_dt = 1
#     solver = "Euler" # Can we do this? Not sure how to pass
#     #etc...
# end

# function run_deposition!(deposition_params,crystals,statevar...)
#     if deposition_params.deposition == "Capacitance"
#         #do capacitance deposition
#         de_params = capacitance(constants,..)
#     elseif deposition_params.deposition == "DepositionCoeff"
#         #do deposition coeff
#         de_params = deposition_coeff(constants,..)
#     elseif deposition_params.deposition == "None"
#         de_params = empty()
#     end

#     if deposition_params.ρ_dep_p == "Constant"
#         ρdep_params = constant_ρdep(de_params,..)
#     elseif deposition_params.ρ_dep_p == "Optimized"
#         ρdep_params = optimized_ρdep(de_params,..)
#     end

#     # Parameters for the ODE solver
#     p = (
#         de_params = de_params,
#         ρdep_params = ρdep_params,
#         constants = deposition_params.constants,
#         #etc...
#     )

#     problem = ODE.ODEProblem(run_deposition, crystals_and_state, (t0,tend), p)
#     return ODE.solve(problem,ODE.Euler(),dt = deposition_params.const_dt)
# end


# function run_deposition(du,u,p,t)
#     de_params = p.de_params
#     ρdep_params = p.ρdep_params
#     constants = p.constants
    
#     crystals, statevar = u #unpacking, would look different obviously

#     du[crystals] = deposition_tendency(de_params,ρdep_params,constants,crystals,statevar)
    
#     dqi = sum(crystals.ξ.*crystals.M)

#     du[qi] = -dqi
#     du[T] = ...
#     du[P] = ...

#     return du
# end





