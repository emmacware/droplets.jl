#Crystals!
export create_crystals, create_statevar, ice_deposition

mutable struct Crystals{FT <: Real, T <: Union{Tuple{FT,FT}, Tuple{FT,FT,FT}}}
    a::FT # equatorial radius (m)
    c::FT # polar radius (m)
    # r::FT # wet radius (m)
    loc::T # location (m)
    vel::T # velocity (m/s)
    ξ::FT # multiplicity can be a float why not
    M::FT # mass (g)
    ρi::FT # ice density (kg/m^3)
    ρd::FT # deposition density (kg/m^3)



end

function create_crystals(Ns,Nx,Ny,Δx,Δy,astart,cstart,ξstart,Mstart,ρistart,ρdstart)
    crystals = []
    for i in 1:Ns
        a = astart[i]
        c = cstart[i]
        loc = (rand(0.003:Nx*Δx), rand(0.003:Ny*Δy))
        vel = (0.0, 0.0)
        ξ = ξstart[i]
        M = Mstart[i]
        ρi = ρistart[i]
        ρd = ρdstart[i]
        crystal = Crystals(a,c,loc,ξ,M,ρi,ρd)
        push!(crystals, crystal)
    end
    return crystals 
end

function create_crystals(Ns, Nx, Ny, Nz, Δx, Δy, Δz, astart, cstart, ξstart, Mstart, ρistart, ρdstart)
    # Another method in 3D
    crystals = []
    for i in 1:Ns
        a = astart[i]
        c = cstart[i]
        loc = (rand(0.003:Nx*Δx), rand(0.003:Ny*Δy), rand(0.003:Nz*Δz))
        vel = (0.0, 0.0, 0.0)
        ξ = ξstart[i]
        M = Mstart[i]
        ρi = ρistart[i]
        ρd = ρdstart[i]
        crystal = Crystals(a, c, loc, ξ, M, ρi, ρd)
        push!(crystals, crystal)
    end
    return crystals 
end

################################################################################
# Imagine a setup like this!!:

Base.@kwdef struct deposition_params{FT} 
    deposition = "Capacitance" # "None", "Capacitance", "DepositionCoeff"
    ρ_dep_p = "Constant" # "Constant", "Optimized"
    constants = constants()
    const_dt = 1
    solver = "Euler" # Can we do this? Not sure how to pass
    #etc...
end

function run_deposition!(deposition_params,crystals,statevar...)
    if deposition_params.deposition == "Capacitance"
        #do capacitance deposition
        de_params = capacitance(constants,..)
    elseif deposition_params.deposition == "DepositionCoeff"
        #do deposition coeff
        de_params = deposition_coeff(constants,..)
    elseif deposition_params.deposition == "None"
        de_params = empty()
    end

    if deposition_params.ρ_dep_p == "Constant"
        ρdep_params = constant_ρdep(de_params,..)
    elseif deposition_params.ρ_dep_p == "Optimized"
        ρdep_params = optimized_ρdep(de_params,..)
    end

    # Parameters for the ODE solver
    p = (
        de_params = de_params,
        ρdep_params = ρdep_params,
        constants = deposition_params.constants,
        #etc...
    )

    problem = ODE.ODEProblem(run_deposition, crystals_and_state, (t0,tend), p)
    return ODE.solve(problem,ODE.Euler(),dt = deposition_params.const_dt)
end


function run_deposition(du,u,p,t)
    de_params = p.de_params
    ρdep_params = p.ρdep_params
    constants = p.constants
    
    crystals, qi = u #unpacking, would look different obviously

    du[crystals] = deposition_tendency(de_params,ρdep_params,constants,crystals,statevar)
    
    dqi = sum(crystals.ξ.*crystals.M)

    du[qi] = -dqi

    # du[T] = ...
    # du[P] = ...

    return du
end

function deposition_tendency(de_params, ρdep_params, constants, crystals, statevar)
    du = []
    for crystal in crystals
        # get the r^2 growth tendency

        # get equivalent radius, fortran code is ((sd_ra**2)*sd_rc)**0.3333333d0
        r = ((crystal.a^2)*crystal.c)^(1/3)

        # get capacitance, for now, assume capacitance = a (approximation around ϕ = 1)
        capacitance = crystal.a

        afn = 0.1 # relative growth rate, should be Gi*single

        dr2 = afn*(capacitance/r)/crystal.ρd

        push!(du, dr2)
    end
    return du
end



# For now, just create a function ice_deposition that takes in a single crystal and statevar, and returns the crystal with the updated radius
function ice_deposition(crystal, Δt, statevar)

    r = ((crystal.a^2)*crystal.c)^(1/3)

    # get capacitance, for now, assume capacitance = a (approximation around ϕ ≈ 1)
    capacitance = crystal.a

    afn = calculate_afn(crystal, statevar)
    println("how about now")

    r = (r^2 + Δt*afn*(capacitance/r)/crystal.ρd)^0.5

    return r
end

mutable struct statevar
    # 4 scalars: Temp, Pressure, qv, rhom_scale
    # Make these all keyword arguments
    T::Float64
    P::Float64
    qv::Float64
    rhom_scale::Float64
end

# function to create statevar object with keyword arguments
function create_statevar(;T = 273.0, P = 101325, qv = 1.0, rhom_scale = 1.695)
    return statevar(T, P, qv, rhom_scale)
end

# function to make afn
function calculate_afn(crystal, statevar)

    # Temperature variables
    t0 = 273.0 # K 
    tdeg_sd = statevar.T - t0

    # get maximum dimension of the crystal
    max_diameter = 2.0 * max(crystal.a, crystal.c)


    # Diffusion coefficent [m^2/s] (Beard and Pruppacher, 1971)
    Diff_C_new = 0.211 * ((statevar.T / t0) ^ 1.94) * (101325 / statevar.P) * 1.0E-4

    # Calculate X = Nsc^(1/3)*Nre^(1/2) for ventilation
    # dynamic viscosity μ [kg/(m*s)] at the location of super-droplets
    # (Pruppacher & Klett,1997)
    if tdeg_sd >= 0.0
        μ = (1.7180 + 4.9e-3 * tdeg_sd) * 1.0e-5
    else
        μ = (1.7180 + 4.9e-3 * tdeg_sd - 1.2e-5 * tdeg_sd^2) * 1.0e-5
    end
    Re = statevar.rhom_scale * max_diameter * crystal.vel[end] / μ
    Sc = μ / (statevar.rhom_scale * Diff_C_new)
    Xvnt = Sc^(1/3)*Re^(1/2)
    


    afn = 0.0 # relative growth rate, should be Gi*si

    return afn
end





