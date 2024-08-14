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
    ρ::FT # crystal ice density (kg/m^3)
    ρdep::FT # deposition density (kg/m^3)



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
    # 4 scalars: Temp, Pressure, qv, ρ
    # Make these all keyword arguments
    T::Float64 #Kelvin
    P::Float64 #Pascal
    qv::Float64
    ρ::Float64 # density of moist air
    esati::Float64 # saturation vapor pressure over ice
    ssi::Float64 # super saturation ice
    qvsi::Float64 # sat mixing ratio over ice
end

# function to create statevar object with keyword arguments
function create_statevar(;T = 273.0, P = 101325, qv = 1.0, rhom_scale = 1.695)
    return statevar(T, P, qv, rhom_scale)
end

# function to make afn
# function calculate_afn(crystal, statevar)

#     # Temperature variables
#     t0 = 273.0 # K 
#     tdeg_sd = statevar.T - t0

#     # get maximum dimension of the crystal
#     max_diameter = 2.0 * max(crystal.a, crystal.c)


#     # Diffusion coefficent [m^2/s] (Beard and Pruppacher, 1971)
#     Diff_C_new = 0.211 * ((statevar.T / t0) ^ 1.94) * (101325 / statevar.P) * 1.0E-4

#     # Calculate X = Nsc^(1/3)*Nre^(1/2) for ventilation
#     # dynamic viscosity μ [kg/(m*s)] at the location of super-droplets
#     # (Pruppacher & Klett,1997)
#     if tdeg_sd >= 0.0
#         μ = (1.7180 + 4.9e-3 * tdeg_sd) * 1.0e-5
#     else
#         μ = (1.7180 + 4.9e-3 * tdeg_sd - 1.2e-5 * tdeg_sd^2) * 1.0e-5
#     end
#     Re = statevar.rhom_scale * max_diameter * crystal.vel[end] / μ
#     Sc = μ / (statevar.rhom_scale * Diff_C_new)
#     Xvnt = Sc^(1/3)*Re^(1/2)
    


#     afn = 0.0 # relative growth rate, should be Gi*si

#     return afn
# end




##########################################################
##########################################################
function termvz_ice_Heymsfield_Westbrook_2010(equatorial_radius,  polar_radius,  crystal_density, 
    constants, statevar)
    """
    Calculate the terminal velocity of ice crystals based on the Heymsfield and Westbrook (2010) parameterization.
    # Arguments
    - `equatorial_radius::Float64`: Equatorial radius of the ice crystal [m]
    - `polar_radius::Float64`: Polar radius of the ice crystal [m]
    - `crystal_density::Float64`: Density of the ice crystal [kg/m^3]
    - `constants::Struct`: Mutable struct containing physical constants
    - `statevar::Struct`: Mutable struct containing local environmental variables
    # Returns
    - `sd_vz::Float64`: Terminal velocity of the ice crystal [m/s]
    # Reference
    Heymsfield, A. J., & Westbrook, C. D. (2010). Advances in the estimation 
    of ice particle fall speeds using laboratory and field measurements. 
    Journal of the Atmospheric Sciences, 67(8), 2469-2482. 
    https://doi.org/10.1175/2010JAS3379.1
    """
    # https://doi.org/10.1175/2010JAS3379.1

    # Local parameters
    C0 = 0.35        # Fitting parameter of Re(X*)
    delta0 = 8.0     # Fitting parameter of Re(X*)
    var_k_coef = 1.0 # parameter to determine the projected area of ice particles, from SCALE-SDM
    # maximum dimension of the particle [m]
    sd_maxD = 2.0 * max(equatorial_radius, polar_radius)

    # dynamic viscosity [kg/(m*s)] at the location of super-droplets
    tdeg = statevar.T - constants.T0  # [K] => [degC]
    if tdeg >= 0.0
    sd_dvisc = (1.7180 + 4.9E-3 * tdeg) * 1E-5
    else
    sd_dvisc = (1.718 + 4.9E-3 * tdeg - 1.2E-5 * tdeg^2) * 1E-5
    end

    # mass of the particle [kg]
    sd_mass = (4 / 3) * π * equatorial_radius^2 * polar_radius * crystal_density

    # aspect ratio
    ϕ = polar_radius / equatorial_radius

    # A: area of the particle projected to the flow direction [m^2]
    var_k = exp(-var_k_coef * ϕ)
    sd_area = π * equatorial_radius * (sd_maxD / 2.0) * (crystal_density / constants.ρ_ice)^var_k

    # Ar: A divided by the area of a circumscribing disk
    sd_area_ratio = sd_area / (π * sd_maxD^2 / 4.0)

    # modified Best (Davies) number of the particle
    n_modX = (crystal_density / sd_dvisc^2) * (8.0 * sd_mass * constants.gconst) / π / sqrt(sd_area_ratio)

    # Reynolds number of the particle
    nre = (delta0^2 / 4.0) * (sqrt(1.0 + 4.0 * sqrt(n_modX / C0) / (delta0^2)) - 1.0)^2

    # terminal velocity [m/s]
    sd_vz = sd_dvisc * max(nre, 0.0) / crystal_density / sd_maxD

    return sd_vz
end

# function Gi(statevar,constants)
#     first_term = constants.R*statevar.T/(constants.Mw*Dprimev()*sati(statevar.T))
#     second_term = (Lsub(statevar.T))/(constants.Mw*kprimeT()*statevar.T)
#     third_term = (Lsub(statevar.T))/(constants.R*statevar.T)-1

#     return 1/(first_term + second_term*third_term)
# end

# function Dprimev()

#     return 
# end

# function kprimeT()

#     return 
# end

function fdvstr(statevar,constants,alen,clen)
    2/3*1/(lc_a(statevar.P,statevar.T)/caplengtha(alen,clen) + capacitancedep(alen,clen)/capnew(statevar.P,statevar.T,alen,clen)) +
    1/3*1/(lc_c(statevar.P,statevar.T)/caplengthc(alen,clen) + capacitancedep(alen,clen)/capnew(statevar.P,statevar.T,alen,clen))
end

function fktster(statevar,constants,alen,clen)
    2/3*1/(lt(statevar)/caplengtha(alen,clen) +capacitancedep(alen,clen)/capnew(statevar.P,statevar.T,alen,clen)) +
    1/3*1/(lt(statevar)/caplengthc(alen,clen) + capacitancedep(alen,clen)/capnew(statevar.P,statevar.T,alen,clen))
end

lc_a(P,T) = Dv(P,T)/max(alpha_a,1e-35)*sqrt(2pi/(constants.Rv*T))
lc_c(P,T) = DiffCoeff(P,T)/max(alpha_c,1e-35)*sqrt(2pi/(constants.Rv*T))
lt(statevar) = HeatCoeff(statevar.T,statevar.qv)/(alphat*statevar.ρ*constants.Cp_dry)*sqrt(2pi/(constants.Rd*statevar.T))

#Beard and Pruppacher 1971, m2/s
Dv(P,T) = 2.11e-5*(T/constants.T0)^1.94*(constants.P0/P)
HeatDry(T) = (5.69+1.68e-2*(T-constants.T0))*1e-5 #Thermal Conductivity of dry air (cal/cm/s/C)
HeatVap(T) = (3.73+2.00e-2*(T-constants.T0))*1e-5 #Thermal Conductivity of water vapor (cal/cm/s/C)
HeatCoeff(T,qv) = 4.184e2 * HeatDry(T)*(1-qv*(1.17-1.02*(HeatVap(T)/HeatDry(T))))#Thermal conductivity of moist air (J/(msK))

alpha_a = 1 #change
alpha_c = 1 #change
alphat = 1 #also says 0.96???

caplengtha(alen,clen) = alen*clen/capacitancedep(alen,clen) # =capa
caplengthc(alen,clen) = alen^2/capacitancedep(alen,clen) # =capc 
# cap = capacitancedep(alen,clen)
delv(P,T) = 6.6e-8*(P/constants.P0)*(T/(constants.T0+20)) #mean free path of ai 
# amn = alen+delv
# cmn = clen+delv
function capnew(P,T,alen,clen)
    amn = alen+delv(P,T)
    cmn = clen+delv(P,T)
    return capacitancedep(amn,cmn)
end

function capacitancedep(a,c)
    phi = c/a
    if phi < 1
        return (a^2-c^2)*0.5/asin((1-phi^2)^0.5)
    elseif phi > 1
        return (c^2-a^2)*0.5/log(phi+phi*(1-phi^(-2))^0.5)
    elseif phi == 1
        return a
    end
end

#Enthalpy
Lv(T) = 2.5008e6*(constants.T0/T)^(0.167+3.67e-4*T)
Lf(T) = 0.3337e6+(2031-10.467*(T-constants.T0))*(T-constants.T0)
Lsub(T) = Lv(T)-Lf(T)

function afn(statevar,constants,fhdum,alen,clen)
    #calculate barf here? or needed elsewhere?
    #calculate fhdum here? instead of passing?
    afn = Dv(statevar.P,statevar.T)*fdvstr(statevar,constants,alen,clen)*barf()*statevar.esati/constants.Rv*statevar.T *
    (statevar.ssi-del(statevar,constants,fhdum,alen,clen)*(Lsub(0)/(constants.Rv*statevar.T))-1)
    return afn^(-1)
end

function del(statevar,constants,fhdum,alen,clen)
    del = statevar.ssi*(((statevar.T/alpha(statevar,fhdum,alen,clen))+((Lsub(0)/(constants.Rv*statevar.T))-1))^(-1))
    return del
end

function alpha(statevar,fhdum,alen,clen)
    numerator = Dv(statevar.P,statevar.T)*fdvstr(statevar,constants,alen,clen)*barf()*statevar.esati*Lsub(0)
    denominator = constants.Rv*HeatCoeff(statevar.T,statevar.qv)*fktster(statevar,constants,alen,clen)*fhdum*statevar.T
    return numerator/denominator
end


#dynamic viscosity [kg/(m*s)]
#(Pruppacher & Klett,1997)
if statvar.T >= constants.T0
    dvisc= ( 1.7180 + 4.9e-3*(statevar.T-constants.T0) ) * 1e-5
else
    dvisc = ( 1.718 + 4.9e-3*(statevar.T-constants.T0)-
        - 1.2e-5*(statevar.T-constants.T0)^2) * 1e-5
end


###########################################


#routineeee for crystals

a = crystal.a
c = crystal.c
phi = c/a
max_diameter = 2.0 * max(a, c)
radius = ((a^2)*c)^(1/3)


# capacitance ######!
if sd_phi > 1.001 #prolate
    ecent = sqrt(1 - 1/(phi^2))
    capaci = ecent*c / log((1+ecent)*phi)            
elseif sd_phi < 0.999 #oblate     
    ecent = sqrt(1 - phi^2)
    capaci = ecent*a / asin(ecent)                  
else #round
    capaci = sd_ra
end

igr = 1 #change inherant growth rate if needed
fs = capaci/radius
alphanr = a/radius^(3/(2+igr))  

cptot = constants.Cp_dry * (1.0-statevar.qv)
cptot = cptot + statevar.qv * constants.Cp_vap

sd_reynolds = statevar.ρ*max_diameter*crystal.vel/dvisc
sd_schmidt = dvisc / statevar.ρ / Dv(statevar.P,statevar.T)
sd_prandtl = dvisc_sd / (HeatCoeff(statevar.T,statevar.qv)/cptot)

n_ventX = sd_schmidt^(1/3) * sqrt(sd_reynolds)

# kinetic correction
ivt = 1 / statevar.T
Fac_ice_kk_invrhoi = ( ivt * Lsub(0) / constants.Rv - 1 ) * (Lsub(0)/HeatCoeff(statevar.T,statevar.qv)) * ivt
# psychrometric correction
Fac_ice_dd_invrhoi = constants.Rv / Dv(statevar.P,statevar.T) * statevar.T / statevar.esati
# ventilation effect
if n_ventX <= 1.0
    barf = 1 + 0.14*n_ventX^2  
else
    barf = 0.86 + 0.28*n_ventX  
end
 
ntherm = sqrt(sd_reynolds)*sd_prandtl^(1/3)

if ntherm < 1.4
    fhdum = 1.0 + 0.108*ntherm^2
else
    fhdum = 0.78d0 + 0.308*ntherm
end

