#Crystals!
export create_crystals, create_statevar, ice_deposition

mutable struct Supercrystal{FT <: Real}#, T <:Tuple{FT,FT}}
    a::FT # equatorial radius (m)
    c::FT # polar radius (m)
    # r::FT # wet radius (m)
    xloc::FT # location (m)
    yloc::FT
    vel::FT # terminal velocity (m/s)
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
        xloc = rand(0.003:Nx*Δx)
        yloc = rand(0.003:Ny*Δy)
        vel = 0.0
        ξ = ξstart[i]
        M = Mstart[i]
        ρi = ρistart[i]
        ρd = ρdstart[i]
        crystal = Crystals(a,c,xloc,yloc,vel,ξ,M,ρi,ρd)
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
        crystal = Crystals(a, c, loc,vel, ξ, M, ρi, ρd)
        push!(crystals, crystal)
    end
    return crystals 
end

################################################################################
# Imagine a setup like this!!:

Base.@kwdef struct deposition_params{FT} 
    deposition = "Capacitance" # "None", "Capacitance", "DepositionCoeff"
    ρ_dep_p = "Constant" # "Constant", "Calculated"
    constants = constants # constants struct
    const_dt = 1
    solver = "Euler" # Can we do this? Not sure how to pass
    mval = 10.0
    ro = 1e-6
    temppoly = 268.15 # temperature for polycrystal approximation
    #etc...
end
#deposition params: ro (smallest polycrystal radius),mval??? calculate/constant rho_dep

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

mutable struct statevar_struct
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

# # function to create statevar object with keyword arguments
# function create_statevar(T,P,qv,ρ,esati,ssi,qvsi)
#     return statevar(T, P, qv, ρ, esati, ssi, qvsi)
# end

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
vol(a,c) = 4/3 * π * a^2 * c 

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




function fdvstr(statevar,constants,alen,clen,alpha_a,alpha_c)
    lc_a = Dv(statevar.P,statevar.T)/max(alpha_a,1e-35)*sqrt(2pi/(constants.Rv*statevar.T))
    lc_c = Dv(statevar.P,statevar.T)/max(alpha_c,1e-35)*sqrt(2pi/(constants.Rv*statevar.T))
    fdv = 2/3*1/(lc_a/caplengtha(alen,clen) + capacitancedep(alen,clen)/capnew(statevar.P,statevar.T,alen,clen)) +
    1/3*1/(lc_c/caplengthc(alen,clen) + capacitancedep(alen,clen)/capnew(statevar.P,statevar.T,alen,clen))
    return fdv
end

function fktster(statevar,constants,alen,clen,alphat)
    lt = HeatCoeff(statevar.T,statevar.qv)/(alphat*statevar.ρ*constants.Cp_dry)*sqrt(2pi/(constants.Rd*statevar.T))
    2/3*1/(lt/caplengtha(alen,clen) +capacitancedep(alen,clen)/capnew(statevar.P,statevar.T,alen,clen)) +
    1/3*1/(lt/caplengthc(alen,clen) + capacitancedep(alen,clen)/capnew(statevar.P,statevar.T,alen,clen))
end


#Beard and Pruppacher 1971, m2/s
Dv(P,T) = 2.11e-5*(T/constants.T0)^1.94*(constants.P0/P)
HeatDry(T) = (5.69+1.68e-2*(T-constants.T0))*1e-5 #Thermal Conductivity of dry air (cal/cm/s/C)
HeatVap(T) = (3.73+2.00e-2*(T-constants.T0))*1e-5 #Thermal Conductivity of water vapor (cal/cm/s/C)
HeatCoeff(T,qv) = 4.184e2 * HeatDry(T)*(1-qv*(1.17-1.02*(HeatVap(T)/HeatDry(T))))#Thermal conductivity of moist air (J/(msK))



caplengtha(alen,clen) = alen*clen/capacitancedep(alen,clen) # =capa
caplengthc(alen,clen) = alen^2/capacitancedep(alen,clen) # =capc 
delv(P,T) = 6.6e-8*(P/constants.P0)*(T/(constants.T0+20)) #mean free path of ai 

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
Lv(T,T0) = 2.5008e6*(T0/T)^(0.167+3.67e-4*T)
Lf(T,T0) = 0.3337e6+(2031-10.467*(T-T0))*(T-T0)
Lsub(T,T0) = Lv(T,T0)-Lf(T,T0)

function calc_afn(statevar,constants,fhdum,barf,alen,clen,alpha_a,alpha_c,alphat)
    alpha_num = Dv(statevar.P,statevar.T)*fdvstr(statevar,constants,alen,clen,alpha_a,alpha_c)*barf*statevar.esati*Lsub(constants.T0,constants.T0)
    alpha_denom = constants.Rv*HeatCoeff(statevar.T,statevar.qv)*fktster(statevar,constants,alen,clen,alphat)*fhdum*statevar.T

    del = statevar.ssi*(((statevar.T/(alpha_num/alpha_denom))+((Lsub(constants.T0,constants.T0)/(constants.Rv*statevar.T))-1))^(-1))

    afn = Dv(statevar.P,statevar.T)*fdvstr(statevar,constants,alen,clen,alpha_a,alpha_c)*barf*statevar.esati/constants.Rv*statevar.T *
    (statevar.ssi-del*(Lsub(constants.T0,constants.T0)/(constants.Rv*statevar.T))-1)
    return afn^(-1)
end




function getscritdep(temp_celsius::Float64)
    """
    This function finds the critical (characteristic) supersaturations
    for the basal (c-axis) and prism (a-axis) facets based on polynomial
    fits to available data.

    JYH Penn State University, May 15, 2019

    Article: Harrington et al. (2019, J. Atmos. Sci.)

    Parameters:
    -----------
    - `temp_celsius`: Temperature in Celsius
    Returns:
    --------
    - `sa`: Critical supersaturation for the prism (a-axis) facet
    - `sc`: Critical supersaturation for the basal (c-axis) facet
    """
    x = max(-70.0, min(temp_celsius, -1.0))
    sc = 0.0
    sa = 0.0
    
    if temp_celsius < -30.0
        # original magee (high scrit) 
        sc = 1.8115 + 0.15585 * x + 0.011569 * x^2
        sa = 1.8115 + 0.15585 * x + 0.011569 * x^2
        
        # with bailey and hallett
        sc = 3.7955 + 0.10614 * x + 0.0075309 * x^2
        sa = sc
        # modified magee
        # sc = 3.3649 + 0.13266 * x + 0.008959 * x^2
        # sa = sc
    elseif temp_celsius >= -30.0 && temp_celsius <= -22.0
        sc = 753.63 + 105.97 * x + 5.5532 * x^2 + 0.12809 * x^3 +
             0.001103 * x^4
    elseif temp_celsius > -22.0 && temp_celsius <= -1.0
        sc = 1.1217 + 0.038098 * x - 0.083749 * x^2 -
             -0.015734 * x^3 - 0.0010108 * x^4 -
             -2.9148e-05 * x^5 - 3.1823e-07 * x^6
    end
    
    if temp_celsius >= -30.0 && temp_celsius <= -22.0
        sa = -0.71057 - 0.14775 * x + 0.0042304 * x^2
    elseif temp_celsius > -22.0 && temp_celsius <= -15.0
        sa = -5.2367 - 1.3184 * x - 0.11066 * x^2 - 0.0032303 * x^3
    elseif temp_celsius > -15.0 && temp_celsius <= -10.0
        # sa = 0.755 + 0.03325 * x + 0.00125 * x^2   # fit to Woods
        sa = 0.34572 - 0.0093029 * x + 0.00030832 * x^2 # fit to Nelson & Knight
    elseif temp_celsius > -10.0 && temp_celsius <= -1.0
        # sa = 0.37445 + 0.029636 * x + 0.011575 * x^2 +  # Wood
        #      0.00069444 * x^3
        sa = 0.34572 - 0.0093029 * x + 0.00030832 * x^2 # Nelson & Knight
    end
    
    sa /= 100.0
    sc /= 100.0
    
    return sa, sc
end


function adjustscritdep(temp_celsius::Float64, ihabit::Int, iequiv::Int, 
                        sc_temp::Float64, sa_temp::Float64)
    """
    This function adjusts the critical supersaturations based on temperature
    and habit type.

    Parameters:
    - temp_celsius: Temperature in Celsius
    - ihabit: Habit type (1, 2, or 3)
    - iequiv: Equivalence flag (1 or other)
    - sc_temp: Supersaturation for the c-axis
    - sa_temp: Supersaturation for the a-axis

    Returns:
    - scrit_c: Adjusted supersaturation for the c-axis
    - scrit_a: Adjusted supersaturation for the a-axis
    """
    scrit_c = sc_temp
    scrit_a = sa_temp

    if temp_celsius <= -40.0
        if ihabit == 1
            scrit_c = sc_temp / 1.5
            # print("making columns ", scrit_c, scrit_a)
        elseif ihabit == 2
            scrit_a = sa_temp / 1.5
        end
    end

    if ihabit == 3
        scrit_a = scrit_c
    end

    if iequiv == 1
        scrit_a = 2.0 / 3.0 * scrit_a + 1.0 / 3.0 * scrit_c
    end

    return scrit_c, scrit_a
end


struct FittingParams
    """
    This struct holds the fitting parameters used in the fits for the
    surface supersaturation based on the chosen value of M.
    """
    slopelg::Float64
    aresid::Float64
    cresid::Float64
    uresid::Float64
    nup1::Float64
    nup2::Float64
    ndwn1::Float64
    ndwn2::Float64
end

function getfittingparams(m::Float64)
    """
    This function calculates the fitting parameters used in the fits for the
    surface supersaturation based on the chosen value of M.

    Parameters:
    - m: Parameter in Nelson and Baker (1996) approximation for the deposition coefficient

    Returns:
    - FittingParams: Struct containing the fitting parameters
    """
    slopelg = 0.832
    slopelg = 0.1532 + 0.49078 * m - 0.18002 * m^2 + 0.04165 * m^3 - 
              0.0061541 * m^4 + 0.00058009 * m^5 - 3.3791e-5 * m^6 + 
              1.1092e-6 * m^7 - 1.5691e-8 * m^8
    aresid = 0.46192 - 0.093575 * m + 0.012798 * m^2 - 
             0.0008756 * m^3 + 2.2466e-5 * m^4
    cresid = aresid
    uresid = 0.25872 - 0.057957 * m + 0.0079602 * m^2 - 
             0.00052905 * m^3 + 1.3145e-5 * m^4
    nup1 = -0.69061 + 10.765 * m - 9.1666 * m^2 + 3.8888 * m^3 - 
           0.91779 * m^4 + 0.12544 * m^5 - 0.0098291 * m^6 + 
           0.00040868 * m^7 - 6.9744e-6 * m^8
    nup2 = 5.1239 - 19.204 * m + 15.236 * m^2 - 6.4802 * m^3 + 
           1.5151 * m^4 - 0.20317 * m^5 + 0.015601 * m^6 - 
           0.00063752 * m^7 + 1.0736e-5 * m^8
    ndwn1 = -0.86527 - 0.03202 * m + 0.027342 * m^2 - 0.04378 * m^3 + 
            0.014734 * m^4 - 0.0021036 * m^5 + 0.0001431 * m^6 - 
            4.317e-6 * m^7 + 3.9242e-8 * m^8
    ndwn2 = -4.692 + 12.47 * m - 9.8936 * m^2 + 4.1881 * m^3 - 
            1.0003 * m^4 + 0.14011 * m^5 - 0.011359 * m^6 + 
            0.00049122 * m^7 - 8.7303e-6 * m^8

    return FittingParams(slopelg, aresid, cresid, uresid, nup1, nup2, ndwn1, ndwn2)
end

function calculate_alpha_ac(alen, clen, caplength,
                            statevar, constants,
                            scrit, m)
    """
    This function calculates the deposition coefficients for spheroids. 
    The theory follows that of Zhang and Harrington (2014), but the iteration 
    to find the alpha values is replaced by a parameterization. This version 
    allows for a choice of m.

    Article: Harrington et al. (2021, J. Atmos. Sci.)

    Parameters:
    -----------
    - `alen`: Spheroid equatorial radius (units: m)
    - `clen`: Spheroid polar radius (units: m)
    - `caplength`: Capacitance length (units: m)
    - `statevar`: Mutable struct containing state variables such as temperature (`T`), pressure (`P`), and fractional ice supersaturation (`ssi`)
    - `constants`: Mutable struct containing constants such as universal gas constant (`R`) and molecular weight (`Mw`)
    - `ei0`: Ice saturation vapor pressure
    - `m`: Parameter in Nelson and Baker (1996) approximation for the deposition coefficient
    - 
    - `scrit`: Critical supersaturation
    Returns:
    --------
    - `alpha`: Deposition coefficient
    """
    ei0 = eimk(statevar.T)
    # Fitting parameters
    fitting_params = getfittingparams(m)
    slopelg = fitting_params.slopelg
    aresid = fitting_params.aresid
    # cresid = fitting_params.cresid
    uresid = fitting_params.uresid
    nup1 = fitting_params.nup1
    nup2 = fitting_params.nup2
    ndwn1 = fitting_params.ndwn1
    ndwn2 = fitting_params.ndwn2

    xk = 2.3823e-2 + 7.1177e-5 * (statevar.T - constants.T0) # Thermal conductivity in W/m/K
    # Effective diffusivity
    gtp1 = constants.Rv * statevar.T / (Dv(statevar.P, statevar.T) * ei0)
    gtp2 = (Lsub(statevar.T,constants.T0)^2 / (xk * constants.Rv * statevar.T^2)) - (Lsub(statevar.T,constants.T0) / (xk * statevar.T))
    gtp_stand = 1.0 / (gtp1 + gtp2)
    # mfp = 6.6e-8 * (1013.25 / statevar.P) * (statevar.T / constants.T0) # Mean free path in m
    vel_molec = sqrt(8.0*constants.R*statevar.T/(π*constants.Mw*1e-3))  # molecular speed(m/s), change constants.Mw to kg/mol

    C3 = (vel_molec*caplength*capacitancedep(alen,clen)/capnew(statevar.T,statevar.P,alen,clen))/
         (4.0*Dv(statevar.P,statevar.T) )

    sisloc = Dv(statevar.P,statevar.T) / (gtp_stand * constants.Rv * statevar.T / ei0) / (1.0 / (1.0 + C3))
    slocal_diff = statevar.ssi / sisloc
    # scratio = min(max(slocal_diff / scrit, scrat_low), scrat_hi)
    slrat_interp = 10.0^(slopelg * log10(slocal_diff / scrit))
    slrat_interp_lim = max(min(slrat_interp, 1.0), 1.0 / sisloc)
    sldiff_upratio = max(0.0, slrat_interp - 1.0) + 1.0
    slrat_interp = max(slrat_interp, 1.0 / sisloc) +
                   aresid * slrat_interp * (slrat_interp / (1.0 / sisloc))^ndwn1 *
                   min((slrat_interp / (1.0 / sisloc))^ndwn2, 1.0) -
                   uresid * slrat_interp * slrat_interp_lim^nup1 *
                   sldiff_upratio^nup2
    slocal_p = slocal_diff / slrat_interp
    alpha = (slocal_p / scrit)^m * tanh((scrit / slocal_p)^m)
    alpha = max(min(alpha, 1.0), 1e-6)
    return alpha
end

function eimk(tmp::Float64)::Float64
    """
    Murphy and Koop (2005) parameterization.
    Parameters:
    - tmp: Temperature in Kelvin
    Returns:
    - eimk
    """
    ci = [9.550426, 5723.265, 3.53068, 0.00728332]
    alogpice = ci[1] - ci[2]/tmp + ci[3]*log(tmp) - ci[4]*tmp
    return exp(alogpice)
end



















###########################################














function crystal_deposition(statevar,crystal,constants,dep_params)
#routineeee for crystals -- will this be a function or a tendency???
# function crystal_tendency(crystal,dep_params,statevar,constants)
a = crystal.a
c = crystal.c
phi = c/a
max_diameter = 2.0 * max(a, c)
radius = ((a^2)*c)^(1/3)


# capacitance ######!
if phi > 1.001 #prolate
    ecent = sqrt(1 - 1/(phi^2))
    capaci = ecent*c / log((1+ecent)*phi)            
elseif phi < 0.999 #oblate     
    ecent = sqrt(1 - phi^2)
    capaci = ecent*a / asin(ecent)                  
else #round
    capaci = a
end

igr = 1 #put somewhere else?
fs = capaci/radius
alphanr = a/radius^(3/(2+igr))  


cptot = constants.Cp_dry * (1.0-statevar.qv)
cptot = cptot + statevar.qv * constants.Cp_vap

#dynamic viscosity [kg/(m*s)]
#(Pruppacher & Klett,1997)
if statevar.T >= constants.T0
    dvisc= ( 1.7180 + 4.9e-3*(statevar.T-constants.T0) ) * 1e-5
else
    dvisc = ( 1.718 + 4.9e-3*(statevar.T-constants.T0)-
        1.2e-5*(statevar.T-constants.T0)^2) * 1e-5
end


sd_reynolds = statevar.ρ*max_diameter*crystal.vel/dvisc
sd_schmidt = dvisc / statevar.ρ / Dv(statevar.P,statevar.T)
sd_prandtl = dvisc / (HeatCoeff(statevar.T,statevar.qv)/cptot)

n_ventX = sd_schmidt^(1/3) * sqrt(sd_reynolds)

# kinetic correction
# ivt = 1 / statevar.T
# Fac_ice_kk_invrhoi = ( ivt * Lsub(0) / constants.Rv - 1 ) * (Lsub(0)/HeatCoeff(statevar.T,statevar.qv)) * ivt
# psychrometric correction
# Fac_ice_dd_invrhoi = constants.Rv / Dv(statevar.P,statevar.T) * statevar.T / statevar.esati
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

if statevar.ssi > 0
    if statevar.T < dep_params.temppoly && dep_params.mval < 2.0
        ihabit  = 3
        iequiv  = 1
    else
        iequiv = 0    # default spheroids
        ihabit = 1    # default columns T < -40
    end
    temp_celsius = statevar.T - constants.T0
    scrit_a, scrit_c = getscritdep(temp_celsius) 
    scrit_c, scrit_a = adjustscritdep(temp_celsius, ihabit, iequiv, scrit_c, scrit_a)
    # calculate alpha_a
    alpha_a = calculate_alpha_ac(a, c, caplengtha(a,c) ,statevar, constants, scrit_c, dep_params.mval)
    # calculate alpha_c
    alpha_c = calculate_alpha_ac(a, c,caplengthc(a,c), statevar, constants, scrit_a, dep_params.mval)
    alphat = 0.96 
else
    alpha_a = 1.0
    alpha_c = 1.0
    alphat = 1.0 #not sure about this one, 0.96 elsewhere
end


if alpha_a > 1e-5
    igr = max(min(alpha_c/alpha_a,100),0.01)

    if max_diameter<1.0e-5 
        igr=1.0 # set inherent growth ratio = 1 if D < 10um (Shima)
    end

    # !! if already too long, limit Gamma_star<=1 for deposition
    if phi>40
       igr = min(igr,1.0_RP)
    end

    if statevar.T < dep_params.temppoly && dep_params.mval < 2 #m = dep coeff change???
        igr = 1.0 # for polycrystal approx.
    end

    # recalculate anything that includes igr
    alphanr = a/radius^(3/(2+igr))  
end


# here calculate afn
afn = calc_afn(statevar,constants,fhdum,barf,a,c,alpha_a,alpha_c,alphat)

da, dc, dρ = dep_tendencies(afn,crystal,statevar,constants,dep_params,igr,fs,radius)
return da, dc, dρ
end



function dep_tendencies(afn,crystal,statevar,constants,dep_params,igr,fs,radius)

    if afn > 0 #deposition


        if statevar.T < dep_params.temppoly && dep_params.mval < 2.0  # Gwenore's polycrystals
            # qv_sd_max=qvs_cm1(k,i,j)
            statevar.qvs #not sure if this should be a statevar??
            satri_max  = statevar.qvs / statevar.qvsi #maybe turn the sat into a function 
            ssi_max  = satri_max - 1
            sirat = max(statevar.ssi/ssi_max,0)
            # if rho_dep is assigned during nucleation
            psat = crystal.ρdep
            if radius > dep_params.ro
            rho_dep = max(50.0,psat)
            else
            rho_dep = constants.ρ_ice
            end
            fs = 1 #change, is this supposed to be here
        else 
            sup =statevar.qv/statevar.qvs - 1
            if sup >= 0
                ssi_max = 1
            elseif ssi >= 0 && statevar.qvsi < statevar.qvs
                ssi_max = ((statevar.ssi+1)*statevar.qvsi)/(statevar.qvs-statevar.qvsi) - statevar.qvsi/(statevar.qvs-statevar.qvsi)
                ssi_max = min(ssi_max,1)
                ssi_max = max(ssi_max,0)
            else
                ssi_max = 0
            end

            if igr < 1  #Planar
            if crystal.vel > 0.0
                if crystal.a > sqrt((Dv(statevar.P,statevar.T)*pi*2*crystal.c)/(crystal.vel))
                    rho_dep = (constants.ρ_ice*igr)*ssi_max + constants.ρ_ice*(1-ssi_max)
                else
                    rho_dep = constants.ρ_ice
                end
            else
                    rho_dep = constants.ρ_ice
            end
            else # Columnar
            rho_dep = (constants.ρ_ice/igr)*ssi_max + constants.ρ_ice*(1-ssi_max)
            end

        end
        if radius > dep_params.ro
            rho_dep = min(rho_dep,constants.ρ_ice)
        end

    else # sublimation
            #.. During sublimation using polynomial removal of density
            rho_dep = crystal.ρ
            vmin  = 4/3 * pi *(dep_params.ro)^3
            if vmin < vol(crystal.a,crystal.c)
            betavol = log(constants.ρ_ice/crystal.ρ)*1/(log(vmin/vol(crystal.a,crystal.c)))
            rho_dep = crystal.ρ*(1+betavol)
            else
            rho_dep = crystal.ρ
            end

            #WTF change
            rho_dep = crystal.ρ    #!!KK : used in ICE-BALL  (comment it otherwise)    
            igr=1.0 # set inherent growth ratio = 1 for sublimation
    end

        #this happens after the dep/sub if statement but I think its bullshit   
        # rho_dep=max(rho_dep,50.0_RP)   !!KK: changed it from 100.
        # rho_dep=min(rho_dep,rhoi_mks)




    #..  characteristic r-axis and a-axis after growth timestep
    # d effective radius squared
    d_efr2 = 2.0*afn*fs/rho_dep # how to make a constraint on this? dep_params.min_vol
    da = 1/2*radius*d_efr2
    dc = 1/2*radius*d_efr2
    #   new_sd_ra = alphanr*rnf**(3.0_RP/(2.0_RP+igr))

    # #.. Do not sublimation change the shape of ice from
    # #.. prolate to oblate or vice versa
    #   sd_phif = phi*(rnf^3/radius^3)^((igr-1)/(igr+2))
    #   if statevar.ssi < 0.0 || afn < 0.0
    #     if phi < 1.0 && sd_phif < 1.0
    #         sd_phif=sd_phi
    #         alphanr = sd_ra/radius
    #         # new_sd_ra = alphanr*rnf
    #     end 
    #     if phi < 1.0 && sd_phif > 1.0 
    #         sd_phif = phi
    #         alphanr = sd_ra/radius
    #         # new_sd_ra = alphanr*rnf    
    #     end

    # #.. Do not let sublimation create extreme shapes
    #      if phi > 1
    #         if sd_phif > phi
    #            sd_phif = phi
    #            alphanr = sd_ra/radius
    #         #    new_sd_ra = alphanr*rnf                  
    #         end
    #      else
    #         if sd_phif < phi
    #            sd_phif = phi
    #            alphanr = sd_ra/radius
    #         #    new_sd_ra = alphanr*rnf                    
    #         end
    #      end
    #   end
    
    # #.. C-axis after vapor growth
    #    new_sd_rc = sd_phif*new_sd_ra   

    #change here if switching to not tendencies
    dvol = 2pi*radius*d_efr2
    dρ = (rho_dep-crystal.ρ)*dvol/radius^3
    #    new_sd_rho = sd_rho*(sd_vol/new_sd_vol) + rho_dep*(1.0_RP-sd_vol/new_sd_vol)
    #    new_sd_rho = min(new_sd_rho,rhoi_mks)
    #    new_sd_mass = new_sd_vol*new_sd_rho


    #   #################################################!
    #   ### limiter to avoid unphysical ice particle ###!
    #   #change, spherical for now..
    #   if crystal.a + da*t <= dep_params.ro  # ! Regard it is spherical with bulk ice density if smaller than 1um
    #      new_sd_rho = constants.ρ_ice
    #      new_sd_vol = new_sd_mass / new_sd_rho
    #      new_sd_ra = (new_sd_vol/F_THRD/ONE_PI)**(1.0d0/3.0d0)
    #      new_sd_rc = new_sd_ra
    #   end if

    #   if( min(new_sd_ra,new_sd_rc) <= 1.0d-9 )then ! Stop sublimation if smaller than 1nm
    #      new_sd_ra = 1.0d-9
    #      new_sd_rc = 1.0d-9
    #      new_sd_rho = rhoi_mks
    #   end if    

    #################################################!
    ### update the super-droplet ###!
    return da, dc, dρ
end