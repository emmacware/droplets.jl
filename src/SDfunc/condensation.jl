#---------------------------------------------------------------------------
#Condensation functions in this file (more description throughout the file)
#---------------------------------------------------------------------------
# using DifferentialEquations
export esat,sat,drdtcondensation1,drdtcondensation2,drdtcondensation3#,condense_and_calc_Sv!
export FK,FD,drkohler,θcondenseupdate!,qvcondenseupdate!,dXkohler_function_of_radius
export dXkohler_function_of_radius_activated,drkohler_activated

# FK(T),            returns FK in the Köhler equation
# FD(T),            returns FD in the Köhler equation
# esat(T),          returns saturation vapor pressure
# sat(qvarray,P),   returns environmental saturation
#
# drdtcondensation1(u,p,t)      ODE function for condensation, p = (a,b,S,M,denom)
# drdtcondensation2(u,p,t)      ODE function for condensation, p = (M,m,T,qv,P)
# drdtcondensation3(u,p,t)      ODE function for condensation, p = (M,m,T,Senv) 
#
# drkohler(R,M,m,T,qv,P), returns drdt (LHS of Köhler equation)
# drkohler(R,M,m,Senv),   returns drdt (LHS of Köhler equation)
#
# eq_radius(m,M,qv,P,T),  returns equilibrium radius for droplet activation (considering subsaturation)
# eq_radius(m,M,Senv,T),  returns equilibrium radius for droplet activation (considering subsaturation)
#
# condense_and_calc_Sv!(qvarray,T,P,ρ,Δtg,ΔgridV,Nx,Ny,grid_dict),  returns Sv,grid_dict
# condense_and_calc_Sv!(R,ξ,X,M,m,qvarray,T,P,ρ,Δtg,ΔV),            returns R,X,Sv
# condense_and_calc_Sv!(R,ξ,X,M,Senv,T,ρ,Δtg,ΔV),                   returns R,X,Sv
#
# θcondenseupdate!(Sv, θ, Δtg,P),           returns θ,T
# qvcondenseupdate!(Sv, qvarray, P,T, Δtg), returns ρ*qvarray,ρ



# using Roots #put this in examples
######################################################################
# Functions in the Denominator of the Köhler Equation

function FK(T) ##Temperature is the only one updating with time I think
    Fk = (constants.L ./(constants.Rv.*T) .-1).*(constants.L*constants.ρl)./(constants.k .*T) 
    return Fk
end

function FD(T) ##Temperature is the only one updating with time I think
    # Fd = constants.ρl*constants.Rv .*(T .+243.04)./(constants.Dv .* esat(T))
    Fd = constants.ρl*constants.Rv .*(T)./(constants.Dv .* esat(T))
    return Fd
end

######################################################################
# Saturation Functions

#saturation vapor pressure
esat(T) = 100*6.1094.*exp.(17.625.*(T.-273)./((T.-273) .+243.04)) ## August–Roche–Magnus approximation, hPa to Pa, T in K, converted to C

#environmntal saturation based on the mixing ratio of water vapor to air
sat(qvarray,P) = qvarray.*P./(0.378.*qvarray .+ 0.622) ##qvarray is the specific humidity: mixing ratio of water vapor to moist air, Pa


######################################################################
# The following ode functions are for the condensation radial growth of a droplet
# the argument u is the radius of the droplet
# there are three options given the available parameters
#drdtcondensation1 p = (a,b,S,M,denom)
#drdtcondensation2 p = (M,m,T,qv,P)
#drdtcondensation3 p = (M,m,T,Senv)
#where S is environmental saturation, M is the mass of the solute,
#m is the molecular weight of the solute, T is the temperature, 
#qv is the mixing ration, and P is the pressure


function drdtcondensation1(u, p, t)
    a, b, S, M,denom = p
    return (S-1-(a ./u) .+b.*M./(u.^3))./(denom.*u)
end

function drdtcondensation2(u,p,t)
    M,m,T,qv,P = p
    du = drkohler(u,M,m,T,qv,P,t)
    return du
end

function drdtcondensation3(u,p,t)
    M,m,T,Senv = p
    du = drkohler(u,M,m,T,Senv,t)
    return du
end

######################################################################
#This function calculate the RHS of the Köhler equation
# drdt = ___
# Can take either the (mixing ratio and pressure) or (environmental saturation) as an argument

function drkohler(R,M,m,T,qv,P,timestep)
    a=3.3*10^(-7)/T #(m K/T)
    b=4.3 *2 ./m ./1e6 #m^3 for NaCL
    S = sat.(qv,P) ./esat(T)
    denom = (FK(T)+FD(T))
    dr = (S-1 .-(a./R) .+b .*M ./(R.^3)) ./(denom.*R)
    return R + dr*timestep > 0 ? dr : -R/timestep
end

function drkohler(R,M,m,T,Senv,timestep)
    a=3.3*10^(-7)/T #(m K/T)
    b=4.3 *2 ./m ./1e6 #m^3 for NaCL
    # S = Senv/esat(T)
    denom = (FK(T)+FD(T))
    dr = (Senv-1 .-(a./R) .+b .*M ./(R.^3)) ./(denom.*R)
    return R + dr*timestep > 0 ? dr : -R/timestep
end

function drkohler_activated(R,T,Senv,timestep)
    # S = Senv/esat(T)
    denom = (FK(T)+FD(T))
    dr = (Senv-1) ./(denom.*R)
    return R + dr*timestep > 0 ? dr : -R/timestep
end

function dXkohler_function_of_radius(R,M,m,T,qv,P,timestep)
    dX = 4*π*R^2*drkohler(R,M,m,T,qv,P)
    return dX
end

function dXkohler_function_of_radius(R,M,m,T,Senv,timestep)
    dX = 4*π*R^2*drkohler(R,M,m,T,Senv,timestep)
    return dX
end

function dXkohler_function_of_radius_activated(R,T,Senv,timestep)
    dX = 4*π*R^2*drkohler_activated(R,T,Senv,timestep)
    return dX
end


function dq_liq_cond_activated(R,M,m,T,Senv,timestep,ρ_air)
    dX = dXkohler_function_of_radius(R,M,m,T,Senv,timestep)
    dql = sum(dX.*ξ.* constants.ρl / ρ_air) #/per volume?
    return dql
end

function dq_liq_cond_activated(R,T,Senv,timestep,ρ_air)
    dX = dXkohler_function_of_radius_activated(R,T,Senv,timestep)
    dql = sum(dX.*ξ.* constants.ρl / ρ_air) #/per volume?
    return dql
end



# function θcondenseupdate!(Sv,θ,Δtg,P,P0)
#     Exner = (P./P0).^(constants.Rd/constants.Cp)
#     θ = θ .+ -Δtg*constants.L.*Sv./(constants.Cp.*Exner)
#     T = θ.*(P./P0).^(constants.Rd/constants.Cp)
#     return θ,T
# end

# function qvcondenseupdate!(Sv, qvarray, P,T,Δtg)
#     qvarray = qvarray .+ Δtg.*Sv
#     ρd =  P./(constants.Rd.*T)
#     ρ = ρd ./(1 .- qvarray) 
#     return qvarray,ρ
# end

