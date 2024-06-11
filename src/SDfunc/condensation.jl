#---------------------------------------------------------------------------
#Condensation functions in this file (more description throughout the file)
#---------------------------------------------------------------------------
using DifferentialEquations
export esat,sat,drdtcondensation1,drdtcondensation2,drdtcondensation3,condense_and_calc_Sv!
export FK,FD,drkohler,eq_radius,θcondenseupdate!,qvcondenseupdate!

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



using Roots
######################################################################
# Functions in the Denominator of the Köhler Equation

function FK(T) ##Temperature is the only one updating with time I think
    Fk = (L ./(Rv.*T) .-1).*(L*ρl)./(K.*T) 
    return Fk
end

function FD(T) ##Temperature is the only one updating with time I think
    Fd = ρl*Rv.*(T .+243.04)./(D.*esat(T))
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
    du = drkohler(u,M,m,T,qv,P)
    return du
end

function drdtcondensation3(u,p,t)
    M,m,T,Senv = p
    du = drkohler(R,M,m,T,Senv)
    return du
end

######################################################################
#This function calculate the LHS of the Köhler equation
# drdt = ___
# Can take either the (mixing ratio and pressure) or (environmental saturation) as an argument

function drkohler(R,M,m,T,qv,P)
    a=3.3*10^(-7)/T #(m K/T)
    b=4.3 *2 ./m ./1e6 #m^3 for NaCL
    S = sat.(qv,P) ./esat(T)
    denom = (FK(T)+FD(T))
    return (S-1-(a./R) .+b .*M ./(R.^3)) ./(denom.*R)
end

function drkohler(R,M,m,T,Senv)
    a=3.3*10^(-7)/T #(m K/T)
    b=4.3 *2 ./m ./1e6 #m^3 for NaCL
    S = Senv/esat(T)
    denom = (FK(T)+FD(T))
    return (S-1-(a./R) .+b .*M ./(R.^3)) ./(denom.*R)
end

######################################################################
# eq_radius calculates the equilibrium radius of a droplet
# This equation caps the environmental saturation to 1 to be able to find the real root,
# making sure that Senv<Sactivation
# Can take either the (mixing ratio and pressure) or (environmental saturation) as an argument
# Tolerances are arbitrary but tested on multiple examples

function eq_radius(m,M,qv,P,T)
    Senv = sat(qv,P)/esat(T)
    a = 3.3*10^(-7)/T #(m)
    b = 4.3 *2 ./m./1e6 #m^3
    # Cap Senv to find real root
    if Senv > 1
        Senv = 1
    end
    FindR(r) = Senv - (1+a/r - (b*M) /r^3)      
    eqr = find_zero(FindR, (2e-8,4e-7), xtol=1e-12) #might need different tolerances?
    return eqr
end

function eq_radius(m,M,Senv,T)
    S = Senv/esat(T)
    a = 3.3*10^(-7)/T #(m)
    b = 4.3 *2 ./m./1e6 #m^3
    # Cap Senv to find real root
    if Senv > 1
        Senv = 1
    end
    FindR(r) = S - (1+a/r - (b*M) /r^3)      
    eqr = find_zero(FindR, (2e-8,4e-7), xtol=1e-12) #might need different tolerances?
    return eqr
end

######################################################################
# condense_and_calc_Sv! calculates the change in volume of droplets in a grid cell
# There are a few options here for arguments, but they all use the same ODE function
# option 1: qvarray,T,P,ρ,Δtg,ΔgridV,Nx,Ny,grid_dict
#       -this option is set up how I was starting to implement it in 2D
#       it takes the droplets as structs and calculates the condensation and
#       Sv for each grid cell
#       -it is not finished and not usefull without the specific organization scheme
#       grid_dicts with lists of mutable structs
# option 2: R,ξ,X,M,m,qvarray,T,P,ρ,Δtg,ΔV
#       -this option is set up to take in the droplet properties as vectors
#       but is meant for a box or single grid cell: P, T, qv, and ρ should be scalars
#       -the difference between this and option 3 is that this takes in the mixing ratio 
#       and pressure as arguments instead of the environmental saturation
# option 3: R,ξ,X,M,Senv,T,ρ,Δtg,ΔV
#       -this option is set up to take in the droplet properties as vectors
#       meant for a box or single grid cell: T, Senv, and ρ should be scalars
#       -this takes in the environmental saturation instead of the mixing ratio and pressure

#Not sure how to implement this in 2D looping over the full Ns superdroplets
#without splitting it up by grid yet.. is that possible?

function condense_and_calc_Sv!(qvarray,T,P,ρ,Δtg,ΔgridV,Nx,Ny,grid_dict)
    Sv = zeros(Nx,Ny)
    for i in 1:Nx
        for j in 1:Ny
            grid = (i,j)
            top = 0
            Ngrid = length(grid_dict[(i,j)])
            a=3.3*10^(-7)./T[i,j] #(m K/T)
            S = sat(qvarray[i,j],P[i,j])/esat(T[i,j])
            denom = (FK(T[i,j])+FD(T[i,j]))
            if isempty(grid_dict[i,j])
                continue
            else

                dropletspre = ([droplet.X*droplet.ξ*ρl for droplet in grid_dict[i,j]]).+0.0

                # condensationpergridcellradius!(grid_dict[i,j],qvarray[i,j],T[i,j],P[i,j],Δtg)
                for droplet in grid_dict[i,j]
                    b=4.3 *2/droplet.m/1e6 #m^3 for NaCL
                    pcondense = (a,b,S,droplet.M,denom)
                    tspan = (0.0,Δtg)
                    Rc = droplet.R+0.0
                    # rODE = ODEProblem{false}(drdtcondensation1,Rc,tspan,pcondense)
                    rODE = ODEProblem(drdtcondensation1,Rc,tspan,pcondense)

                    # droplet.R = solve(rODE,AutoTsit5(Rosenbrock23()), dt=Δtg).u[end]
                    droplet.R = solve(rODE,ImplicitEuler()).u[end]
                    if droplet.R <=0 #failsafe, should change to dry radius
                        droplet.R = 1e-7
                    end
                    droplet.X = 4/3 * π * droplet.R^3
                end
                

                dropletspost = ([droplet.X*droplet.ξ*ρl for droplet in grid_dict[i,j]]).+0.0
                top = sum(dropletspost.-dropletspre)/Δtg
                # print(top)
            end
            Sv[i,j] = -top/(ρ[i,j] * ΔgridV)
            
        end
    end
    return Sv
end


function condense_and_calc_Sv!(R,ξ,X,M,m,qvarray,T,P,ρ,Δtg,ΔV)

        a=3.3*10^(-7)./T[i,j] #(m K/T)
        b=4.3 *2 ./m./1e6 #m^3 for NaCL
        S = sat(qvarray,P)/esat(T)
        denom = (FK(T)+FD(T))
        dropletspre = (X.*ξ.*ρl).+0.0
        pcondense = (a,b,S,M,denom)
        tspan = (0.0,Δtg)
        Rc = R.+0.0

        #still using drdtcondensation1 because calculating S, M, and denom outside is more efficient
        #becuase we are doing this for every droplet in the vectors

        rODE = ODEProblem{false}(drdtcondensation1,Rc,tspan,pcondense) 
        R = solve(rODE,AutoTsit5(Rosenbrock23()), dt=Δtg).u[end]
        # R = solve(rODE,ImplicitEuler()).u[end]
        #update X
        X = 4/3 * π .* R.^3
        dropletspost = (X.*ξ.*ρl).+0.0

        Sv = sum(dropletspost.-dropletspre)/(Δtg* ρ * ΔV)
    return R,X,Sv
end

function condense_and_calc_Sv!(R,ξ,X,M,m,Senv,T,ρ,Δtg,ΔV)

    a=3.3*10^(-7)./T #(m K/T)
    b=4.3 *2 ./m./1e6 #m^3 for NaCL
    S = Senv/esat(T)
    denom = (FK(T)+FD(T))
    dropletspre = (X.*ξ.*ρl).+0.0
    pcondense = (a,b,S,M,denom)
    tspan = (0.0,Δtg)
    Rc = R.+0.0 #necessary?
    #still using drdtcondensation1 because calculating S, M, and denom outside is more efficient
    #becuase we are doing this for every droplet in the vectors
    rODE = ODEProblem{false}(drdtcondensation1,Rc,tspan,pcondense)
    R = solve(rODE,AutoTsit5(Rosenbrock23()), dt=Δtg).u[end]
    # R = solve(rODE,ImplicitEuler()).u[end]
    X = 4/3 * π .* R.^3

    dropletspost = (X.*ξ.*ρl).+0.0
    Sv = sum(dropletspost.-dropletspre)/(Δtg* ρ * ΔV)
    return R,X,Sv
end

######################################################################
# θcondenseupdate! and qvcondenseupdate! are functions to update the 
# temperature and mixing ratio after condensation, using Sv
# both can take the arguments as vectors or scalars (Δtg needs to be scalar)

function θcondenseupdate!(Sv,θ,Δtg,P,P0,constants)
    Exner = (P./P0).^(Constants.Rd/constants.Cp)
    θ = θ .+ -Δtg*constants.L.*Sv./(constants.Cp.*Exner)
    T = θ.*(P./P0).^(constants.Rd/constants.Cp)
    return θ,T
end

function qvcondenseupdate!(Sv, qvarray, P,T,constants,Δtg)
    qvarray = qvarray .+ Δtg.*Sv
    ρd =  P./(constants.Rd.*T)
    ρ = ρd ./(1 .- qvarray) 
    return qvarray,ρ
end

