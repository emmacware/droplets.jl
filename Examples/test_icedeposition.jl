##################################################
# Ice Depositional Growth example
##################################################
using Random
using Distributions
using Plots
using Interpolations
using Revise # use this when developing Droplets, no reload required
using Droplets

# set random seed
Random.seed!(1234)


# This file has examples of ice depositional growth using 
# the vectorized structure and the crystal structure
# using deposition_timestep! from deposition.jl

# Settings 
#-------------------------------------------------
Ns = 2
Nx = 1
Ny = 1
Nz = 1
Δx = 100
Δy = 100
Δz = 100
n0 = 2^13
R0 = 30.531e-6 #meters
M0 = 1e-16 #kg

Δt = 1
ΔV = 10^6
deposition_type = "DepositionCoeff" # specify the type of deposition
# density_param = "param_here" # specify the density parameter
# constants = "constants_here" # specify the constants
# solver_type = "solver_here" # specify the solver type

radius_bins_edges = 10 .^ range(log10(10*1e-6), log10(5e3*1e-6), length=128) 
smooth = true

ξstart,Rstart,Xstart,Mstart = init_ξ_const(Ns,ΔV,n0,R0,M0)


# initialize the crystal density to be pure ice
ρistart =  ones(Ns)*917.0 #kg/m^3

# initialize the deposition density
ρdstart = max.(Random.rand(Ns)*917.0,65) #kg/m^3

# Initialize the crystals
crystals = create_crystals(Ns,Nx,Ny,Nz,Δx,Δy,Δz,
                           Rstart,Rstart,ξstart,Mstart,ρistart,ρdstart)

# initialize the state variables
state = create_statevar()

# run ice deposition
Δt = 1
rnew = ice_deposition(crystals[1],Δt,state)

print("initial radius: ", crystals[1].a, "\n")
print("new radius: ", rnew, "\n")