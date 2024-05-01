using Random
using Combinatorics
using Distributions
using Plots
using Interpolations
using DifferentialEquations
include("../SDfunc/coalescence.jl")
include("../SDfunc/DSDvis.jl")
include("../SDfunc/setup.jl")
include("../SDfunc/constants.jl")
include("../SDfunc/condensation.jl")
include("../SDfunc/updateposition.jl")


# ##################################################
# Sedimentation example, combining condensation and coalescence,
# allowing for vertical movement of droplets

# ##################################################

# #Settings 

#spatial domain (single column)
#-------------------------------------------------
Nx = 1
Ny = 5
Δx = 100
Δy = 100
qv = [0.01, 0.01, 0.01, 0.01, 0.01] #kg/kg
P = [101.3, 101.3, 101.3, 101.3, 101.3] #hPa

#-------------------------------------------------

# Create the grid with lower and upper bounds of dx and dy
grid_box = Array{Tuple{Float64, Float64,Float64,Float64}}(undef, Nx, Ny)
grid_box_mids_x = zeros(Nx, Ny)
grid_box_mids_y = zeros(Nx, Ny)
for i in 1:Nx
    for j in 1:Ny
        dx_lower = (i-1) * Δx
        dx_upper = i * Δx
        dy_lower = (j-1) * Δy
        dy_upper = j * Δy
        mid_x = (dx_lower + dx_upper)/2
        mid_y = (dy_lower + dy_upper)/2
        grid_box[i, j] = (dx_lower, dx_upper, dy_lower, dy_upper)
        grid_box_mids_x[i, j] = mid_x
        grid_box_mids_y[i, j] = mid_y
    end
end

# #-------------------------------------------------
Sv = zeros(Nx,Ny)
P0 = 101325 #Pa
T = zeros(Nx,Ny).+305
for j in 1:Ny
    height = grid_box_mids_y[1,j]/1000
    T[:,j] .= 305-9.5*height
end

qvarray = [0.0295 0.028 0.027 0.026 0.025] #kg/kg, specific humidity or mixing ratio of vapor/moist air
P = P0.*exp.(-gconst.*grid_box_mids_y./(T.*Rd))
ρu = zeros(Nx,Ny) #just sedimentation
ρv = zeros(Nx,Ny+1) #just sedimentation
ρd =  P./(Rd.*T)
ρ = ρd ./(1 .- qvarray) #kg/m^3, density of moist air

#superdroplets
# #-------------------------------------------------
Ns = 2^12
n0 = 2^15
R0 = 30.531e-6 #meters
M0 = 1e-16 #kg
Δt = 1
ΔV = 10^2*Nx*Δx*Ny*Δy #depth of 100m, full system
kernel = golovin
radius_bins_edges = 10 .^ range(log10(10*1e-6), log10(5e3*1e-6), length=128) 
smooth = true

ξstart,Rstart,Xstart,Mstart = init_ξ_const(Ns,ΔV,n0,R0,M0)

superdroplets = create_NaCl_superdroplets(Ns,Nx,Ny,Δx,Δy,Rstart,Xstart,ξstart,Mstart)
grid_dict = Dict{Tuple{Int, Int}, Vector{Superdroplet}}()
droplet_gridbox(superdroplets,Nx,Ny,Δx,Δy,grid_dict)

# #-------------------------------------------------
anim = Animation()
plot_grid_with_droplets(grid_box, superdroplets)
ylims!(0,500)
frame(anim)
sleep(1)

for t in 1:15
    for _ in 1:60
        for i in 1:Ny
            droplets = grid_dict[(1,i)]
            ngrid = length(droplets)
            coalescence_timestep!(droplets,ngrid,Δt,ΔV,kernel=kernel)
        end
        condense_and_calc_Sv!(qvarray,T,P,ρ,Δt,ΔV,Nx,Ny,grid_dict)
        qvcondenseupdate!(Sv, qvarray, P,T, Δt)
        update_position!(superdroplets,Nx,Ny,Δt,ρu,ρv,ρ, grid_dict,grid_box)
    end
    plot_grid_with_droplets(grid_box, superdroplets)
    ylims!(0,500)
    frame(anim)
    sleep(1)
end
gif(anim, "sediment.gif", fps = 1)

# plot_grid_with_droplets(grid_box, superdroplets)
# plot(collect(drop.R for drop in superdroplets))