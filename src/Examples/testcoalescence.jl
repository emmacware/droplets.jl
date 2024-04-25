##################################################
# Coalescence example
##################################################
using Random
using Combinatorics
using Distributions
using Plots
using Interpolations
include("../SDfunc/coalescence.jl")
include("../SDfunc/DSDvis.jl")
include("../SDfunc/setup.jl")

# This file has examples of box collision coalescence using 
# the vectorized structure and the droplet structure
# using coalescence_timestep! from coalescence.jl


# ##################################################
# # Coalescence example, with vectorized structure
# ##################################################
# #Settings 
# #-------------------------------------------------
Ns = 2^15
n0 = 2^23
R0 = 30.531e-6 #meters
Δt = 1
ΔV = 10^6
kernel = golovin
radius_bins_edges = 10 .^ range(log10(10*1e-6), log10(5e3*1e-6), length=128) 
smooth = true

ξstart,Rstart,Xstart = init_ξ_const(Ns,ΔV,n0,R0)




function run_small_alpha!(ξ,R,X,Ns,Δt,ΔV;smooth=true,kernel=golovin,radius_bins_edges=radius_bins_edges)
    println("Running simulation...")
    seconds = 0

    plot1 = plot()
    coal_func_time = 0
    bins = zeros(length(radius_bins_edges)-1)

    simtime = @elapsed begin
    while seconds <= 3600
    
        if seconds >= 0 && mod(seconds,1200) == 0
            println("Time: ", seconds, " seconds")

            xx,yy = binning_dsd(X,ξ,radius_bins_edges,seconds,smooth=smooth,scope_init=2)
            bins = hcat(bins,yy)

            plot1 = plot!(xx,yy,lc=:red,
                # linetype=:steppost,
                # linewidth=2,
                xaxis=:log,xlims=(10,10000), 
                label="t = "*string(seconds)*" sec",legend=:outerright,
                legendtitle = string(kernel)*" kernel, "*s"$N_s$=2^"*string(Int(log2(Ns))),guidefontsize=10)
            title!("Time Evolution of the Mass Density Distribution")
            xlabel!("Droplet Radius R (μm)")
            ylabel!("Mass Density Distribution g(lnR)(g/m^3/unit ln R)")
        
        end

        time = @elapsed begin
        ξ,R,X = coalescence_timestep!(ξ,R,X,Ns,Δt,ΔV,kernel=kernel)
        end # coalescence_timestep function timing
        coal_func_time += time
        X = 4/3*π * R.^3

        seconds = seconds + Δt
    end
    end # simulation time
    println("simtime =",simtime)
    println("coal_func_time =",coal_func_time)

    return plot1
end

plot1 = run_small_alpha!(ξstart,Rstart,Xstart,Ns,Δt,ΔV,smooth=smooth,kernel=golovin,radius_bins_edges=radius_bins_edges)











# ##################################################
# # Coalescence example, with droplet structure
# ##################################################

# #Settings
# #-------------------------------------------------
# Ns = 2^15
# n0 = 2^23
# R0 = 30.531e-6 #meters
# M0 = 1e-16 #kg
# Δt = 1
# ΔV = 10^6
# Nx = Ny = 100
# Δx = Δy = 10
# kernel = golovin
# radius_bins_edges = 10 .^ range(log10(10*1e-6), log10(5e3*1e-6), length=128) 
# smooth = true

# ξstart,Rstart,Xstart,Mstart = init_ξ_const(Ns,ΔV,n0,R0,M0)
# superdrops = create_NaCl_superdroplets(Ns,Nx,Ny,Δx,Δy,Rstart,Xstart,ξstart,Mstart)

function run_small_alpha!(droplets,Ns,Δt,ΔV;smooth=true,
    kernel=golovin,radius_bins_edges=10 .^ range(log10(10*1e-6), log10(5e3*1e-6), length=128))
    
    println("Running simulation...")
    seconds = 0

    plot1 = plot()

    coal_func_time = 0
    bins = zeros(length(radius_bins_edges)-1)

    simtime = @elapsed begin
    while seconds <= 3600
   
        
        if seconds >= 0 && mod(seconds,1200) == 0
            println("Time: ", seconds, " seconds")
            # droplets = sort(droplets, by = droplet -> droplet.R)

            xx,yy = binning_dsd(droplets,radius_bins_edges,seconds,smooth=smooth,scope_init=2)
            bins = hcat(bins,yy)

                plot1 = plot!(xx,yy,lc=:red,
                    # linetype=:steppost,
                    # linewidth=2,
                    xaxis=:log,xlims=(10,10000), 
                    label="t = "*string(seconds)*" sec",legend=:outerright,
                    legendtitle = string(kernel)*" kernel, "*s"$N_s$=2^"*string(Int(log2(Ns))),guidefontsize=10)

                title!("Time Evolution of the Mass Density Distribution")
                xlabel!("Droplet Radius R (μm)")
                ylabel!("Mass Density Distribution g(lnR)(g/m^3/unit ln R)")

                plot!(plot1)
        end

        time = @elapsed begin
        droplets = coalescence_timestep!(droplets,Ns,Δt,ΔV,kernel=kernel)
        # println("Time: ", seconds)
        end # coalescence_timestep function timing
        coal_func_time += time


        seconds = seconds + Δt
    end
    end # simulation time
    println("simtime =",simtime)
    println("coal_func_time =",coal_func_time)

    bins = bins[:,2:end]

    return plot1
end


# plot1 = run_small_alpha!(superdrops,Ns,Δt,ΔV,smooth=smooth,kernel=golovin,radius_bins_edges=radius_bins_edges)