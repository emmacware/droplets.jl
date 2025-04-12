#---------------------------------------------------------
#Collision Coalescence Functions
#---------------------------------------------------------
#currently, this collision-coalescence is only set up to handle radius, 
#multiplicity, volume, and solute mass as superdroplet attributes
using Distributions
using Combinatorics
using Random
using Interpolations

#types
export Serial, Parallel,Adaptive,none,coag_settings,droplet_attributes
#initialization methods
export init_ξ_const,init_logarithmic,init_uniform_sd, init_monodisperse
#kernels
export terminal_v,hydrodynamic,golovin,coalescence_timestep!
#coagulation timestep
export coagulation_run
# SDM logic
export adaptive_pαdt!, pair_Ps_adaptive!,compute_pαdt!,split_highest_multiplicity!,test_pairs!,pair_Ps!,sdm_update!


# fix and export: collision_efficiency
struct Serial end
struct Parallel end
struct Adaptive end
struct none end

"""
    coag_settings{FT<:AbstractFloat}

A struct representing the settings for coalescence.

# Fields
- `FT`: AbstractFloat.
- `Δt`: The time step for the coalescence simulation.
- `ΔV`: The volume of the coalescence domain.
- `Ns`: The number of superdroplets.
- `scale`: The scaling factor that allows for linear sampling.
- `R_min`: The minimum radius of sampled droplets, used in some initializations. Radius in meters
- `R_max`: The maximum radius of sampled droplets, used in some initializations. Radius in meters.
- `golovin_kernel_coeff`: The Golovin kernel coefficient, 1/seconds.
- `hydrodynamic_collision_eff_func`: A boolean indicating whether to use the
        hydrodynamic collision efficiency function. Currently not implemented.
- `kernel`: The kernel function used for coalescence. Implemented options are `golovin` or `hydrodynamic`,
    can take any user implemented function with input type (::droplet_attributes, pairindex::Tuple{Int,Int}, ::coag_settings{FT}).
- `n0`: The initial real world droplet concentration.
- `R0`: The initial seed radius of the droplets.
"""
Base.@kwdef struct coag_settings{FT<:AbstractFloat}
    Δt::FT = FT(1.0)
    ΔV::FT = FT(1e6)
    Ns::Int = 2^15# number of superdroplets
    scale::FT = Ns * (Ns - 1) / 2 / (Ns / 2)
    R_min::FT = FT(1e-9)
    R_max::FT = FT(1e-3)
    golovin_kernel_coeff::FT = FT(1.5e3)
    hydrodynamic_collision_eff_func::Bool = false
    kernel::Function = golovin # golovin, hydrodynamic
    n0::FT = FT(2^23) # initial droplet concentration
    R0::FT = FT(30.531e-6) # initial radius
end

"""
    struct droplet_attributes{FT<:AbstractFloat}

A struct representing the attributes of a droplet.

Fields:
- FT: AbstractFloat.
- ξ: Vector of integers representing the multiplicity of each droplet.
- R: Vector of floats representing the radius of each droplet.
- X: Vector of floats representing the volume of each droplet.

"""
struct droplet_attributes{FT<:AbstractFloat}
    ξ::Vector{Int}
    R::Vector{FT}
    X::Vector{FT}
end


"""
    struct coagulation_run{FT<:AbstractFloat}

Struct initialising temp memory used for coagulation.

Input:
- Ns::Int : Number of superdroplets.

Fields:
- I::Vector{Int} : Vector of indexes to be shuffled for permutations.
- pαdt::Vector{FT} : Vector of coalescence probabilities for each pair.
- ϕ::Vector{FT} : Random numbers to be used in Monte Carlo coalescence.
- lowest_zero::Ref{Bool} : Reference to a boolean indicating if the lowest multiplicity is zero.
- deficit::Ref{FT} : Reference to a float representing the deficit in mass or volume.

"""
struct coagulation_run{FT<:AbstractFloat}
    Ns::Int
    I::Vector{Int}
    pαdt::Vector{FT}
    ϕ::Vector{FT}
    lowest_zero::Ref{Bool}
    deficit::Ref{FT}

    function coagulation_run{FT}(Ns::Int) where FT<:AbstractFloat
        I = collect(1:Ns)
        pαdt = zeros(FT, div(Ns, 2))
        ϕ = zeros(FT, div(Ns, 2))
        lowest_zero = Ref(false)
        deficit = Ref(zero(FT))
        new{FT}(Ns, I, pαdt, ϕ, lowest_zero, deficit)
    end
end


#---------------------------------------------------------
# INITIALIZATION
#---------------------------------------------------------

"""
    init_ξ_const(settings::coag_settings{FT}) where FT<:AbstractFloat

init_ξ_const initializes droplets based on the constant-multiplicity initialization method,
using an exponential distribution around the initial volume of the droplets. The method is
as described in Shima et al. (2009) https://doi.org/10.5194/acp-9-4491-2009

Arguments
- `settings`: Coagulation settings, type ::coag_settings.

"""
function init_ξ_const(settings::coag_settings{FT}) where FT<:AbstractFloat
    Ns = settings.Ns
    ΔV = settings.ΔV
    n0 = settings.n0
    R0 = settings.R0
    ξstart::Vector{Int} = (div(n0*ΔV,Ns)*ones(Ns))
    X0 = radius_to_volume(R0)# initial volume m3    
    Xstart::Vector{FT} = (rand(Exponential(X0), Ns))
    Rstart::Vector{FT} = volume_to_radius.(Xstart)

    droplets = droplet_attributes(ξstart, Rstart, Xstart)
    return droplets
end

"""
    init_logarithmic(settings::coag_settings{FT}) where FT<:AbstractFloat

init_logarithmic initializes the superdroplets, using a logarithmically spaced droplet radius spectrum,
initializing the multiplicities so that the water volume in the system forms an exponential distribution
around the initial volume.

Arguments
- `settings`: Coagulation settings.

"""

function init_logarithmic(settings::coag_settings{FT})where FT<:AbstractFloat
    Ns = settings.Ns
    ΔV = settings.ΔV
    n0 = settings.n0
    R0 = settings.R0
    R_min = settings.R_min
    R_max = settings.R_max

    X0 = radius_to_volume(R0)
    exp_dist = Exponential(X0)
    boundaries_found = false
    while boundaries_found == false
        radius_bins = 10 .^ range(log10(R_min), log10(R_max), length=Ns+1)
        # Calculate the volume for each bin
        vmin = radius_to_volume(radius_bins[1])
        vmax = radius_to_volume(radius_bins[end-1])
        
        pdf_min = pdf(exp_dist, vmin)
        pdf_max = pdf(exp_dist, vmax)

        first_bin_dr = radius_bins[2]-radius_bins[1]
        last_bin_dr = radius_bins[end]-radius_bins[end-1]
        dvdr_first = 4 * π * radius_bins[1].^2
        dvdr_last = 4 * π * radius_bins[end-1].^2

        ξ_first = pdf_min * first_bin_dr * dvdr_first * (n0*ΔV)
        ξ_last = pdf_max * last_bin_dr * dvdr_last * (n0*ΔV)

        if ξ_first >= 1 && ξ_last >=1
            boundaries_found = true
        else 
            if ξ_first <= 1
                R_min *= 1.01
            end
            if ξ_last <= 1
                R_max *= 0.99
            end
        end
    end

    radius_bins_new = 10 .^ range(log10(R_min), log10(R_max), length=Ns+1)
    sd_radii::Vector{FT} = [rand(Uniform(radius_bins_new[i], radius_bins_new[i+1])) for i in 1:Ns]
    volumes::Vector{FT} = radius_to_volume.(sd_radii)
    pdf_values = pdf.(exp_dist, volumes)
    bin_widths_new = diff(radius_bins_new)
    dvdr = 4 * π .* sd_radii.^2
    multiplicities = pdf_values .* bin_widths_new .* dvdr * (n0*ΔV)
    ξstart::Vector{Int} = floor.(multiplicities.+0.5)

    droplets = droplet_attributes(ξstart, sd_radii, volumes)
    return droplets
end

"""
init_uniform_sd(settings::coag_settings{FT}) where FT<:AbstractFloat

init_uniform_sd initializes the superdroplets, using an evenly spaced droplet radius spectrum,
initializing the multiplicities so that the water volume in the system forms an exponential distribution
around the initial volume.

Arguments
- `settings`: Coagulation settings.

"""

function init_uniform_sd(settings::coag_settings{FT})where FT<:AbstractFloat
    Ns = settings.Ns
    ΔV = settings.ΔV
    n0 = settings.n0
    R0 = settings.R0
    R_min = settings.R_min
    R_max = settings.R_max

    X0 = radius_to_volume(R0)
    exp_dist = Exponential(X0)
    boundaries_found = false
    while boundaries_found == false
        radius_bins = range(R_min, R_max, length=Ns+1)
        # Calculate the volume for each bin
        vmin = radius_to_volume.(radius_bins[1])
        vmax = radius_to_volume.(radius_bins[end-1])
        
        pdf_min = pdf(exp_dist, vmin)
        pdf_max = pdf(exp_dist, vmax)

        first_bin_dr = radius_bins[2]-radius_bins[1]
        last_bin_dr = radius_bins[end]-radius_bins[end-1]
        dvdr_first = 4 * π * radius_bins[1].^2
        dvdr_last = 4 * π * radius_bins[end-1].^2

        ξ_first = pdf_min * first_bin_dr * dvdr_first * (n0*ΔV)
        ξ_last = pdf_max * last_bin_dr * dvdr_last * (n0*ΔV)

        if ξ_first >= 1 && ξ_last >=1
            boundaries_found = true
        else 
            if ξ_first <= 1
                R_min *= 1.01
            end
            if ξ_last <= 1
                R_max *= 0.99
            end
        end
    end

    radius_bins_new = range(R_min, R_max, length=Ns+1)
    sd_radii::Vector{FT} = [rand(Uniform(radius_bins_new[i], radius_bins_new[i+1])) for i in 1:Ns]
    volumes::Vector{FT} = radius_to_volume.(sd_radii)
    pdf_values = pdf.(exp_dist, volumes)
    bin_widths_new = diff(radius_bins_new)
    dvdr = 4 * π .* sd_radii.^2
    multiplicities = pdf_values .* bin_widths_new .* dvdr * (n0*ΔV)
    ξstart::Vector{Int} = floor.(multiplicities.+0.5)

    droplets = droplet_attributes(ξstart, sd_radii, volumes)
    return droplets
end

"""
    init_monodisperse(settings::coag_settings{FT}) where FT<:AbstractFloat

init_monodisperse initializes the superdroplets so that all droplets have the same attributes.
Arguments
- `settings`: Coagulation settings.

"""

function init_monodisperse(settings::coag_settings{FT})where FT<:AbstractFloat
    Ns = settings.Ns
    ΔV = settings.ΔV
    n0 = settings.n0
    R0 = settings.R0
    ξstart::Vector{Int} = (div(n0*ΔV,Ns)*ones(Ns))
    Rstart::Vector{FT} = (R0*ones(Ns))
    Xstart::Vector{FT} = radius_to_volume.(Rstart)

    droplets = droplet_attributes(ξstart, Rstart, Xstart)
    return droplets
end



#---------------------------------------------------------
# HYDRODYNAMIC KERNEL
#---------------------------------------------------------

# terminal velocity of droplets
"""
    terminal_v(r::FT) where FT<:AbstractFloat

Compute the terminal velocity of a droplet of radius r, using tables from
Gunn and Kinzer, (1949), https://doi.org/10.1175/1520-0469(1949)006<0243:TTVOFF>2.0.CO;2

# Arguments
- `r::FT`: The radius of the droplet in meters.

# Returns
The terminal velocity of the droplet in meters/second.

"""
function terminal_v(r::FT)::FT where FT<:AbstractFloat  # terminal velocity 


    if 2*r*100<0.0078
        tv=1.2*10e6*(r*100)^2
        tv = tv/100
    else    
        d_table = [0.078,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,
            2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0,5.2,5.4,5.6,5.8]./10

        v_table = [18,27,72,117,162,206,247,287,327,367,403,464,517,565,609,649,690,
        727,757,782,806,826,844,860,872,883,892,898,903,907,909,912,914,916,917]

        interpo_extrapo = linear_interpolation(d_table,v_table,extrapolation_bc=Line())
        tv = interpo_extrapo(2*r*100)/100 # radius in meters to diameter in cm, then velocity cm->m
    end
    return tv
end

#not working correctly:
# #collision efficiency function
# function collision_efficiency(R1::FT,R2::FT)::FT where FT<:AbstractFloat
#     #Parameterization from Berry 1967
#     #https://doi.org/10.1175/1520-0469(1967)024<0688:CDGBC>2.0.CO;2
#     r = max(R1,R2)*1e6
#     rs = min(R1,R2)*1e6

#     p = rs/r
#     D = (-27)/(r^1.65)
#     E = (-58)/(r^1.9)
#     F = (15/r)^4 +1.13 
#     G = (16.7/r)^8 +1 +0.004*r

#     Y = 1+p+D/(p^F)+E/((1-p)^G)
#     if Y<0
#         Y=0
#     end

#     return Y
# end

#---------------------------------------------------------
# Coalescence Kernels
#---------------------------------------------------------
"""
    hydrodynamic(droplets::droplet_attributes, (j,k)::Tuple{Int,Int}, settings::coag_settings{FT})::FT where FT<:AbstractFloat

Hydronynamic kernel function for droplet coalescence, used to calculate the probabilities of two 
    droplets colliding.
    K(r,r') = π * (r+r')^2 * |v(r)-v(r')| * E(r,r')
    where E(r,r') is the collision efficiency function, (currently not implemented), and v is the terminal velocity of 
    a droplet of radius size r.

# Arguments
- `droplets::droplet_attributes`: The attributes of the droplets.
- `(j,k)::Tuple{Int,Int}`: The indices of the droplets.
- `settings::coag_settings{FT}`: The coagulation settings.

"""
@inline function hydrodynamic(droplets::droplet_attributes, (j,k)::Tuple{Int,Int}, settings::coag_settings{FT})::FT where FT<:AbstractFloat
    Rj, Rk = droplets.R[j], droplets.R[k]
    if settings.hydrodynamic_collision_eff_func == true
        E = collision_efficiency(Rj,Rk)
    else    
        E = 1
    end
    Rsum = (Rj+Rk) # sum of radius is in meters
    vdif = abs(terminal_v(Rj)-terminal_v(Rk))
    return E *π * Rsum^2 *vdif
end

"""
    golovin(droplets, (j, k), settings)

Golovin (or additive) kernel function for droplet coalescence, used to calculate the probabilities of two droplets colliding.
    K(x,x') = b(x+x'),where b is the Golovin kernel coefficient, and x is the volume of the droplet.
Taken from Golovin (1963), this kernel has an analytic solution for the Smulochowski Coagulation Equation.

# Arguments
- `droplets`: droplet attributes
- `(j, k)`: indices of the droplets
- `settings`: coagulation settings

"""
@inline function golovin(droplets::droplet_attributes, (j,k)::Tuple{Int,Int}, settings::coag_settings{FT})::FT where FT<:AbstractFloat
    return settings.golovin_kernel_coeff *(droplets.X[j] + droplets.X[k])# Xsum
end

#---------------------------------------------------------
# PROBABILITIES
#---------------------------------------------------------

"""
    pair_Ps_adaptive!(α::Int, pair::Tuple{Int,Int}, droplets::droplet_attributes, coag_data::coagulation_run, t_max::Vector{FT})

This function calculates the coalescence probability for a given pair of droplets, using adaptive timestepping logic (Bartman et al. 2021). 
It takes the following arguments:

- `α::Int`: The index of the coalescence model.
- `pair::Tuple{Int,Int}`: The indices of the droplets in the pair.
- `droplets::droplet_attributes`: The attributes of the droplets.
- `coag_data::coagulation_run`: The coagulation run data.
- `t_max::Vector{FT}`: The maximum timestep allowed given multiple sampling.

"""
function pair_Ps_adaptive!(α::Int,pair::Tuple{Int,Int}, droplets::droplet_attributes, coag_data::coagulation_run,t_max::Vector{FT},
    t_left::Ref{FT},kernel::Function,settings::coag_settings{FT}) where FT<:AbstractFloat
    if droplets.ξ[pair[1]] < droplets.ξ[pair[2]]
        k = pair[1]
        j = pair[2]
    else
        k = pair[2]
        j = pair[1]
    end

    ξj, ξk = droplets.ξ[j], droplets.ξ[k]

    tmp = ξj * kernel(droplets, pair, settings) * settings.scale / settings.ΔV

    coag_data.pαdt[α] = tmp
    if tmp*t_left[] > coag_data.ϕ[α]
        t_max[α] = (div(ξj, ξk)) / tmp # div gets the floor
    else
        t_max[α] = t_left[]
    end

end


"""
    compute_pαdt!(L::Vector{Tuple{Int,Int}}, droplets::droplet_attributes, coag_data::coagulation_run, kernel::Function, coagsettings::coag_settings{FT}) where FT<:AbstractFloat

Map the probability function over the list of droplet pairs, L, and update the coagulation data in place.

# Arguments
- `L::Vector{Tuple{Int,Int}}`: List of droplet indices to be considered for coalescence.
- `droplets::droplet_attributes`: Droplet attributes.
- `coag_data::coagulation_run`: Coagulation data.
- `kernel::Function`: Coalescence kernel function.
- `coagsettings::coag_settings{FT}`: Coagulation settings.

"""
@inline function compute_pαdt!(L::Vector{Tuple{Int,Int}}, droplets::droplet_attributes,coag_data::coagulation_run,kernel::Function,coagsettings::coag_settings{FT}) where FT<:AbstractFloat
    map(i -> pair_Ps!(i, L[i], droplets,coag_data,kernel, coagsettings), eachindex(coag_data.pαdt))
    coag_data.pαdt .*=  coagsettings.scale * coagsettings.Δt / coagsettings.ΔV
end

@inline function pair_Ps!(α::Int, (j,k)::Tuple{Int,Int}, droplets::droplet_attributes,coag_data::coagulation_run,kernel::Function,coagsettings::coag_settings{FT}) where FT<:AbstractFloat
    coag_data.pαdt[α] = max(droplets.ξ[j], droplets.ξ[k]) * kernel(droplets,(j,k), coagsettings)
end


"""
    adaptive_pαdt!(L::Vector{Tuple{Int,Int}}, droplets::droplet_attributes, coag_data::coagulation_run, t_left::Ref{FT}, kernel::Function, settings::coag_settings{FT}) where FT<:AbstractFloat

Map the probability function over the list of droplet pairs, L, using adaptive timestepping logic and update the coagulation data in place.

# Arguments
- `L::Vector{Tuple{Int,Int}}`: List of droplet indices.
- `droplets::droplet_attributes`: Droplet attributes.
- `coag_data::coagulation_run`: Coagulation data.
- `t_left::Ref{FT}`: Timestep after adaptive substep.
- `kernel::Function`: Coalescence kernel function.
- `settings::coag_settings{FT}`: Coagulation settings.

"""
function adaptive_pαdt!(L::Vector{Tuple{Int,Int}},droplets::droplet_attributes,coag_data::coagulation_run,t_left::Ref{FT},kernel::Function, settings::coag_settings{FT}) where FT<:AbstractFloat
    t_max = Vector{FT}(undef, length(L))
    map(α -> pair_Ps_adaptive!(α,L[α],droplets,coag_data,t_max,t_left,kernel,settings),eachindex(coag_data.pαdt))
    Δtm = minimum(t_max)
    Δt = min(Δtm,t_left[])
    coag_data.pαdt .*= Δt
    t_left[] -= Δt
    return nothing
end


#----------------------------------------------------------
# COALESCENCE
#----------------------------------------------------------

"""
    coalescence_timestep!(run::backend, scheme::scheme_type, droplets::droplet_attributes)

Perform a coalescence timestep for the given droplets using the Superdroplet Method (SDM) 
Shima et al. (2009)
when the lowest multiplicity of superdroplets is less than 1, the 
largest superdroplet is split into two equal parts, as proposed by
Dziekan and Pawlowska (ACP, 2017) https://doi.org/10.5194/acp-17-13509-2017

# Arguments
- `run::backend`: Threading over Linear Sampling option
- `scheme::schemetype`: adaptive or none
- `droplets::droplet_attributes`: The superdroplets.

"""
function coalescence_timestep!(run::Union{Serial, Parallel},scheme::none, droplets::droplet_attributes,
    coag_data::coagulation_run,settings::coag_settings{FT}) where FT<:AbstractFloat
    Ns::Int = settings.Ns
    
    shuffle!(coag_data.I)
    L = [(coag_data.I[l-1], coag_data.I[l]) for l in 2:2:Ns]

    compute_pαdt!(L, droplets,coag_data,settings.kernel, settings)

    rand!(coag_data.ϕ)

    test_pairs!(run,Ns,L,droplets,coag_data)

end 


function coalescence_timestep!(run::Union{Serial, Parallel},scheme::Adaptive,droplets::droplet_attributes{FT},
    coag_data::coagulation_run,settings::coag_settings{FT}) where FT<:AbstractFloat
    Ns = settings.Ns
    t_left = Ref(settings.Δt)

    while t_left[] > 0
        shuffle!(coag_data.I)
        L = [(coag_data.I[l-1], coag_data.I[l]) for l in 2:2:Ns]
        rand!(coag_data.ϕ)

        adaptive_pαdt!(L,droplets,coag_data,t_left,settings.kernel,settings)

        test_pairs!(run,Ns,L,droplets,coag_data)

    end
    return nothing
end 


#----------------------------------------------------------
# UPDATE
#----------------------------------------------------------
"""
    test_pairs!(scheme::backend, Ns::Int, L::Vector{Tuple{Int,Int}}, droplets::droplet_attributes, coag_data::coagulation_run)

Perform the SDM coalescence update for the superdroplets. Update the droplet attributes in place in case of 
    coalescence event.

# Arguments
- `scheme::Union{Serial, Parallel}`: The scheme to use for the update.
- `Ns::Int`: The number of superdroplets.
- `L::Vector{Tuple{Int,Int}}`: The list of droplet pairs.
- `droplets::droplet_attributes`: The droplet attributes.
- `coag_data::coagulation_run`: The coagulation data.
"""

function test_pairs!(scheme::Serial,Ns::Int,L::Vector{Tuple{Int,Int}},droplets::droplet_attributes{FT},coag_data::coagulation_run) where FT<:AbstractFloat
    
    coag_data.lowest_zero[] = false
    for α::Int in 1:div(Ns, 2)
            
        if coag_data.ϕ[α] >= coag_data.pαdt[α]
            continue
        end
        sdm_update!(L[α], α, droplets,coag_data)
    end
    if coag_data.lowest_zero[] == true
        split_highest_multiplicity!(droplets)
    end
end

function test_pairs!(scheme::Parallel,Ns::Int,L::Vector{Tuple{Int,Int}},droplets::droplet_attributes{FT},coag_data::coagulation_run) where FT<:AbstractFloat
    
    coag_data.lowest_zero[] = false
    Threads.@threads for α in 1:div(Ns, 2)
        if coag_data.ϕ[α] >= coag_data.pαdt[α]
            continue
        end
        sdm_update!(L[α], α, droplets,coag_data)
    end
    if coag_data.lowest_zero[] == true
        split_highest_multiplicity!(droplets)
    end
end

"""
    sdm_update!(pair::Tuple{Int,Int}, α::Int, droplets::droplet_attributes, coag_data::coagulation_run)
Perform the SDM coalescence update for the superdroplets when a coalesence event is determined, updating
    the droplet attributes in place.

# Arguments
- `pair::Tuple{Int,Int}`: The pair of droplets to be updated.
- `α::Int`: The index of the coalescence model.
- `droplets::droplet_attributes`: The droplet attributes.
- `coag_data::coagulation_run`: The coagulation data.
"""

function sdm_update!(pair::Tuple{Int,Int},α::Int, droplets::droplet_attributes{FT},coag_data::coagulation_run) where FT<:AbstractFloat

    if droplets.ξ[pair[1]] < droplets.ξ[pair[2]]
        k = pair[1]
        j = pair[2]
    else
        k = pair[2]
        j = pair[1]
    end

    ξj, ξk = droplets.ξ[j], droplets.ξ[k]
    pα =  coag_data.pαdt[α]

    pα_floor::FT = @fastmath floor(pα)
    γ::FT  = if coag_data.ϕ[α] < pα - pα_floor
        pα_floor + 1
    else
        pα_floor
    end

    floor_ξj_div_ξk = floor(ξj / ξk)
    if γ >= floor_ξj_div_ξk
        coag_data.deficit[] += (γ-floor_ξj_div_ξk)*ξk
        γ = floor_ξj_div_ξk
    end
    # γ  = min(γ , floor(ξj / ξk))
    ξ_j_minus_γ_tilde_ξ_k = ξj - γ  * ξk

    if ξ_j_minus_γ_tilde_ξ_k > 0
        droplets.ξ[j] = ξ_j_minus_γ_tilde_ξ_k
        volume = γ * droplets.X[j] + droplets.X[k]
        droplets.R[k] = volume_to_radius(volume)
        droplets.X[k] = volume 
    elseif ξ_j_minus_γ_tilde_ξ_k == 0
        half_ξ_k = floor(ξk / 2)
        droplets.ξ[j] = half_ξ_k
        droplets.ξ[k] -= half_ξ_k

        volume = γ *droplets.X[j] + droplets.X[k]
        droplets.R[j] = droplets.R[k] = volume_to_radius(volume)
        droplets.X[k] = droplets.X[j] = volume
        if half_ξ_k == 0
            coag_data.lowest_zero[] = true
        end
    else
        println("nooooo")
    end
    return nothing
end

"""
    split_highest_multiplicity!(droplets::droplet_attributes{FT}) where FT<:AbstractFloat
Split the superdroplet with the highest multiplicity into two equal parts, as proposed by
    Dziekan and Pawlowska (ACP, 2017) https://doi.org/10.5194/acp-17-13509-2017
    This function is called when the lowest multiplicity of superdroplets is less than 1.
# Arguments
- `droplets::droplet_attributes{FT}`: The droplet attributes.
"""


function split_highest_multiplicity!(droplets::droplet_attributes{FT}) where FT<:AbstractFloat
    if maximum(droplets.ξ) > 1
        while (minimum(droplets.ξ) <= 0 && maximum(droplets.ξ) > 1)
            argmin_i = argmin(droplets.ξ)
            argmax_i = argmax(droplets.ξ)
            droplets.ξ[argmin_i] = floor(droplets.ξ[argmax_i]/2)
            droplets.R[argmin_i] = droplets.R[argmax_i]
            droplets.X[argmin_i] = droplets.X[argmax_i]
            # M[argmin_i] = M[argmax_i]
            droplets.ξ[argmax_i] -= floor(droplets.ξ[argmax_i]/2)
        end
    elseif (maximum(droplets.ξ) <= 1)

        println("Highest superdroplet cannot be split")
        #right now, break the model until this situation gets handled
        if (maximum(droplets.ξ) < 1)
            error("Highest and Lowest Superdroplet have ξ==0")
        end

        #Later:remove superdroplet.. how to handle between cells?

        # # Cannot split highest multiplicity superdroplet, have to remove superdroplet from system
        # println("Superdroplet ", argmin(ξ), " has multiplicity of ", ξ[argmin(ξ)], ", removing from system")
        # deleteat!(R,argmin(ξ))
        # # deleteat!(M,argmin(ξ))
        # deleteat!(X,argmin(ξ))
        # deleteat!(ξ,argmin(ξ))
        # # Ns=Ns-1
    end
    return nothing
end




