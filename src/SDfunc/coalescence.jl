#---------------------------------------------------------
#Collision Coalescence Functions
#---------------------------------------------------------
#currently, this collision-coalescence is only set up to handle radius, 
#multiplicity, volume, and solute mass as superdroplet attributes
using Distributions
using Combinatorics
using Random
export init_ξ_const,init_logarithmic,terminal_v,collision_efficiency,hydrodynamic,golovin,calc_Ps,coalescence_timestep!,all_or_nothing!,coalescence_unittest_graph!
export Serial, Parallel,Adaptive,none,coag_settings,droplets_allocations


struct Serial end
struct Parallel end
struct Adaptive end
struct none end

Base.@kwdef struct coag_settings{FT<:AbstractFloat}
    Δt::FT = FT(1.0)
    ΔV::FT = FT(1e6)
    Ns::Int # number of superdroplets
    scale::FT
    R_min::FT = FT(1e-9)
    R_max::FT = FT(1e-3)
    golovin_kernel_coeff::FT = FT(1.5e3)
    hydrodynamic_collision_eff_func::Bool = false
    kernel::Function = golovin # golovin, hydrodynamic
    n0::FT = FT(2^23) # initial droplet concentration
    R0::FT = FT(30.531e-6) # initial radius
end

struct droplets_allocations{FT<:AbstractFloat}
    ξ::Vector{Int}
    R::Vector{FT}
    X::Vector{FT}
    I::Vector{Int}
    pαdt::Vector{FT}
    ϕ::Vector{FT}
end

#functions:

#INITIALIZATION
    #init_ξ_const(Ns,ΔV,n0,R0),         returns ξstart, Rstart, Xstart
    #init_ξ_const(Ns,ΔV,n0,R0,M0),      returns ξstart, Rstart, Xstart, Mstart

#KERNEL FUNCTIONS
    #terminal_v(r),                     returns tv
    #collision_efficiency(R1,R2),       returns efficiency
    #hydrodynamic(droplet1,droplet2),   returns hydrodynamic kernel calc
    #hydrodynamic(R1,R2,X,Xs),          returns hydrodynamic kernel calc
    #golovin(droplet1,droplet2),        returns golovin kernel calc
    #golovin(R1,R2,X,Xs),               returns golovin kernel calc

#SDM FUNCTIONS
    #calc_Ps(droplet1,droplet2,Δt,ΔV ),                       returns Ps
    #calc_Ps(R1,R2,X,Xs,Δt,ΔV,ξj,ξk ),                        returns Ps
    #coalescence_timestep!(droplets,Ns,Δt,ΔV ),               returns droplets
    #coalescence_timestep!(ξ,R,X,Ns,Δt,ΔV ),                  returns ξ,R,X
    #coalescence_timestep!(ξ,R,M,X,Ns,Δt,ΔV ),                returns ξ,R,M,X
    #all_or_nothing!(R1,R2,X1,X2,ξ1,ξ2,M1,M2,scale,Δt,ΔV ),   returns R1,R2,X1,X2,ξ1,ξ2,M1,M2
    #all_or_nothing!(R1,R2,X1,X2,ξ1,ξ2,scale,Δt,ΔV ),         returns R1,R2,X1,X2,ξ1,ξ2
    #all_or_nothing!(droplet1,droplet2,scale,Δt,ΔV ),         returns droplet1,droplet2


#SMALL ALPHA STUDY
    #coalescence_timestep_small_alpha!(droplets,Ns,y,Δt,ΔV ), returns droplets

#UNIT TEST RUN, RETURNS GRAPH 0:1200:3600 seconds
    #coalescence_unittest_graph!(droplets,Ns,Δt,ΔV;smooth=true,label=true,kernel=golovin,
        #radius_bins_edges=10 .^ range(log10(10*1e-6), log10(5e3*1e-6), length=128))
    #coalescence_unittest_graph!(ξstart,Rstart,Xstart,Ns,Δt,ΔV,smooth=smooth,label=true,kernel=golovin,
        #radius_bins_edges=10 .^ range(log10(10*1e-6), log10(5e3*1e-6), length=128))


#---------------------------------------------------------
# INITIALIZATION
#---------------------------------------------------------
# init_ξ_const initializes arrays based on the constant-multiplicity initialization method
# described in Shima et al. (2009) https://doi.org/10.5194/acp-9-4491-2009
# n0 is the number concentration of droplets, R0 an initial radius of the droplets
# ΔV is the volume of the domain, Ns is the number of superdroplets
# the second option adds mass of solute as an attribute, M0 is the initial mass of solute
# the function draws from an exponential distribution
# MOVE

function init_ξ_const(settings::coag_settings{FT}) where FT<:AbstractFloat
    Ns = settings.Ns
    ΔV = settings.ΔV
    n0 = settings.n0
    R0 = settings.R0
    ξstart::Vector{Int} = (n0*ΔV/Ns*ones(Ns))
    # R0 = Float64(30.531e-6) # meters
    X0 = Float64(4*π/3*R0^3) # initial volume m3    
    Xstart::Vector{FT} = (rand(Exponential(X0), Ns))
    Rstart::Vector{FT} = ((3 .*Xstart./(4*π)).^(1/3))
    return ξstart, Rstart, Xstart
end





function init_logarithmic(settings::coag_settings{FT})where FT<:AbstractFloat
    Ns = settings.Ns
    ΔV = settings.ΔV
    n0 = settings.n0
    R0 = settings.R0
    R_min = settings.R_min
    R_max = settings.R_max

    X0 = (4/3) * π * R0^3
    exp_dist = Exponential(X0)
    boundaries_found = false
    while boundaries_found == false
        radius_bins = 10 .^ range(log10(R_min), log10(R_max), length=Ns+1)
        # Calculate the volume for each bin
        volumes = (4/3) * π .* radius_bins.^3
        vmin = (4/3) * π .* radius_bins[1]^3
        vmax = (4/3) * π .* radius_bins[end-1]^3
        
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
    volumes::Vector{FT} = (4/3) * π .* sd_radii.^3
    pdf_values = pdf.(exp_dist, volumes)
    bin_widths_new = diff(radius_bins_new)
    dvdr = 4 * π .* sd_radii.^2
    multiplicities = pdf_values .* bin_widths_new .* dvdr * (n0*ΔV)
    ξstart::Vector{Int} = floor.(multiplicities.+0.5)

return ξstart, sd_radii, volumes
end

function init_uniform_sd(settings::coag_settings{FT})where FT<:AbstractFloat
    Ns = settings.Ns
    ΔV = settings.ΔV
    n0 = settings.n0
    R0 = settings.R0
    R_min = settings.R_min
    R_max = settings.R_max

    X0 = (4/3) * π * R0^3
    exp_dist = Exponential(X0)
    boundaries_found = false
    while boundaries_found == false
        radius_bins = range(R_min, R_max, length=Ns+1)
        # Calculate the volume for each bin
        # volumes = (4/3) * π .* radius_bins.^3
        vmin = (4/3) * π * radius_bins[1]^3
        vmax = (4/3) * π * radius_bins[end-1]^3
        
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
    volumes::Vector{FT} = (4/3) * π .* sd_radii.^3
    pdf_values = pdf.(exp_dist, volumes)
    bin_widths_new = diff(radius_bins_new)
    dvdr = 4 * π .* sd_radii.^2
    multiplicities = pdf_values .* bin_widths_new .* dvdr * (n0*ΔV)
    ξstart::Vector{Int} = floor.(multiplicities.+0.5)

return ξstart, sd_radii, volumes
end



#---------------------------------------------------------
# HYDRODYNAMIC KERNEL
#---------------------------------------------------------

# terminal velocity of droplets
function terminal_v(r)::FT # terminal velocity 
    # Tables from
    # THE TERMINAL VELOCITY OF FALL FOR WATER DROPLETS IN STAGNANT AIR
    # Gunn, R., Kinzer, G. D. (1949), https://doi.org/10.1175/1520-0469(1949)006<0243:TTVOFF>2.0.CO;2

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

#collision efficiency function
function collision_efficiency(R1,R2)
    #Parameterization from Berry 1967
    #https://doi.org/10.1175/1520-0469(1967)024<0688:CDGBC>2.0.CO;2
    r = max(R1,R2)*1e6
    rs = min(R1,R2)*1e6

    p = rs/r
    D = (-27)/(r^1.65)
    E = (-58)/(r^1.9)
    F = (15/r)^4 +1.13 
    G = (16.7/r)^8 +1 +0.004*r

    Y = 1+p+D/p^F+E/(1-p)^G
    if Y<0
        Y=0
    end

    return Y
end

###########################################################
# The hydrodynamic kernel function two methods: one for superdroplet 
# structures and one for two radius values -- it takes dummy volumes so 
# that the call could be generalized in the coalescence_timestep function


@inline function hydrodynamic(droplets::droplets_allocations, (j,k)::Tuple{Int,Int}, settings::coag_settings{FT})::FT where FT<:AbstractFloat
    Rj, Rk = droplets.R[j], droplets.R[k]
    E = collision_efficiency(Rj,Rk)
    Rsum = (Rj+Rk) # sum of radius is in meters
    vdif = abs(terminal_v(Rj)-terminal_v(Rk))
    return E *π * Rsum^2 *vdif
end

#---------------------------------------------------------
# GOLOVIN KERNEL (Golovin,1963)
#---------------------------------------------------------
###########################################################
# The golovin kernel function two methods: one for superdroplet 
# structures and one for two volume values -- it takes dummy radius values so 
# that the call could be generalized in the coalescence_timestep function



@inline function golovin(droplets::droplets_allocations, (j,k)::Tuple{Int,Int}, settings::coag_settings{FT})::FT where FT<:AbstractFloat
    return settings.golovin_kernel_coeff *(droplets.X[j] + droplets.X[k])# Xsum
end

#---------------------------------------------------------
# PROBABILITIES
#---------------------------------------------------------
#calc_PS calculates the probability of two superdroplets colliding in a 
#given timestep and volume,
#kernel is a kwarg that can be set to either golovin or hydrodynamic
#the function has two methods: one for superdroplet structures and one for
#free radius, volume, and multiplicity values



function pair_Ps_adaptive!(α::Int,pair::Tuple{Int,Int}, droplets::droplets_allocations, t_max,kernel::Function,settings::coag_settings{FT}) where FT<:AbstractFloat
    if droplets.ξ[pair[1]] < droplets.ξ[pair[2]]
        k = pair[1]
        j = pair[2]
    else
        k = pair[2]
        j = pair[1]
    end

    ξj, ξk = droplets.ξ[j], droplets.ξ[k]

    tmp = ξj * kernel(droplets, pair, settings) * settings.scale / settings.ΔV
    droplets.pαdt[α] = tmp
    t_max[α] = (div(ξj, ξk)) / tmp # div gets the floor 

end

#----------------------------------------------------------
# COALESCENCE
#----------------------------------------------------------
#the coalescence_timestep function takes superdroplets and runs them
#through an All-or-nothing update timestep based on the SDM method 
#described in Shima et al. (2009)
#when the lowest multiplicity of superdroplets is less than 1, the 
#largest superdroplet is split into two equal parts, as proposed by
#Dziekan and Pawlowska (ACP, 2017) https://doi.org/10.5194/acp-17-13509-2017

#currently all can be parallelized on CPUs with Threads.@threads


@inline function compute_pαdt!(L::Vector{Tuple{Int,Int}}, droplets::droplets_allocations, kernel::Function,coagsettings::coag_settings{FT}) where FT<:AbstractFloat
    map(i -> pair_Ps!(i, L[i], droplets, kernel, coagsettings), eachindex(droplets.pαdt))
    droplets.pαdt .*=  coagsettings.scale * coagsettings.Δt / coagsettings.ΔV
end

@inline function pair_Ps!(α::Int, (j,k)::Tuple{Int,Int}, droplets::droplets_allocations, kernel::Function,coagsettings::coag_settings{FT}) where FT<:AbstractFloat
    droplets.pαdt[α] = max(droplets.ξ[j], droplets.ξ[k]) * kernel(droplets,(j,k), coagsettings)
end

function coalescence_timestep!(run::Serial,scheme::none, droplets::droplets_allocations,
    settings::coag_settings{FT}) where FT<:AbstractFloat
    Ns::Int = settings.Ns
    
    shuffle!(droplets.I)
    L = [(droplets.I[l-1], droplets.I[l]) for l in 2:2:Ns]

    compute_pαdt!(L, droplets,settings.kernel, settings)

    rand!(droplets.ϕ)

    for α::Int in 1:div(Ns, 2)
        
        if droplets.ϕ[α] >= droplets.pαdt[α]
            continue
        end
        sdm_update!(L[α], α, droplets)

    end


end 


function coalescence_timestep!(run::Parallel, scheme::none,droplets::droplets_allocations,
    settings::coag_settings{FT}) where FT<:AbstractFloat

    Ns = settings.Ns
    # this is quicker outside of the parallelization....
    # probably can change???
    L = [(droplets.I[l-1], droplets.I[l]) for l in 2:2:Ns]
    compute_pαdt!(L, droplets,settings.kernel, settings)(L, droplets,settings.kernel, settings)
    rand!(droplets.ϕ)

    Threads.@threads for α in 1:div(Ns, 2)
        if droplets.ϕ[α] >= droplets.pαdt[α]
            continue
        end
        sdm_update!(L[α], α, droplets)
    end
    return nothing
end


function adaptive_pαdt!(L::Vector{Tuple{Int,Int}},droplets::droplets_allocations,t_left::Ref{FT},kernel::Function, settings::coag_settings{FT}) where FT<:AbstractFloat
    t_max = Vector{FT}(undef, length(L))
    map(α -> pair_Ps_adaptive!(α,L[α],droplets,t_max,kernel,settings),eachindex(droplets.pαdt))
    Δtm = minimum(t_max)
    Δt = min(Δtm,t_left[])
    droplets.pαdt .*= Δt
    t_left[] -= Δt
    return nothing
end


function coalescence_timestep!(run::Serial,scheme::Adaptive,droplets::droplets_allocations{FT},
    settings::coag_settings{FT}) where FT<:AbstractFloat
    Ns = settings.Ns
    t_left = Ref(settings.Δt)

    while t_left[] > 0
        shuffle!(droplets.I)
        L = [(droplets.I[l-1], droplets.I[l]) for l in 2:2:Ns]
        rand!(droplets.ϕ)

        adaptive_pαdt!(L,droplets,t_left,settings.kernel,settings)
        
        for α::Int in 1:div(Ns, 2)
        
            if droplets.ϕ[α] >= droplets.pαdt[α]
                continue
            end
            sdm_update!(L[α], α, droplets)
    
        end

        #this takes a lot of searching.. could we put the failure condition somewhere else?
        if minimum(droplets.ξ) <= 0
            println("still happening")
            split_highest_multiplicity!(droplets.ξ,droplets.R,droplets.X)
        end
    end
    return nothing
end 

function coalescence_timestep!(run::Parallel,scheme::Adaptive,droplets::droplets_allocations{FT},
    settings::coag_settings{FT}) where FT<:AbstractFloat
    Ns = settings.Ns
    t_left = Ref(settings.Δt)

    while t_left[] > 0
        shuffle!(droplets.I)
        L = [(droplets.I[l-1], droplets.I[l]) for l in 2:2:Ns]
        rand!(droplets.ϕ)

        #surely there is a way to parallelize this
        adaptive_pαdt!(L,droplets,t_left,settings.kernel,settings)
        
        Threads.@threads for α in 1:div(Ns, 2)
            if droplets.ϕ[α] >= droplets.pαdt[α]
                continue
            end
            sdm_update!(L[α], α, droplets)
        end

        #this takes a lot of searching.. could we put the failure condition somewhere else?
        if minimum(droplets.ξ) <= 0
            println("still happening")
            split_highest_multiplicity!(droplets.ξ,droplets.R,droplets.X)
        end
    end
    return nothing
end 


function sdm_update!(pair::Tuple{Int,Int},α::Int, droplets::droplets_allocations{FT}) where FT<:AbstractFloat

    if droplets.ξ[pair[1]] < droplets.ξ[pair[2]]
        k = pair[1]
        j = pair[2]
    else
        k = pair[2]
        j = pair[1]
    end

    ξj, ξk = droplets.ξ[j], droplets.ξ[k]
    pα =  droplets.pαdt[α]

    pα_floor::FT = @fastmath floor(pα)
    γ::FT  = if droplets.ϕ[α] < pα - pα_floor
        pα_floor + 1
    else
        pα_floor
    end

    γ  = min(γ , floor(ξj / ξk))
    ξ_j_minus_γ_tilde_ξ_k = ξj - γ  * ξk

    if ξ_j_minus_γ_tilde_ξ_k > 0
        droplets.ξ[j] = ξ_j_minus_γ_tilde_ξ_k
        R_k_cubed = γ * droplets.R[j]^3 + droplets.R[k]^3
        droplets.R[k] = R_k_cubed^(1/3)
        droplets.X[k] = 4/3 * π * R_k_cubed
    elseif ξ_j_minus_γ_tilde_ξ_k == 0
        half_ξ_k = floor(ξk / 2)
        droplets.ξ[j] = half_ξ_k
        droplets.ξ[k] -= half_ξ_k
        R_k_cubed = γ* droplets.R[j]^3 + droplets.R[k]^3
        droplets.R[j] = droplets.R[k] = R_k_cubed^(1/3)
        droplets.X[k] = droplets.X[j] = 4/3 * π * R_k_cubed
    else
        println("nooooo")
    end
    return nothing
end



function split_highest_multiplicity!(ξ,R,X)
    if maximum(ξ) > 1
        while (minimum(ξ) <= 0 && maximum(ξ) > 1)
            println("Superdroplet ", argmin(ξ), " has multiplicity of ", ξ[argmin(ξ)])
            # Split superdroplet with highest multiplicity, half goes to argmin(ξ) and half stays
            # Idea from Dziekan and Pawlowska (ACP, 2017) https://doi.org/10.5194/acp-17-13509-2017
            argmin_i = argmin(ξ)
            argmax_i = argmax(ξ)
            ξ[argmin_i] = floor(ξ[argmax_i]/2)
            R[argmin_i] = R[argmax_i]
            X[argmin_i] = X[argmax_i]
            # M[argmin_i] = M[argmax_i]
            ξ[argmax_i] = floor(ξ[argmax_i]/2)
        end
    elseif (maximum(ξ) <= 1)
        #right now, break the model until this situation gets handled
        error("Highest and Lowest Superdroplet have ξ==0")

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




