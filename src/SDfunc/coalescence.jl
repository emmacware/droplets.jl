#---------------------------------------------------------
#Collision Coalescence Functions
#---------------------------------------------------------
#currently, this collision-coalescence is only set up to handle radius, 
#multiplicity, volume, and solute mass as superdroplet attributes
using Distributions
using Combinatorics
using Random
using StaticArrays
export init_ξ_const,init_logarithmic,terminal_v,collision_efficiency,hydrodynamic,golovin,calc_Ps,coalescence_timestep!,all_or_nothing!,coalescence_unittest_graph!
export Serial, Parallel,Adaptive,none,coag_settings


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


# function init_ξ_const(Ns,ΔV,n0,R0,M0)
#     ξstart = n0*ΔV/Ns*ones(Float64,Ns)
#     # R0 = Float64(30.531e-6) # meters
#     X0 = Float64(4*π/3*R0^3) # initial volume m3
#     Xstart = rand(Exponential(X0), Ns)
#     Rstart = (3 .*Xstart./(4*π)).^(1/3)
#     Mstart = rand(Exponential(M0), Ns)
#     return ξstart, Rstart, Xstart,Mstart
# end


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

function hydrodynamic(droplet1,droplet2,settings::coag_settings{FT})::FT where FT<:AbstractFloat
    E = collision_efficiency(droplet1.R,droplet2.R)
    Rsum = (droplet1.R+droplet2.R) # sum of radius is in meters
    vdif = abs(terminal_v(droplet1.R)-terminal_v(droplet2.R))
    return E*π*Rsum^2*vdif
end

function hydrodynamic(R,X,j,k,settings::coag_settings{FT})::FT where FT<:AbstractFloat
    E = collision_efficiency(R[j],R[k])
    Rsum = (R[j]+R[k]) # sum of radius is in meters
    vdif = abs(terminal_v(R[j])-terminal_v(R[k]))
    return E *π * Rsum^2 *vdif
end

#---------------------------------------------------------
# GOLOVIN KERNEL (Golovin,1963)
#---------------------------------------------------------
###########################################################
# The golovin kernel function two methods: one for superdroplet 
# structures and one for two volume values -- it takes dummy radius values so 
# that the call could be generalized in the coalescence_timestep function

# function golovin(droplet1,droplet2,coag_settings)
#     b = 1.5e3 # seconds^-1
#     Xsum = droplet1.X+droplet2.X # m3
#     return b*Xsum
# end

# function golovin(R,X,j,k,settings::coag_settings{FT})where FT<:AbstractFloat
#     # b = 1.5e3 # seconds^-1
#     Xsum = X[j]+X[k] # m3
#     return settings.golovin_kernel_coeff*Xsum
# end

# @inline function golovin(R::Vector{FT}, X::Vector{FT}, j::Int, k::Int, settings::coag_settings{FT})::FT where FT<:AbstractFloat
#     Xsum = X[j] + X[k] # m3
#     return settings.golovin_kernel_coeff * Xsum
# end

@inline function golovin(X::Vector{FT}, (j,k)::Tuple{Int,Int}, settings::coag_settings{FT})::FT where FT<:AbstractFloat
    # Xsum::FT = X[j] + X[k] # m3
    return settings.golovin_kernel_coeff *(X[j] + X[k])# Xsum
    # return FT(1500)*Xsum
end

#---------------------------------------------------------
# PROBABILITIES
#---------------------------------------------------------
#calc_PS calculates the probability of two superdroplets colliding in a 
#given timestep and volume,
#kernel is a kwarg that can be set to either golovin or hydrodynamic
#the function has two methods: one for superdroplet structures and one for
#free radius, volume, and multiplicity values


function calc_Ps(droplet1,droplet2,Δt,ΔV,settings::coag_settings{FT} )::FT where FT<:AbstractFloat
    Pjk = settings.kernel(droplet1,droplet2,settings)*Δt/ΔV
    Ps = max(droplet1.ξ,droplet2.ξ)*Pjk 
    return Ps
end

# function calc_Ps(R,X,ξ,j,k,settings::coag_settings{FT} )where FT<:AbstractFloat
#     return max(ξ[j],ξ[k])*settings.kernel(R,X,j,k,settings)*settings.Δt/settings.ΔV
# end

function calc_Ps(R::Vector{FT}, X::Vector{FT}, ξ::Vector{Int}, j::Int, k::Int, settings::coag_settings{FT})::FT where FT<:AbstractFloat
    return max(ξ[j], ξ[k]) * settings.kernel(R, X, j, k, settings) * settings.Δt / settings.ΔV
end


# function pair_Ps_adaptive(pair,ξ,R,X,settings::coag_settings{FT})where FT<:AbstractFloat
#     ξj = max(ξ[pair[1]],ξ[pair[2]])
#     ξk = min(ξ[pair[1]],ξ[pair[2]])
#     # if ξk==0
#     #     println("oop")
#     # end
#     pαdt = ξj*settings.kernel(R,X,pair[1],pair[2],settings)*settings.scale/settings.ΔV
#     Δtmax = (div(ξj,ξk))/pαdt #div gets the floor 
#     return pαdt,Δtmax
# end

function pair_Ps_adaptive(pair::Tuple{Int,Int}, ξ::Vector{Int}, X::Vector{FT}, settings::coag_settings{FT})::Tuple{FT,FT} where FT<:AbstractFloat
    ξk::FT, ξj::FT = if ξ[pair[1]] < ξ[pair[2]]
        ξ[pair[1]], ξ[pair[2]]
    else
        ξ[pair[2]], ξ[pair[1]]
    end
    # # ξj::FT = max(ξ[pair[1]], ξ[pair[2]])
    # ξk::FT = min(ξ[pair[1]], ξ[pair[2]])
    
    pαdt::FT = ξj * settings.kernel(X, pair, settings) * settings.scale / settings.ΔV
    Δtmax::FT = (div(ξj, ξk)) / pαdt # div gets the floor 
    
    return pαdt, Δtmax
end


# function pair_Ps(pair::Tuple{Int,Int}, ξ::Vector{Int}, R::Vector{FT}, X::Vector{FT}, settings::coag_settings{FT})::FT where FT<:AbstractFloat
#     return max(ξ[pair[1]], ξ[pair[2]]) * settings.kernel(R, X, pair[1], pair[2], settings)
# end

# @inline function pair_Ps((j,k)::Tuple{Int,Int}, ξ::Vector{Int}, R::Vector{FT}, X::Vector{FT}, settings::coag_settings{FT})::FT where FT<:AbstractFloat
#     return max(ξ[j], ξ[k]) * settings.kernel(R, X, j, k, settings)
# end
@inline function pair_Ps((j,k)::Tuple{Int,Int}, ξ::Vector{Int}, X::Vector{FT}, settings::coag_settings{FT})::FT where FT<:AbstractFloat
    return max(ξ[j], ξ[k]) * settings.kernel(X,(j,k), settings)
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




# @inline function compute_pαdt(L::Vector{Tuple{Int,Int}}, ξ::Vector{Int}, R::Vector{FT}, X::Vector{FT}, coagsettings::coag_settings{FT})::Vector{FT} where FT<:AbstractFloat
#     return map(pair -> pair_Ps(pair, ξ, R, X, coagsettings), L) .* coagsettings.scale * coagsettings.Δt / coagsettings.ΔV
# end
@inline function compute_pαdt(L::Vector{Tuple{Int,Int}}, ξ::Vector{Int}, X::Vector{FT}, coagsettings::coag_settings{FT})::Vector{FT} where FT<:AbstractFloat
    return map(pair -> pair_Ps(pair, ξ, X, coagsettings), L) .* coagsettings.scale * coagsettings.Δt / coagsettings.ΔV
end

# #####
function coalescence_timestep!(run::Serial,scheme::none, ξ::Vector{Int}, R::Vector{FT}, 
    X::Vector{FT}, I::Vector{Int},ϕ::Vector{FT},
    settings::coag_settings{FT}) where FT<:AbstractFloat
    Ns::Int = settings.Ns
    
    shuffle!(I)
    L::Vector{Tuple{Int, Int}} = [(I[l-1], I[l]) for l in 2:2:Ns]

    # pαdt = map(pair -> pair_Ps(pair,ξ,R,X,settings),L).*settings.scale*settings.Δt/settings.ΔV
    pαdt = compute_pαdt(L, ξ, X, settings)
    rand!(ϕ)
    # ϕ = rand(FT, div(settings.Ns, 2))

    # for (α, pair) in enumerate(L)
    for α::Int in 1:div(Ns, 2)
        
        if ϕ[α] >= pαdt[α]
            continue
        end
        sdm_update!(L[α],ξ,R,X,ϕ[α], pαdt[α])
    end

end 

function construct_L(I::Vector{Int}, Ns::Int)::Vector{Tuple{Int, Int}}
    L = Vector{Tuple{Int, Int}}(undef, div(Ns, 2))
    shuffle!(I)
    Threads.@threads for l in 2:2:Ns
        L[div(l, 2)] = (I[l-1], I[l])
    end
    return L
end

function coalescence_timestep!(run::Parallel, scheme::none,ξ::Vector{Int}, R::Vector{FT}, 
    X::Vector{FT}, I::Vector{Int},ϕ::Vector{FT},
    settings::coag_settings{FT}) where FT<:AbstractFloat

    Ns = settings.Ns
    L = construct_L(I, Ns)

    Threads.@threads for α in 1:div(Ns, 2)
        pαdt = pair_Ps(L[α], ξ, X, settings) * settings.scale * settings.Δt / settings.ΔV
        ϕ[α] = rand(FT)
        # rand!(ϕ)
        if ϕ[α] >= pαdt
            continue
        end
        sdm_update!(L[α], ξ, R, X, ϕ[α], pαdt)
    end
end

function adaptive_pαdt(L::Vector{Tuple{Int,Int}}, ξ::Vector{Int}, X::Vector{FT},t_left::FT, settings::coag_settings{FT}) where FT<:AbstractFloat
    deficit_opt = map(pair -> pair_Ps_adaptive(pair,ξ,X,settings),L)
    Δtm = minimum([last(pair) for pair in deficit_opt])
    Δt = min(Δtm,t_left)
    pαdt = [first(pair) for pair in deficit_opt].*Δt
    return pαdt, Δt
end

function adaptive_pαdt_parallel(L::Vector{Tuple{Int,Int}}, ξ::Vector{Int}, X::Vector{FT},t_left::FT, settings::coag_settings{FT}) where FT<:AbstractFloat
    pd = Vector{FT}(undef, length(L))
    tm = Vector{FT}(undef, length(L))
    Threads.@threads for α in 1:length(L)
        pd[α],tm[α] = pair_Ps_adaptive(L[α],ξ,X,settings)
    end
    Δtm = minimum(tm)
    Δt = min(Δtm,t_left)
    pαdt = pd.*Δt
    return pαdt, Δt
end

# #####
function coalescence_timestep!(run::Serial,scheme::Adaptive,ξ::Vector{Int}, R::Vector{FT}, 
    X::Vector{FT}, I::Vector{Int},ϕ::Vector{FT},
    settings::coag_settings{FT}) where FT<:AbstractFloat
    Ns = settings.Ns
    t_left = settings.Δt #+ 0

    while t_left > 0
        shuffle!(I)
        L::Vector{Tuple{Int, Int}} = [(I[l-1], I[l]) for l in 2:2:Ns]

        rand!(ϕ)
        # deficit_opt = map(pair -> pair_Ps_adaptive(pair,ξ,X,settings),L)
        # Δtm = minimum([last(pair) for pair in deficit_opt])
        # Δt = min(Δtm,t_left)
        # pαdt = [first(pair) for pair in deficit_opt].*Δt
        pαdt, Δt = adaptive_pαdt(L, ξ, X, t_left ,settings)
        
        # for (α, pair) in enumerate(L)
            
        #     if ϕ[α] >= pαdt[α]
        #         continue
        #     end
        #     sdm_update!(pair,ξ,R,X,ϕ[α], pαdt[α])
        # end

        for α::Int in 1:div(Ns, 2)
        
            if ϕ[α] >= pαdt[α]
                continue
            end
            sdm_update!(L[α],ξ,R,X,ϕ[α], pαdt[α])
        end

        #this takes a lot of searching.. could we put the failure condition somewhere else?
        if minimum(ξ) <= 0
            split_highest_multiplicity!(ξ,R,X)
        end

        t_left -= Δt
    end
end 

#still dev
# function coalescence_timestep!(run::Parallel,scheme::Adaptive,ξ::Vector{Int}, R::Vector{FT}, 
#     X::Vector{FT}, I::Vector{Int},ϕ::Vector{FT},
#     settings::coag_settings{FT}) where FT<:AbstractFloat

#     Ns = settings.Ns
#     t_left = settings.Δt

#     while t_left > 0
#         L = construct_L(I, Ns)

#         # ϕ = rand(div(Ns, 2))

#         pαdt,Δt = adaptive_pαdt_parallel(L, ξ, X, t_left, settings)
        
#         Threads.@threads for α=1:Int(floor(Ns/2))
            
#             if ϕ[α] >= pαdt[α]
#                 continue
#             end
#             sdm_update!(L[α],ξ,R,X,ϕ[α], pαdt[α])
#         end

#         #this takes a lot of searching.. could we put the failure condition somewhere else?
#         if minimum(ξ) <= 0
#             split_highest_multiplicity!(ξ,R,X)
#         end

#         t_left -= Δt
#     end
# end 



function sdm_update!(pair::Tuple{Int,Int}, ξ::Vector{Int}, R::Vector{FT}, X::Vector{FT}, ϕ::FT, pα::FT) where FT<:AbstractFloat
        k, j = if ξ[pair[1]] < ξ[pair[2]]
        pair[1], pair[2]
    else
        pair[2], pair[1]
    end

    pα_floor::FT = floor(pα)
    γ::FT = if ϕ < pα - pα_floor
        pα_floor + 1
    else
        pα_floor
    end

    γ_tilde::FT = min(γ, floor(ξ[j] / ξ[k]))
    ξ_j_minus_γ_tilde_ξ_k = ξ[j] - γ_tilde * ξ[k]

    if ξ_j_minus_γ_tilde_ξ_k > 0
        ξ[j] = ξ_j_minus_γ_tilde_ξ_k
        R_k_cubed = γ_tilde * R[j]^3 + R[k]^3
        R[k] = R_k_cubed^(1/3)
        X[k] = 4/3 * π * R[k]^3
    elseif ξ_j_minus_γ_tilde_ξ_k == 0
        half_ξ_k = floor(ξ[k] / 2)
        ξ[j] = half_ξ_k
        ξ[k] -= half_ξ_k
        R_k_cubed = γ_tilde * R[j]^3 + R[k]^3
        R[j] = R[k] = R_k_cubed^(1/3)
        X[k] = X[j] = 4/3 * π * R[k]^3
    else
        println("nooooo")
    end
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
end










######### OLD DROPLETS structure

# function coalescence_timestep!(droplets,Ns,Δt,ΔV )
#     shuffle!(I)
#     L= [(I[l-1],I[l]) for l=2:2:length(I)]
#     scale = Ns*(Ns-1)/2/(Ns/2)

#     Threads.@threads for α=1:Int(floor(Ns/2))
#     # for α in 1:Ns/2

#         jj=L[α][1] #struggled with choosing j>k so its hard coded
#         kk=L[α][2]

#         if droplets[jj].ξ<droplets[kk].ξ
#             k=jj
#             j=kk
#         else
#             j=jj
#             k=kk
#         end			

#         ϕ = rand()
#         pα = scale*calc_Ps(droplets[j],droplets[k],Δt,ΔV,kernel=kernel)


#         if ϕ < pα - floor(pα)
#             γ = floor(pα)+1
#         elseif ϕ >= pα - floor(pα)
#             γ = floor(pα)
#         end

#         if γ == 0
#             continue
#         end

#         γ_tilde = min(γ,floor(droplets[j].ξ/droplets[k].ξ))
#         # if γ_tilde == floor(ξ[j]/ξ[k])
#         #     println("limit exceeded")
#         # end

#         if droplets[j].ξ - γ_tilde*droplets[k].ξ > 0
#             droplets[j].ξ = droplets[j].ξ-γ_tilde*droplets[k].ξ
#             droplets[k].R = (γ_tilde*droplets[j].R^3+droplets[k].R^3)^(1/3)
#             droplets[k].X = 4/3*π * droplets[k].R^3
#             droplets[k].M = γ_tilde*droplets[j].M+droplets[k].M

#         elseif droplets[j].ξ - γ_tilde*droplets[k].ξ == 0
#             droplets[j].ξ = floor(droplets[k].ξ/2)
#             droplets[k].ξ = droplets[k].ξ - floor(droplets[k].ξ/2)
#             droplets[j].R = droplets[k].R = (γ_tilde*droplets[j].R^3+droplets[k].R^3)^(1/3)
#             droplets[j].X = droplets[k].X = 4/3*π * droplets[j].R^3
#             droplets[j].M = droplets[k].M = γ_tilde*droplets[j].M+droplets[k].M
#         elseif droplets[j].ξ - droplets[k].ξ < 0
#             print("nooooo")
#         end 
#     end

#     #this takes a lot of searching.. could we put the failure condition somewhere else/
#     # min_ξ = findmin(droplet -> droplet.ξ, droplets)
#     # max_ξ = findmax(droplet -> droplet.ξ, droplets)
#     # 
#     # if (min_ξ[1] <= 0 && max_ξ[1] > 1)
#     #     while (min_ξ[1] <= 0 && max_ξ[1] > 1)
#     #         # println("Superdroplet ", argmin(ξ), " has multiplicity of ", ξ[argmin(ξ)])
#     #         # Split superdroplet with highest multiplicity, half goes to argmin(ξ) and half stays
#     #         # Idea from Dziekan and Pawlowska (ACP, 2017) https://doi.org/10.5194/acp-17-13509-2017
#     #         min_ξ = findmin(droplet -> droplet.ξ, droplets)
#     #         max_ξ = findmax(droplet -> droplet.ξ, droplets)
            
#     #         argmin_i = min_ξ[2]
#     #         argmax_i = max_ξ[2]
#     #         droplets[argmin_i].ξ = floor(max_ξ[1]/2)
#     #         droplets[argmin_i].R = droplets[argmax_i].R
#     #         droplets[argmin_i].X = droplets[argmax_i].X
#     #         droplets[argmin_i].M = droplets[argmax_i].M
#     #         droplets[argmax_i].ξ = floor(max_ξ[1]/2)
#     #     end
#     # elseif (min_ξ[1] <= 0 && max_ξ[1] <= 1)
#     #     # Cannot split highest multiplicity superdroplet, have to remove superdroplet from system
#     #     println("Superdroplet ", min_ξ[2], " has multiplicity of ", min_ξ[1], ", removing from system")
#     #     deleteat!(droplets,droplets[min_ξ[2]])
#     #     Ns=Ns-1 #do I need to return this??
#     # end
# return droplets
# end



# function coalescence_timestep_small_alpha!(droplets,Ns,y,Δt,ΔV )#,M)
#     I =(sample(1:Ns, (y*2), replace = false))
#     L= [(I[l-1],I[l]) for l=2:2:length(I)]
#     scale = Ns*(Ns-1)/2/(y)


#     Threads.@threads for α=1:y
#     # for α in 1:y

#         jj=L[α][1] #struggled with choosing j>k so its hard coded
#         kk=L[α][2]

#         if droplets[jj].ξ<droplets[kk].ξ
#             k=jj
#             j=kk
#         else
#             j=jj
#             k=kk
#         end			

#         ϕ = rand()
#         pα = scale*calc_Ps(droplets[j],droplets[k],Δt,ΔV,kernel=kernel)


#         if ϕ < pα - floor(pα)
#             γ = floor(pα)+1
#         elseif ϕ >= pα - floor(pα)
#             γ = floor(pα)
#         end

#         if γ == 0
#             continue
#         end

#         γ_tilde = min(γ,floor(droplets[j].ξ/droplets[k].ξ))
#         # if γ_tilde == floor(ξ[j]/ξ[k])
#         #     println("limit exceeded")
#         # end

#         if droplets[j].ξ - γ_tilde*droplets[k].ξ > 0
#             droplets[j].ξ = droplets[j].ξ-γ_tilde*droplets[k].ξ
#             droplets[k].R = (γ_tilde*droplets[j].R^3+droplets[k].R^3)^(1/3)
#             droplets[k].X = 4/3*π * droplets[k].R^3
#             droplets[k].M = γ_tilde*droplets[j].M+droplets[k].M

#         elseif droplets[j].ξ - γ_tilde*droplets[k].ξ == 0
#             droplets[j].ξ = floor(droplets[k].ξ/2)
#             droplets[k].ξ = droplets[k].ξ - floor(droplets[k].ξ/2)
#             droplets[j].R = droplets[k].R = (γ_tilde*droplets[j].R^3+droplets[k].R^3)^(1/3)
#             droplets[j].X = droplets[k].X = 4/3*π * droplets[j].R^3
#             droplets[j].M = droplets[k].M = γ_tilde*droplets[j].M+droplets[k].M
#         elseif droplets[j].ξ - droplets[k].ξ < 0
#             print("nooooo")
#         end 
#     end
#         #this takes a lot of searching.. could we put the failure condition somewhere else/
#         # min_ξ = findmin(droplet -> droplet.ξ, droplets)
#         # max_ξ = findmax(droplet -> droplet.ξ, droplets)
        
#         # if (min_ξ[1] <= 0 && max_ξ[1] > 1)
#         #     while (min_ξ[1] <= 0 && max_ξ[1] > 1)
#         #         # println("Superdroplet ", argmin(ξ), " has multiplicity of ", ξ[argmin(ξ)])
#         #         # Split superdroplet with highest multiplicity, half goes to argmin(ξ) and half stays
#         #         # Idea from Dziekan and Pawlowska (ACP, 2017) https://doi.org/10.5194/acp-17-13509-2017
#         #         min_ξ = findmin(droplet -> droplet.ξ, droplets)
#         #         max_ξ = findmax(droplet -> droplet.ξ, droplets)
                
#         #         argmin_i = min_ξ[2]
#         #         argmax_i = max_ξ[2]
#         #         droplets[argmin_i].ξ = floor(max_ξ[1]/2)
#         #         droplets[argmin_i].R = droplets[argmax_i].R
#         #         droplets[argmin_i].X = droplets[argmax_i].X
#         #         droplets[argmin_i].M = droplets[argmax_i].M
#         #         droplets[argmax_i].ξ = floor(max_ξ[1]/2)
#         #     end
#         # elseif (min_ξ[1] <= 0 && max_ξ[1] <= 1)
#         #     # Cannot split highest multiplicity superdroplet, have to remove superdroplet from system
#         #     println("Superdroplet ", min_ξ[2], " has multiplicity of ", min_ξ[1], ", removing from system")
#         #     deleteat!(droplets,droplets[min_ξ[2]])
#         #     Ns=Ns-1 #do I need to return this??
#         # end
# return droplets
# end


