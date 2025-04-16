using Distributions
using Combinatorics
using Random
using Interpolations
using StaticArrays

export coalescence_timestep!, static_droplet_attributes, coagulation_run, coag_settings, run_settings
export KiD

struct KiD end #<: scheme_type end

struct static_droplet_attributes{FT<:AbstractFloat,NSD}
    ξ::SVector{NSD,FT}
    X::SVector{NSD,FT}
end


# function coalescence_timestep!(run::Serial,scheme::KiD, droplets::static_droplet_attributes,
#     coag_data::coagulation_run,settings::coag_settings{FT}) where FT<:AbstractFloat
    
#     Ns::Int = settings.Ns
#     HealthyIdx = findfirst(iszero, droplets.ξ)
#     HealthyIdx = HealthyIdx === nothing ? Ns : HealthyIdx - 1
#     #or findfirst(istrue, healthy)

    
#     idx = shuffle(1:HealthyIdx)
#     #coagdata.I[1:HealthyIdx] = shuffle(1:HealthyIdx)
#     L = [(idx[l-1], idx[l]) for l in 2:2:HealthyIdx]
#     # L = [(coag_data.I[l-1], coag_data.I[l]) for l in 2:2:Ns]
#     # L = collect(partition(idx, 2)) IterTools...

#     compute_pαdt!(L, droplets,coag_data,settings.kernel, settings)

#     rand!(coag_data.ϕ)

#     test_pairs!(run,Ns,L,droplets,coag_data)

# end 

# function coalescence_timestep!(run::Serial,scheme::KiD, ξFT::SVector,X::SVector,
#     coag_data::coagulation_run,settings::coag_settings{FT}) where FT<:AbstractFloat

#     R = (volume_to_radius.(X))
#     ξ = SVector{length(ξFT), Int}(map(Int, ξFT))

#     droplets = static_droplet_attributes(ξ, R, X)
    
#     Ns::Int = settings.Ns
#     HealthyIdx = findfirst(iszero, droplets.ξ)
#     HealthyIdx = HealthyIdx === nothing ? Ns : HealthyIdx - 1

#     if HealthyIdx == 0
#         return (;SD_Vol = X,SD_Mult = ξFT)
#     end
    
#     idx = shuffle(1:HealthyIdx)
#     #coagdata.I[1:HealthyIdx] = shuffle(1:HealthyIdx)
#     L = [(idx[l-1], idx[l]) for l in 2:2:HealthyIdx]
#     # L = [(coag_data.I[l-1], coag_data.I[l]) for l in 2:2:Ns]
#     # L = collect(partition(idx, 2)) IterTools...

#     compute_pαdt!(L, droplets,coag_data,settings.kernel, settings)

#     rand!(coag_data.ϕ)

#     ξint, X = test_pairs!(run,Ns,L,droplets,coag_data)

#     ξFT = SVector{length(ξFT), FT}(map(FT, ξint))
#     # X = droplets.X
#     return (;SD_Vol = X,SD_Mult = ξFT)
# end 

function coalescence_timestep!(run::Serial,scheme::KiD, ξFT::SVector,X::SVector,
    coag_data::coagulation_run,settings::coag_settings{FT}) where FT<:AbstractFloat
    
    first_healthy_idx = findfirst(iszero, ξFT)
    if first_healthy_idx === 1
        return (;SD_Vol = X,SD_Mult = ξFT)
    end
    Ns::Int = first_healthy_idx-1
    droplets = droplet_attributes{FT}(Int.(Vector(ξFT[1:Ns])), Vector((volume_to_radius.(X[1:Ns]))), Vector(X[1:Ns]))
    
    I = shuffle!(1:Ns)
    L = [(I[l-1], I[l]) for l in 2:2:Ns]

    compute_pαdt!(L, droplets,coag_data,settings.kernel, settings)

    rand!(coag_data.ϕ)

    test_pairs!(run,Ns,L,droplets,coag_data)
    ξFT = SVector{length(ξFT), FT}(FT.(droplets.ξ))
    X = SVector{length(ξFT), FT}(droplets.X)
    return (;SD_Vol = X,SD_Mult = ξFT)
end 


function size_thresh_separate_droplets(SD_Vol,SD_Mult,ρd)
    FT = eltype(SD_Vol)
    # Define the threshold for separating droplets
    threshold = 40*1e-6 # meters
    threshold_volume = (4/3) * π * (threshold^3) # m^3

    # Separate the droplets based on the threshold
    N_liq = sum(SD_Mult[SD_Vol .< threshold])
    N_rai = sum(SD_Mult[SD_Vol .>= threshold])

    q_liq = 1000*sum(SD_Vol[SD_Vol .< threshold].* SD_Mult[SD_Vol .< threshold])/ρd
    q_rai = 1000*sum(SD_Vol[SD_Vol .>= threshold].* SD_Mult[SD_Vol .>= threshold])/ρd

    return (N_liq=FT(N_liq), N_rai=FT(N_rai), q_liq=FT(q_liq), q_rai=FT(q_rai))
end

# function coalescence_timestep!(run::Serial,scheme::KiD, droplets::static_droplet_attributes{FT,NSD},
#     coag_data::coagulation_run,settings::coag_settings{FT}) where {FT, NSD}
    
#     Ns::Int = settings.Ns
#     droplets = static_to_droplets_tmp(droplets)
    
#     shuffle!(coag_data.I)
#     L = [(coag_data.I[l-1], coag_data.I[l]) for l in 2:2:Ns]

#     compute_pαdt!(L, droplets,coag_data,settings.kernel, settings)

#     rand!(coag_data.ϕ)

#     test_pairs!(run,Ns,L,droplets,coag_data)
#     ξFT = SVector{length(ξFT), FT}(FT.(droplets.ξ))
#     X = SVector{length(ξFT), FT}(droplets.X)

#     return static_droplet_attributes{FT,NSD}(ξFT, X)
# end 



@inline function compute_pαdt!(L::Vector{Tuple{Int,Int}}, droplets::static_droplet_attributes,coag_data::coagulation_run,kernel::Function,coagsettings::coag_settings{FT}) where FT<:AbstractFloat
    map(i -> pair_Ps!(i, L[i], droplets,coag_data,kernel, coagsettings), eachindex(L))#coag_data.pαdt))
    coag_data.pαdt .*=  coagsettings.scale * coagsettings.Δt / coagsettings.ΔV
end

@inline function pair_Ps!(α::Int, (j,k)::Tuple{Int,Int}, droplets::static_droplet_attributes,coag_data::coagulation_run,kernel::Function,coagsettings::coag_settings{FT}) where FT<:AbstractFloat
    coag_data.pαdt[α] = max(droplets.ξ[j], droplets.ξ[k]) * kernel(droplets,(j,k), coagsettings)
end

@inline function golovin(droplets::static_droplet_attributes, (j,k)::Tuple{Int,Int}, settings::coag_settings{FT})::FT where FT<:AbstractFloat
    return settings.golovin_kernel_coeff *(droplets.X[j] + droplets.X[k])# Xsum
end

function test_pairs!(backend::Serial,Ns::Int,L::Vector{Tuple{Int,Int}},droplets::static_droplet_attributes{FT},coag_data::coagulation_run) where FT<:AbstractFloat
    
    coag_data.lowest_zero[] = false
    for α::Int in 1:length(L)#div(Ns, 2)
            
        if coag_data.ϕ[α] >= coag_data.pαdt[α]
            continue
        end
        droplets = sdm_update!(L[α], α, droplets,coag_data)
    end
    # if coag_data.lowest_zero[] == true
    #     split_highest_multiplicity!(droplets)
    # end
    ξint, X = droplets.ξ, droplets.X
end

function sdm_update!(pair::Tuple{Int, Int}, α::Int, droplets::static_droplet_attributes{FT}, coag_data::coagulation_run) where FT<:AbstractFloat
    if droplets.ξ[pair[1]] < droplets.ξ[pair[2]]
        k = pair[1]
        j = pair[2]
    else
        k = pair[2]
        j = pair[1]
    end

    ξj, ξk = droplets.ξ[j], droplets.ξ[k]
    pα = coag_data.pαdt[α]

    pα_floor::FT = @fastmath floor(pα)
    γ::FT = if coag_data.ϕ[α] < pα - pα_floor
        pα_floor + 1
    else
        pα_floor
    end

    floor_ξj_div_ξk = floor(ξj / ξk)
    if γ >= floor_ξj_div_ξk
        coag_data.deficit[] += (γ - floor_ξj_div_ξk) * ξk
        γ = floor_ξj_div_ξk
    end

    ξ_j_minus_γ_tilde_ξ_k = ξj - γ * ξk
    volume = γ * droplets.X[j] + droplets.X[k]


    tmpξ = droplets.ξ
    tmpR = droplets.R
    tmpX = droplets.X

    if ξ_j_minus_γ_tilde_ξ_k > 0
        tmpξ = setindex(tmpξ, ξ_j_minus_γ_tilde_ξ_k, j)
        tmpR = setindex(tmpR, volume_to_radius(volume), k)
        tmpX = setindex(tmpX, volume, k)
    elseif ξ_j_minus_γ_tilde_ξ_k == 0
        half_ξ_k = floor(ξk / 2)
        tmpξ = setindex(tmpξ, half_ξ_k, j)
        tmpξ = setindex(tmpξ, tmpξ[k] - half_ξ_k, k)
        tmpR = setindex(tmpR, volume_to_radius(volume), j)
        tmpR = setindex(tmpR, volume_to_radius(volume), k)
        tmpX = setindex(tmpX, volume, j)
        tmpX = setindex(tmpX, volume, k)
        if half_ξ_k == 0
            coag_data.lowest_zero[] = true
        end
    else
        println("nooooo")
    end
    droplets = static_droplet_attributes(tmpξ, tmpR, tmpX)
    return droplets
end

function binning_func(droplets::static_droplet_attributes, t::FT,
    runsettings::run_settings{FT},coagsettings::coag_settings{FT}) where FT<:AbstractFloat

    weights = runsettings.binning_method(droplets,coagsettings)
    numdens = binning_1d(droplets.R,weights,runsettings)

    if t != 0 && runsettings.smooth == true
        numdens = smoothbins!(numdens,runsettings)
    end

    return numdens
end

function mass_density_lnr(droplets::static_droplet_attributes,coagsettings::coag_settings{FT}) where FT<:AbstractFloat

    tot_vol = droplets.ξ .* droplets.X
    weights_kilograms = tot_vol*constants.ρl # convert to mass
    kg_per_vol = weights_kilograms/coagsettings.ΔV 

    return kg_per_vol
end


function split_highest_multiplicity(droplets::static_droplet_attributes)
    ξ = droplets.ξ
    R = droplets.R
    X = droplets.X

    if maximum(ξ) > 1
        while (minimum(ξ) <= 0 && maximum(ξ) > 1)
            argmin_i = argmin(ξ)
            argmax_i = argmax(ξ)

            # Compute new values for ξ, R, and X
            new_ξ = setindex(ξ, floor(ξ[argmax_i] / 2), argmin_i)
            new_ξ = setindex(new_ξ, ξ[argmax_i] - floor(ξ[argmax_i] / 2), argmax_i)

            new_R = setindex(R, R[argmax_i], argmin_i)
            new_X = setindex(X, X[argmax_i], argmin_i)

            # Update the droplets with the new values
            ξ, R, X = new_ξ, new_R, new_X
        end
    elseif maximum(ξ) <= 1
        println("Highest superdroplet cannot be split")
        if maximum(ξ) < 1
            error("Highest and Lowest Superdroplet have ξ == 0")
        end
    end

    # Return a new instance of static_droplet_attributes with updated values
    droplets = static_droplet_attributes{FT}(ξ, R, X)
end

function binning_1d(values_unsorted::SVector{},weights_unsorted::SVector{},runsettings::run_settings{FT}) where FT<:AbstractFloat
    bin_edges = runsettings.radius_bins_edges

    numdens::Vector{FT} = zeros(runsettings.num_bins)
    idx = sortperm(values_unsorted)
    values = values_unsorted[idx]
    weights = weights_unsorted[idx]

    droplet_idx = 1
    for j in 1:runsettings.num_bins
        bin_edge_high = (bin_edges[j+1])
        for i in droplet_idx:length(values)
            if values[i] < bin_edge_high
                numdens[j]=numdens[j]+weights[i]
                droplet_idx += 1
            else
                break
            end
        end
    end

    if runsettings.normalize_bins_dlnr == true
        numdens = numdens ./ diff(log.(bin_edges))
    end
    return numdens
end