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

function coalescence_timestep!(run::Serial,scheme::KiD, ξFT::SVector,X::SVector,
    coag_data::coagulation_run,settings::coag_settings{FT}) where FT<:AbstractFloat
    
    first_healthy_idx = findfirst(iszero, ξFT)
    if first_healthy_idx === 1
        return (;SD_Vol = X,SD_Mult = ξFT)
    end
    Ns::Int = first_healthy_idx-1
    droplets = droplet_attributes{FT}(Int.(Vector(ξFT[1:Ns])), Vector(X[1:Ns]))
    
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

