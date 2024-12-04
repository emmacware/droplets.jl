#---------------------------------------------------------
#Visualization functions for droplet size distribution
#---------------------------------------------------------





#####

Base.@kwdef struct run_settings{FT<:AbstractFloat} #??
    num_bins::Int = Int(128)
    radius_bins_edges = 10 .^ range(log10(10*1e-6), log10(5e3*1e-6), length=num_bins+1) 
    smooth::Bool = true #bin smoothing
    smooth_scope::Int = Int(2)
    init_random_seed::Int = Int(30) 
    coag_threading = Parallel() # Serial(),use Julia NThreads for coalescence
    scheme = none() #Adaptive,Small_Alpha
    output_steps::Vector{FT} = [0,1200,2400,3600]
    init_method = init_logarithmic #init_ξ_const,init_uniform_sd
    # spectrum = exponential #only this for now
    binning_method = mass_density_lnr #number_density
    normalize_bins_dlnr = true
    # condensation::string = "none" #not yet
    #PySDM does steps instead of seconds+=dt... think about it
    #decide output here? bins,error,timing,...
    #should constants go in here? Or a constructor? 
end


#integrated error function, as in PySDM
function error_measure_old(y,ytrue,x)
    errors = ytrue .- y
    errors1 = errors[1:end-1] .+ errors[2:end]
    dx = diff(x)
    errors1 = errors1 .* dx ./2
    err = sum(abs.(errors1))
    return err
end

function error_measure(y,ytrue)
    return sqrt(mean((y.-ytrue).^2))
end

# Gaussian kernel density estimator function
# Useless honestly just bin it, smoothing function is there if you want it
function gnormal(R,X,ξ,Ns)
    σ₀ = 0.62
    σ = σ₀*Ns^(-1/5)
    W(Y) = 1/(2*π).^0.5/σ * exp.(-Y.^2 ./(2*σ^2))
    Y = [log(R[i]) .- log.(R) for i = 1:lastindex(R)]
    g = 1/ΔV * sum(ξ.*X.*ρl.*1e3.*W.(Y)) # 1e3 converts kg to g
    return g
end


function plot_dsd(bins,runsettings::run_settings{FT};color="black") where FT<:AbstractFloat   
    radius_bins_edges = runsettings.radius_bins_edges
    mids = 0.5*(radius_bins_edges[1:end-1] + radius_bins_edges[2:end])*1e6
    plot1 = plot!(mids,bins,lc=color,xaxis=:log,legend=false)
    return plot1
end

function ppmc_dsd(bins,runsettings::run_settings{FT};color="black") where FT<:AbstractFloat
    radius_bins_edges = runsettings.radius_bins_edges
    mids = 0.5*(radius_bins_edges[1:end-1] + radius_bins_edges[2:end])
    bins = replace(bins, 0 => 0.0001)
    plot1 = plot!(mids,bins,lc=color,xaxis=:log,yaxis=:log,legend=false,ylims = [1e9,2e14])
    return plot1
end



function number_density(droplets::droplets_allocations, t::FT,
    runsettings::run_settings{FT},coagsettings::coag_settings{FT}) where FT<:AbstractFloat

    i = sortperm(droplets.R)
    R::Vector{FT} = droplets.R[i]
    ξ::Vector{} = droplets.ξ[i]

    weights = ξ/coagsettings.ΔV

    numdens = binning_1d(R,weights,runsettings)

    if t != 0 && runsettings.smooth == true
        numdens = smoothbins!(numdens,runsettings)
    end

    return numdens
end


function mass_density_lnr(droplets::droplets_allocations, t::FT,
    runsettings::run_settings{FT},coagsettings::coag_settings{FT}) where FT<:AbstractFloat
    
    i = sortperm(droplets.R)
    Rsort::Vector{FT} = droplets.R[i]
    Xsort::Vector{FT} = droplets.X[i]
    ξsort::Vector{} = droplets.ξ[i]

    weights = ξsort .* Xsort
    weights *= 1e3 # convert to mass
    weights *= 1e3 # convert from kg to grams
    weights /= coagsettings.ΔV 

    numdens = binning_1d(Rsort,weights,runsettings)

    if t != 0 && runsettings.smooth == true
        numdens = smoothbins!(numdens,runsettings)
    end

    return numdens
end



function binning_1d(values::Vector{},weights::Vector{},runsettings::run_settings{FT}) where FT<:AbstractFloat
    bin_edges = runsettings.radius_bins_edges

    numdens::Vector{FT} = zeros(runsettings.num_bins)
    
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



function smoothbins!(numdens::Vector{FT},runsettings::run_settings{FT}) where FT<:AbstractFloat
    scope = runsettings.smooth_scope

    new_numdens = copy(numdens)
    for _ in 1:2

        for i in scope+1:length(numdens)-scope
            new_numdens[i] = mean(numdens[i-scope:i+scope])
        end
        scope = 1

        for i in scope+1:length(numdens)-scope
            numdens[i] = mean(new_numdens[i-scope:i+scope])
        end
    end
    return numdens
end