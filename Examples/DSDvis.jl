#---------------------------------------------------------
#Visualization functions for droplet size distribution
#---------------------------------------------------------





#####

Base.@kwdef struct run_settings{FT<:AbstractFloat} #??
    num_bins::Int = Int(128)
    radius_bins_edges = 10 .^ range(log10(10*1e-6), log10(5e3*1e-6), length=num_bins) 
    smooth::Bool = true #bin smoothing
    smooth_scope::Int = Int(2)
    init_random_seed::Int = Int(30) 
    coag_threading = Parallel() # Serial(),use Julia NThreads for coalescence
    scheme = none() #Adaptive,Small_Alpha
    output_steps::Vector{FT} = [0,1200,2400,3600]
    init_method = init_logarithmic #init_ξ_const,init_uniform_sd
    # spectrum = exponential #only this for now
    binning_method = mass_density_lnr #number_density
    # condensation::string = "none" #not yet
    #PySDM does steps instead of seconds+=dt... think about it
    #decide output here? bins,error,timing,...
    #should constants go in here? Or a constructor? 
end



function coag_runtime(randseed::Int,ξ::Vector{Int},R::Vector{FT},X::Vector{FT},
    coag_settings::coag_settings{FT},run_settings::run_settings{FT}) where FT<:AbstractFloat
    
    Random.seed!(randseed)
    println("Running simulation...")

    coal_func_time::FT = 0.0
    bins::Matrix{FT} = zeros(FT, run_settings.num_bins - 1, 4)

    ϕ = Vector{FT}(undef, div(coag_settings.Ns, 2))
    I::Vector{Int} = shuffle(1:coag_settings.Ns)
    simtime::FT = @elapsed begin
        for i  in  1:length(run_settings.output_steps)
            # if i,seconds in enumerate(run_settings.output_steps)
            
            if i ==1
                bins[:,i] = run_settings.binning_method(X, ξ,run_settings.output_steps[i],run_settings)
                println("Time: ", run_settings.output_steps[i], " seconds")
                continue
            end

            timestepper = (run_settings.output_steps[i]-run_settings.output_steps[i-1])/coag_settings.Δt
            ctime::FT = @elapsed begin
                for _ in 1:timestepper
                    coalescence_timestep!(run_settings.coag_threading, run_settings.scheme, ξ, R, X,I,ϕ,coag_settings)
                end
            end
            coal_func_time += ctime
            bins[:,i] = mass_density_lnr(X, ξ,run_settings.output_steps[i],run_settings)
            println("Time: ", run_settings.output_steps[i], " seconds")
        end
    end
    println("simtime =", simtime)
    println("coal_func_time =", coal_func_time)

    return bins, coal_func_time
end


#binning_dsd function for superdroplets, superdroplet struct input
#becuase of what it was used for, time is an input so that it would not smooth at t=0 
function binning_dsd(droplets,bins,t;smooth = true,scope_init = 1)
    droplets = sort(droplets, by = droplet -> droplet.R)
    bin_edges = bins
    # mids = 0.5*(bin_edges[1:end-1] + bin_edges[2:end])*1e6 #um
    mids = bin_edges[1:end-1]*1e6 #um
    water = zeros(length(mids))
    droplet_idx = 1
    for j in 1:length(mids)

        lo = 4*pi .*((bin_edges[j]) .^ 3) ./3 
        hi = 4*pi .*((bin_edges[j+1]) .^ 3) ./3
        leng = log(bin_edges[j+1])-log(bin_edges[j])
        for i in droplet_idx:length(droplets)

            if droplets[i].X < hi
                water[j]=water[j]+droplets[i].X*droplets[i].ξ/leng
                droplet_idx += 1
            else
                break
            end

        end

    end
    scope = scope_init
    if t != 0 && smooth == true
        new_water = copy(water)
        for j in 1:2
            # scope = scope_init
            for i in scope+1:length(water)-scope
                new_water[i] = mean(water[i-scope:i+scope])
            end
            scope = 1
            for i in scope+1:length(water)-scope
                water[i] = mean(new_water[i-scope:i+scope])
            end
        end
    end
    return mids, water
end


# function mass_density_lnr(Xunsorted,ξunsorted,t,settings)
#     bin_edges = settings.radius_bins_edges
#     i = sortperm(Xunsorted)
#     X = Xunsorted[i]
#     ξ = ξunsorted[i]

#     # mids = 0.5*(bin_edges[1:end-1] + bin_edges[2:end])*1e6 #um
#     mids = bin_edges[1:end-1]*1e6 #um

#     water = zeros(length(mids))
#     droplet_idx = 1
#     for j in 1:length(mids)

#         # lo = 4*pi .*((bin_edges[j]) .^ 3) ./3 
#         hi = 4*pi .*((bin_edges[j+1]) .^ 3) ./3
#         leng = log(bin_edges[j+1])-log(bin_edges[j])
#         for i in droplet_idx:length(X)

#             if X[i] < hi
#                 water[j]=water[j]+X[i]*ξ[i]/leng
#                 droplet_idx += 1
#             else
#                 break
#             end

#         end

#     end

#     scope = settings.smooth_scope

#     if t != 0 && settings.smooth == true
#         new_water = copy(water)
#         for _ in 1:2
#             # scope = scope_init
#             # println(scope)
#             for i in scope+1:length(water)-scope
#                 new_water[i] = mean(water[i-scope:i+scope])
#             end
#             scope = 1
#             # println(scope)
#             for i in scope+1:length(water)-scope
#                 water[i] = mean(new_water[i-scope:i+scope])
#             end
#         end


#     end

#     return mids, water

# end

function mass_density_lnr(Xunsorted::Vector{FT}, ξunsorted::Vector{Int}, 
    t::FT,settings::run_settings{FT}) where FT<:AbstractFloat
    
    bin_edges = settings.radius_bins_edges
    i = sortperm(Xunsorted)
    X::Vector{FT} = Xunsorted[i]
    ξ::Vector{} = ξunsorted[i]

    # mids = 0.5*(bin_edges[1:end-1] + bin_edges[2:end])*1e6 #um
    mids = bin_edges[1:end-1] * 1e6 # um

    water = zeros(length(mids))
    droplet_idx = 1
    for j in 1:length(mids)
        # lo = 4 * π * ((bin_edges[j]) ^ 3) / 3 
        hi = 4 * π * ((bin_edges[j+1]) ^ 3) / 3
        leng = log(bin_edges[j+1]) - log(bin_edges[j])
        for i in droplet_idx:length(X)
            if X[i] < hi
                water[j] = water[j] + X[i] * ξ[i] / leng
                droplet_idx += 1
            else
                break
            end
        end
    end

    scope = settings.smooth_scope

    if t != 0 && settings.smooth == true
        new_water = copy(water)
        for _ in 1:2
            # scope = scope_init
            # println(scope)
            for i in scope+1:length(water)-scope
                new_water[i] = mean(water[i-scope:i+scope])
            end
            scope = 1
            # println(scope)
            for i in scope+1:length(water)-scope
                water[i] = mean(new_water[i-scope:i+scope])
            end
        end
    end

    return water
end



# Smoothing function for bins, moving average 
function smoothbins!(bins,scope_init=2)
    scope = scope_init
    new_bins = copy(bins)
    for j in 1:2
        # scope = scope_init
        for i in scope+1:length(bins)-scope
            new_bins[i] = mean(bins[i-scope:i+scope])
        end
        scope = 1
        for i in scope+1:length(bins)-scope
            bins[i] = mean(new_bins[i-scope:i+scope])
        end
    end
    return bins
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


function plot_dsd(bins;color="black",radius_bins_edges = 10 .^ range(log10(10*1e-6), log10(5e3*1e-6), length=128))

     
    mids = 0.5*(radius_bins_edges[1:end-1] + radius_bins_edges[2:end])*1e6
    plot1 = plot!(mids,bins,lc=color,xaxis=:log,legend=false)

    return plot1
end

function ppmc_dsd(bins;color="black",radius_bins_edges = 10 .^ range(log10(1e-7), log10(1e-3), length=50))

     
    mids = 0.5*(radius_bins_edges[1:end-1] + radius_bins_edges[2:end])
    bins = replace(bins, 0 => 0.0001)
    plot1 = plot!(mids,bins,lc=color,xaxis=:log,yaxis=:log,legend=false,ylims = [1e9,2e14])



    return plot1
end





#binning_dsd function for superdroplets, vector input
#becuase of what it was used for, time is an input so that it would not smooth at t=0 
# function number_density_dsd(Xunsorted,ξunsorted,bins,t;smooth = true,scope_init = 1)
function number_density_dsd(Xunsorted::Vector{FT}, ξunsorted::Vector{Int}, 
    t::FT,settings::run_settings{FT}) where FT<:AbstractFloat

    bin_edges = settings.radius_bins_edges
    i = sortperm(Xunsorted)
    X = Xunsorted[i]
    ξ = ξunsorted[i]

    # mids = 0.5*(bin_edges[1:end-1] + bin_edges[2:end])*1e6 #um
    mids = bin_edges[1:end-1] #m

    numdens = zeros(length(mids))
    droplet_idx = 1
    for j in 1:length(mids)

        lo = 4*pi .*((bin_edges[j]) .^ 3) ./3 
        hi = 4*pi .*((bin_edges[j+1]) .^ 3) ./3
        leng = (bin_edges[j+1])-(bin_edges[j])
        for i in droplet_idx:length(X)

            if X[i] < hi
                numdens[j]=numdens[j]+ξ[i]/leng
                droplet_idx += 1
            else
                break
            end

        end

    end

    scope = scope_init

    if t != 0 && smooth == true
        new_numdens = copy(numdens)
        for j in 1:2
            # scope = scope_init
            # println(scope)
            for i in scope+1:length(numdens)-scope
                new_numdens[i] = mean(numdens[i-scope:i+scope])
            end
            scope = 1
            # println(scope)
            for i in scope+1:length(numdens)-scope
                numdens[i] = mean(new_numdens[i-scope:i+scope])
            end
        end


    end

    return mids, numdens

end


# #binning_dsd function for superdroplets, vector input
# #becuase of what it was used for, time is an input so that it would not smooth at t=0 
# function number_density_dsd(Xunsorted,ξunsorted,bins,t;smooth = true,scope_init = 1)
#     bin_edges = bins
#     i = sortperm(Xunsorted)
#     X = Xunsorted[i]
#     ξ = ξunsorted[i]

#     # mids = 0.5*(bin_edges[1:end-1] + bin_edges[2:end])*1e6 #um
#     mids = bin_edges[1:end-1] #m

#     numdens = zeros(length(mids))
#     droplet_idx = 1
#     for j in 1:length(mids)

#         lo = 4*pi .*((bin_edges[j]) .^ 3) ./3 
#         hi = 4*pi .*((bin_edges[j+1]) .^ 3) ./3
#         leng = (bin_edges[j+1])-(bin_edges[j])
#         for i in droplet_idx:length(X)

#             if X[i] < hi
#                 numdens[j]=numdens[j]+ξ[i]/leng
#                 droplet_idx += 1
#             else
#                 break
#             end

#         end

#     end

#     scope = scope_init

#     if t != 0 && smooth == true
#         new_numdens = copy(numdens)
#         for j in 1:2
#             # scope = scope_init
#             # println(scope)
#             for i in scope+1:length(numdens)-scope
#                 new_numdens[i] = mean(numdens[i-scope:i+scope])
#             end
#             scope = 1
#             # println(scope)
#             for i in scope+1:length(numdens)-scope
#                 numdens[i] = mean(new_numdens[i-scope:i+scope])
#             end
#         end


#     end

#     return mids, numdens

# end

#binning_dsd function for superdroplets, vector input
#becuase of what it was used for, time is an input so that it would not smooth at t=0 
function number_density(Xunsorted,ξunsorted,t,settings)
    bin_edges = settings.radius_bins_edges
    i = sortperm(Xunsorted)
    X = Xunsorted[i]
    ξ = ξunsorted[i]
    # mids = 0.5*(bin_edges[1:end-1] + bin_edges[2:end])*1e6 #um
    mids = bin_edges[1:end-1] #m
    numdens = zeros(length(mids))
    droplet_idx = 1
    for j in 1:length(mids)
        # lo = 4*pi .*((bin_edges[j]) .^ 3) ./3 
        hi = 4*pi .*((bin_edges[j+1]) .^ 3) ./3
        leng = (bin_edges[j+1])-(bin_edges[j])
        for i in droplet_idx:length(X)
            if X[i] < hi
                numdens[j]=numdens[j]+ξ[i]/leng
                droplet_idx += 1
            else
                break
            end
        end
    end
    scope = settings.smooth_scope
    if t != 0 && settings.smooth == true
        new_numdens = copy(numdens)
        for _ in 1:2
            # scope = scope_init
            # println(scope)
            for i in scope+1:length(numdens)-scope
                new_numdens[i] = mean(numdens[i-scope:i+scope])
            end
            scope = 1
            # println(scope)
            for i in scope+1:length(numdens)-scope
                numdens[i] = mean(new_numdens[i-scope:i+scope])
            end
        end
    end
    return mids, numdens
end
