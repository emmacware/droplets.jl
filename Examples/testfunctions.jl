using Plots

# function build_box_run(coag_settings::coag_settings{FT}, run_settings::run_settings{FT}) where FT<:AbstractFloat
#     ξ, R, X = run_settings.init_method(coag_settings)
#     Ns::Int = coag_settings.Ns
#     I::Vector{Int} = shuffle(1:Ns)
    
#     function coag_runtime(randseed::Int,ξ::Vector{Int},R::Vector{FT},X::Vector{FT},
#         I::Vector{Int},coag_settings::coag_settings{FT},run_settings::run_settings{FT}) where FT<:AbstractFloat
#         Random.seed!(randseed)
        
#         println("Running simulation...")
#         seconds::FT = 0.0

#         coal_func_time::FT = 0.0
#         bins::Matrix{FT} = zeros(FT, run_settings.num_bins - 1, 0)

#         ϕ = Vector{FT}(undef, div(coag_settings.Ns, 2))
#         simtime::FT = @elapsed begin
#             while seconds <= run_settings.output_steps[end]
#                 if seconds in run_settings.output_steps
#                     println("Time: ", seconds, " seconds")
#                     xx, yy = run_settings.binning_method(X, ξ,seconds, run_settings)
#                     bins = hcat(bins, yy)
#                 end

#                 ctime::FT = @elapsed begin
#                     coalescence_timestep!(run_settings.coag_threading, run_settings.scheme, ξ, R, X, I,ϕ,coag_settings)
#                 end
#                 coal_func_time += ctime
#                 seconds += coag_settings.Δt
#                 # GC.gc()
#             end
#         end
#         println("simtime =", simtime)
#         println("coal_func_time =", coal_func_time)

#         return bins, coal_func_time
#     end
#     return coag_runtime
# end

function coag_runtime(randseed::Int,droplets::droplets_allocations,
    coag_settings::coag_settings{FT},run_settings::run_settings{FT}) where FT<:AbstractFloat
    
    Random.seed!(randseed)
    println("Running simulation...")

    coal_func_time::FT = 0.0
    bins::Matrix{FT} = zeros(FT, run_settings.num_bins, length(run_settings.output_steps))
    threading,scheme = run_settings.coag_threading, run_settings.scheme
    simtime::FT = @CPUelapsed begin
        for i  in  1:length(run_settings.output_steps)
            # if i,seconds in enumerate(run_settings.output_steps)
            
            if i !=1
                timestepper = (run_settings.output_steps[i]-run_settings.output_steps[i-1])/coag_settings.Δt
                ctime::FT = @CPUelapsed begin
                    for _ in 1:timestepper
                        coalescence_timestep!(threading,scheme,droplets,coag_settings)
                    end
                end
                coal_func_time += ctime
            end
            bins[:,i] = run_settings.binning_method(droplets,run_settings.output_steps[i],run_settings,coag_settings)
            println("Time: ", run_settings.output_steps[i], " seconds")
        end
    end
    println("simtime =", simtime)
    println("coal_func_time =", coal_func_time)

    return bins, coal_func_time
end


function coag_runtime_log_deficit(randseed::Int,droplets::deficit_allocations,
    coag_settings::coag_settings{FT},run_settings::run_settings{FT}) where FT<:AbstractFloat
    deficit_droplets = zeros(FT, Int(run_settings.output_steps[end]/coag_settings.Δt))
    deficit_pairs = zeros(FT, Int(run_settings.output_steps[end]/coag_settings.Δt))
    
    Random.seed!(randseed)
    println("Running simulation...")

    coal_func_time::FT = 0.0
    bins::Matrix{FT} = zeros(FT, run_settings.num_bins, 4)
    threading,scheme = run_settings.coag_threading, run_settings.scheme

    simtime::FT = @CPUelapsed begin
        for i  in  1:length(run_settings.output_steps)
            if i !=1
                timestepper = (run_settings.output_steps[i]-run_settings.output_steps[i-1])/coag_settings.Δt
                ctime::FT = @CPUelapsed begin
                    for t in 1:timestepper
                        coalescence_timestep!(threading,log_deficit(),droplets,coag_settings)
                        deficit_droplets[Int((run_settings.output_steps[i-1]/coag_settings.Δt)+t)] = sum(droplets.deficit)
                        deficit_pairs[Int((run_settings.output_steps[i-1]/coag_settings.Δt)+t)] = count(x -> x != 0, droplets.deficit)
                    end
                end
                coal_func_time += ctime
            end
            bins[:,i] = run_settings.binning_method(droplets.droplets,run_settings.output_steps[i],run_settings,coag_settings)
            println("Time: ", run_settings.output_steps[i], " seconds")
        end
    end
    println("simtime =", simtime)
    println("coal_func_time =", coal_func_time)

    return bins, coal_func_time,deficit_droplets,deficit_pairs
end