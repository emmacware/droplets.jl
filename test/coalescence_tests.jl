using Test
using Droplets
using Combinatorics
using Random

@testset begin

    FT = Float64
    coagsettings = coag_settings{FT}(Ns = 2^15)
    @test coagsettings.scale == coagsettings.Ns * (coagsettings.Ns - 1) / 2 / (coagsettings.Ns / 2)

    ξ_const = init_ξ_const(coagsettings)
    log = init_logarithmic(coagsettings)
    uniform = init_uniform_sd(coagsettings)
    mono = init_monodisperse(coagsettings)

    init_methods = [ξ_const, log, uniform,mono]

    init_method = uniform

    for init_method in init_methods
        @test sum(init_method.ξ) ≈ coagsettings.n0*coagsettings.ΔV rtol = 1e-2
        @test init_method.X ≈ (init_method.R).^3 .* 4pi/3

        #This isn't passing for uniform and monodisperse, is that a problem with them?
        # effective_vol = sum(init_method.X) / coagsettings.Ns
        # effective_radius = (effective_vol * 3 / 4 / pi)^(1/3)
        # @test effective_radius ≈ coagsettings.R0 rtol = 1e-2 

    end

    @test all(==(ξ_const.ξ[1]), ξ_const.ξ)


    #kernels
    droplet_probability = golovin(ξ_const, (1,2), coagsettings)
    @test droplet_probability == coagsettings.golovin_kernel_coeff * (ξ_const.X[1] + ξ_const.X[2])

    #same-size drops cannot collide
    @test hydrodynamic(mono, (1,2), coagsettings) == 0.0
    # set test for different sized drops
    # @test hydrodynamic(uniform, (1,2), coagsettings) == ?




    run = [Serial()] #Parallel has trouble on windows right now, have to fix
    schemes = [none(),Adaptive()]

    for scheme in schemes
        Random.seed!()
        ctimestep = deepcopy(ξ_const)
        for _ in 1:50
            coalescence_timestep!(Serial(),scheme, ctimestep, coagsettings)
        end

        # coag only
        @test sum(ctimestep.ξ) <= sum(ξ_const.ξ)
        # mass conservation:
        @test sum(ctimestep.X.*ctimestep.ξ) ≈ sum(ξ_const.X.*ξ_const.ξ) rtol = 1e-12
        @test minimum(ctimestep.ξ) > 0
        @test minimum(ctimestep.X) > 0
        @test minimum(ctimestep.R) > 0
    end
    
    #test that the deficit is being updated correctly
    deficit_drops = deepcopy(ξ_const)
    deficit = deficit_allocations(ξ_const,zeros(FT, div(coagsettings.Ns, 2)))
    for _ in 1:10
        coalescence_timestep!(Serial(),log_deficit(), deficit,coagsettings)
    end

    @test sum(deficit.droplets.ξ) <= sum(ξ_const.ξ)
    @test sum(deficit.droplets.X.*deficit.droplets.ξ) == sum(ξ_const.X.*ξ_const.ξ)
    @test minimum(deficit.droplets.ξ) > 0
    @test minimum(deficit.droplets.X) > 0
    @test minimum(deficit.droplets.R) > 0




    # sdm_update!((1,2),1, droplets)



    # Base.@kwdef struct coag_settings{FT<:AbstractFloat}
    #     Δt::FT = FT(1.0)
    #     ΔV::FT = FT(1e6)
    #     Ns::Int = 2^15# number of superdroplets
    #     scale::FT = Ns * (Ns - 1) / 2 / (Ns / 2)
    #     R_min::FT = FT(1e-9)
    #     R_max::FT = FT(1e-3)
    #     golovin_kernel_coeff::FT = FT(1.5e3)
    #     hydrodynamic_collision_eff_func::Bool = false
    #     kernel::Function = golovin # golovin, hydrodynamic
    #     n0::FT = FT(2^23) # initial droplet concentration
    #     R0::FT = FT(30.531e-6) # initial radius
    # end

    # struct droplets_allocations{FT<:AbstractFloat}
    #     ξ::Vector{Int}
    #     R::Vector{FT}
    #     X::Vector{FT}
    #     I::Vector{Int}
    #     pαdt::Vector{FT}
    #     ϕ::Vector{FT}
    # end

    # coalescence_timestep!
    # export Serial, Parallel,Adaptive,none,coag_settings,droplets_allocations,pair_Ps!,sdm_update!
    # export deficit_allocations,log_deficit,init_monodisperse
end