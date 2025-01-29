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

    init_methods = [ξ_const, log,mono]

    init_method = uniform

    for init_method in init_methods
        @test sum(init_method.ξ) ≈ coagsettings.n0*coagsettings.ΔV rtol = 1e-2
        @test init_method.X ≈ radius_to_volume.(init_method.R)

        #This isn't passing for uniform, I feel like thats a problem
        effective_vol = sum(init_method.X) / coagsettings.Ns
        effective_radius = volume_to_radius(effective_vol)
        @test effective_radius ≈ coagsettings.R0 rtol = 1e-2 
    
    end

    @test all(==(ξ_const.ξ[1]), ξ_const.ξ)


    #kernels
    droplet_probability = golovin(ξ_const, (1,2), coagsettings)
    @test droplet_probability == coagsettings.golovin_kernel_coeff * (ξ_const.X[1] + ξ_const.X[2])

    big_golovin_coefficient_settings = coag_settings{FT}(Ns = 2^15, golovin_kernel_coeff = 1700)
    small_golovin_coefficient_settings = coag_settings{FT}(Ns = 2^15, golovin_kernel_coeff = 1300)
    @test golovin(ξ_const, (1,2), big_golovin_coefficient_settings) > golovin(ξ_const, (1,2), small_golovin_coefficient_settings)

    #same-size drops cannot collide
    @test hydrodynamic(mono, (1,2), coagsettings) == 0.0
    # set test for different sized drops
    # @test hydrodynamic(uniform, (1,2), coagsettings) == ?



    run = [Serial()] #Parallel has trouble on windows right now, have to fix
    schemes = [none(),Adaptive()]
    coag_data = coagulation_run{FT}(coagsettings.Ns)

    for scheme in schemes
        Random.seed!()
        ctimestep = deepcopy(ξ_const)
        for _ in 1:50
            coalescence_timestep!(Serial(),scheme, ctimestep,coag_data, coagsettings)
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
    deficit = deficit_allocations{FT}(deficit_drops)
    for _ in 1:10
        coalescence_timestep!(Serial(),log_deficit(), deficit,coag_data,coagsettings)
    end

    @test sum(deficit.droplets.ξ) <= sum(ξ_const.ξ)
    @test sum(deficit.droplets.X.*deficit.droplets.ξ) ≈ sum(ξ_const.X.*ξ_const.ξ) rtol = 1e-12
    @test minimum(deficit.droplets.ξ) > 0
    @test minimum(deficit.droplets.X) > 0
    @test minimum(deficit.droplets.R) > 0

    #test that deficit for Adaptive is zero
    coag_data = coagulation_run{FT}(coagsettings.Ns)    
    for _ in 1:10
        coalescence_timestep!(Serial(),Adaptive(), deficit_drops,coag_data, coagsettings)
    end
    @test coag_data.deficit[] == 0

    # test sdm update for extreme cases: 
    # multiple collisions, split highest multiplicity, lowest_zero
    two_drops_coagsettings = coag_settings{FT}(Ns = 2)
    R = [FT(1.5e-6), FT(1.5e-6)]
    X = radius_to_volume.(R)
    ξ = [2,1]
    two_drops = droplet_attributes{FT}(ξ.+0,R.+0,X.+0)
    two_drops_coag_data = coagulation_run{FT}(two_drops_coagsettings.Ns)
    two_drops_coag_data.ϕ[1] = 0.5
    two_drops_coag_data.pαdt[1] = 2
    @test two_drops_coag_data.lowest_zero[] == false
    test_pairs!(Serial(),2,[(1,2)],two_drops,two_drops_coag_data)
    @test two_drops_coag_data.lowest_zero[] == true
    @test two_drops.ξ == [0,1]
    @test two_drops.X[1] == two_drops.X[2] == 2*X[1]+X[2]
    @test two_drops.R == volume_to_radius.(two_drops.X)

    # test_pairs ξ_j_minus_γ_tilde_ξ_k > 0
    two_drops.ξ .= [6,2]
    two_drops.R .= [FT(1.5e-6), FT(1.5e-6)]
    two_drops.X .= X
    two_drops_coag_data.ϕ[1] = 0.5
    two_drops_coag_data.pαdt[1] = 2
    test_pairs!(Serial(),2,[(1,2)],two_drops,two_drops_coag_data)
    @test two_drops.ξ == [2,2]
    @test two_drops.X[2] == 2*X[1]+X[2]
    @test two_drops.X[1] == X[1]
    @test two_drops.R ≈ volume_to_radius.(two_drops.X)

    # test deficit
    begin
        # arrange
        two_drops.ξ .= [8,4]
        two_drops.R .= [FT(1.5e-6), FT(1.5e-6)]
        two_drops.X .= X
        two_drops_coag_data.ϕ[1] = 0.5
        two_drops_coag_data.pαdt[1] = 3

        # act
        test_pairs!(Serial(),2,[(1,2)],two_drops,two_drops_coag_data)

        # assert
        @test two_drops_coag_data.deficit[] == 4
    end

    two_drops_coag_data.deficit[] = 0

    #test split_highest_multiplicity
    ξ = [0,1,2,4]
    R = [FT(1.5e-6), FT(1.5e-6), FT(1.5e-6), FT(1.5e-6)]
    X = radius_to_volume.(R)
    four_drops = droplet_attributes{FT}(ξ.+0,R.+0,X.+0)
    split_highest_multiplicity!(four_drops)
    @test four_drops.ξ == [2,1,2,2]
    @test four_drops.X[1] == four_drops.X[4] == X[4]


    # pair_Ps
    L = [(1,2),(3,4)]
    ξ = [2,3,4,5]
    R = [1.5e-6,1.5e-6,1.5e-6,1.5e-6]
    X = radius_to_volume.(R)
    drops = droplet_attributes{FT}(ξ.+0,R.+0,X.+0)
    coag_data = coagulation_run{FT}(4)
    coagsettings = coag_settings{FT}(Ns = 4)
    compute_pαdt!(L,drops,coag_data,golovin,coagsettings)
    @test coag_data.pαdt[1] == 3*coagsettings.golovin_kernel_coeff*(X[1]+X[2])*(coagsettings.scale * coagsettings.Δt / coagsettings.ΔV)
    @test coag_data.pαdt[2] == 5*coagsettings.golovin_kernel_coeff*(X[3]+X[4])*(coagsettings.scale * coagsettings.Δt / coagsettings.ΔV)

    ξ = [2e9,3e10,4e9,5e10,(7e13-1),7e13]
    R = [3e-6,3e-6,3e-6,3e-6,3e-5,3e-5]
    X = radius_to_volume.(R)
    drops = droplet_attributes{FT}(ξ.+0,R.+0,X.+0)
    coag_data = coagulation_run{FT}(6)
    coagsettings = coag_settings{FT}(Ns = 6)
    coag_data.ϕ .= [0.5,1e-9,0.2]
    t_left = Ref(100.0)
    t_max = Vector{FT}(undef, 3)

    pair_Ps_adaptive!(1,(1,2), drops, coag_data,t_max,t_left,golovin,coagsettings)
    pair_Ps_adaptive!(2,(3,4), drops, coag_data,t_max,t_left,golovin,coagsettings)
    pair_Ps_adaptive!(3,(5,6), drops, coag_data,t_max,t_left,golovin,coagsettings)
    tmp_pair2 = ξ[4] * golovin(drops, (3,4), coagsettings)*(coagsettings.scale / coagsettings.ΔV)
    tmp_pair3 = ξ[6] * golovin(drops, (5,6), coagsettings)*(coagsettings.scale / coagsettings.ΔV)
    @test t_max[1] == t_left[]
    @test t_max[2] ≈ (div(5e10, 4e9))/tmp_pair2
    @test t_max[3] ≈ (div(ξ[6], ξ[5]))/tmp_pair3

    adaptive_pαdt!([(1,2),(3,4),(5,6)],drops,coag_data,t_left,golovin,coagsettings)
    @test coag_data.pαdt[2] == tmp_pair2 * (div(ξ[6], ξ[5]))/tmp_pair3
    @test coag_data.pαdt[3] == tmp_pair3 * (div(ξ[6], ξ[5]))/tmp_pair3
    @test t_left[] == 100 - (div(ξ[6], ξ[5]))/tmp_pair3
    

end

